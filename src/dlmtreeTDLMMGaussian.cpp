/**
 * @file dlmtreeTDLMMGaussian.cpp
 * @author Seongwon Im (seongwonim.github.io)
 * @brief Heterogeneous distributed lag mixture model 
 * @version 1.0
 * @date 2022-10-20
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "RcppEigen.h"
#include "Node.h"
#include "NodeStruct.h"
#include "modDat.h"
#include "exposureDat.h"
#include "Fncs.h"
#include "modelCtr.h"
using namespace Rcpp;

// MCMC updated
// 1. Tree pair: Tree1 & Tree 2
// 2. Exp object is a vector due to mixture setting
void dlmtreeTDLMMGaussian_TreeMCMC(int t, 
                                   Node* modTree, 
                                   Node* dlmTree1, Node* dlmTree2,
                                   dlmtreeCtr* ctr, dlmtreeLog *dgn,
                                   modDat* Mod, std::vector<exposureDat*> Exp);

// MHR updated
// 1. Tree pair: Tree1 & Tree 2
// 2. Exposure-variance: m1Var, m2Var, mixVar                      
treeMHR dlmtreeTDLMM_MHR(std::vector<Node*> modTerm,
                         std::vector<Node*> dlmTerm1, 
                         std::vector<Node*> dlmTerm2,
                         dlmtreeCtr* ctr, Eigen::VectorXd ZtR, 
                         double treeVar, double m1Var, double m2Var, double mixVar);

// [[Rcpp::export]]
Rcpp::List dlmtreeTDLMMGaussian(const Rcpp::List model){ 
  // Rcout << "*** Initializing dlmtreeTDLMMGaussian ***\n";
  // **** Set up general control variables ****
  // Model control(ctr) object with pointer
  dlmtreeCtr *ctr = new dlmtreeCtr;

  // MCMC parameters
  // Rcout << "MCMC parameters \n";
  ctr->iter = as<int>(model["nIter"]);        // Number of iterations
  ctr->burn = as<int>(model["nBurn"]);        // Number of burn-in
  ctr->thin = as<int>(model["nThin"]);        // Number of thinning 
  ctr->nRec = floor(ctr->iter / ctr->thin);   // Number of actual recording after thinning & burn-in
  ctr->nTrees = as<int>(model["nTrees"]);     // Number of ensemble tree

  // Data setup
  // Rcout << "Data setup \n";
  ctr->Y = as<Eigen::VectorXd>(model["Y"]);   // Response variable
  ctr->n = (ctr->Y).size();                   // Sample size

  // Modifier tree hyperparameter
  ctr->modZeta = as<double>(model["zeta"]);   // Hyperparameter for modifier tree
  ctr->modKappa = 100;                        // Hyperparameter for Dirichlet-Categorial 

  // Model estimation initialization
  // Rcout << "Model pre-calculation \n";
  ctr->Z = as<Eigen::MatrixXd>(model["Z"]);         // Fixed effect model matrix
  ctr->Zw = ctr->Z;                                 // Zw in case of logistic/zinb (w = omega)
  ctr->pZ = (ctr->Z).cols();                        // Number of fixed effect covariates

  ctr->VgInv = (ctr->Z).transpose() * (ctr->Z);     // V_gamma initial calculation
  ctr->VgInv.diagonal().array() += 1.0 / 100000.0;  // V_gamma + prior
  ctr->Vg = ctr->VgInv.inverse();                   // V_gamma inverse
  ctr->VgChol = (ctr->Vg).llt().matrixL();          // V_gamma inverse cholesky
  
  // Binomial flag set as false                     
  ctr->binomial = 0;                                // Binomial flag (Never used but may be later)

  // Diagnostics & messages
  ctr->verbose = bool(model["verbose"]);            // Technical messages
  ctr->diagnostics = bool(model["diagnostics"]);   // Technical diagnostics

  // Rcout << "Step probabilities \n";
  // Tree prior probabilities for the modifier tree and DLM tree
  ctr->stepProbMod = as<std::vector<double>>(model["stepProbMod"]);   // modifier tree step prob.
  ctr->stepProb = as<std::vector<double>>(model["stepProbTDLM"]);     // dlm tree step prob.
  ctr->treePriorMod = as<std::vector<double>>(model["treePriorMod"]); // modifier tree alpha, beta
  ctr->treePrior = as<std::vector<double>>(model["treePriorTDLM"]);   // dlm tree alpha, beta

  // [Mixture & Shrinkage] 
  // Store shrinkage: 
  // 3 = all (mu(s), tau(t))
  // 2 = trees (tau(t))
  // 1 = exposures/interactions (mu(s)) - (default)
  // 0 = none
  ctr->shrinkage = as<int>(model["shrinkage"]);

  // Different notation: modKappa in TDLM = mixKappa in HDLM
  // Rcout << "Mixture: mixKappa setup \n";
  ctr->mixKappa = as<double>(model["mixPrior"]);  // Positive scalar hyperparameter for sparsity of exposures (Linero): Set as 1
  bool updateKappa = false;                       // Check to update mixKappa: if negative, set to 1 as default.
  if (ctr->mixKappa < 0) {                        // This "if" condition will not be used as default is 1 
    updateKappa = true;
    ctr->mixKappa = 1;
  }

  // **** Setup modifier data ****
  // Rcout << "modifier data/parameter \n";
  modDat *Mod = new modDat(as<std::vector<int>>(model["modIsNum"]),   // Boolean vector of continuous vs categorical
                           as<Rcpp::List>(model["modSplitIdx"]),      // If we have wanted split points, read row indice in  
                           as<std::vector<int>>(model["fullIdx"]));   // 0:(nrow(data) - 1)

  NodeStruct *modNS;            // NodeStruct object for the modifier tree
  modNS = new ModStruct(Mod);   // Take Mod info and construct the modifier tree
  ctr->pM = Mod->nMods;         // Number of modifiers (Simulation has three)

  // **** Pre-calculate single node tree pair matrices ****
  // Rcout << "Multiple exposure parameter control \n";
  std::vector<exposureDat*> Exp;                      // Exposure structure
  Rcpp::List exp_dat = as<Rcpp::List>(model["X"]);    // List of exposures
  ctr->nExp = exp_dat.size();                         // Number of exposures
  for (int i = 0; i < ctr->nExp; i++) {
    Exp.push_back(
      new exposureDat(
        as<Eigen::MatrixXd>(
          as<Rcpp::List>(exp_dat[i])["Tcalc"]), ctr->Z, ctr->Vg));
  }

  // Rcout << "Exposure lag & interaction \n";
  ctr->pX = Exp[0]->pX;     // Total number of lags (columns of exposure matrix X)
  ctr->nSplits = 0;         // Never used

  // **** Mixture / interaction management ****
  ctr->interaction = as<int>(model["interaction"]);   // Interaction: TDLMMadd(1) / TDLMMns(2) / TDLMMall(3)
  ctr->nMix = 0;                                      // Interaction flag with a default of none
  if (ctr->interaction) { 
    ctr->nMix += int (ctr->nExp * (ctr->nExp - 1.0) / 2.0); // (M choose 2) for HDLMMns
    if (ctr->interaction == 2) { // HDLMMall
      ctr->nMix += ctr->nExp;                               // (M choose 2) + M for HDLMMall
    }
  }

  // **** Create trees ****
  // Rcout << "Creating trees \n";
  NodeStruct *expNS;   // NodeStruct object
  // DLNMStruct takes (xmin_in, xmax_in, 
  //                   tmin_in, tmax_in, 
  //                   Xp_in, Tp_in)
  expNS = new DLNMStruct(0, ctr->nSplits + 1,  // xmin_in, xmax_in 
                         1, int (ctr->pX),     // tmin_in, tmax_in
                         as<Eigen::VectorXd>(model["splitProb"]), // Xprob_in
                         as<Eigen::VectorXd>(model["timeProb"])); // Tprob_in

  // Exposure information vectors
  ctr->expProb = as<Eigen::VectorXd>(model["expProb"]);  // Exposure selection probability

  (ctr->dlmTree1Exp).resize(ctr->nTrees);    
  (ctr->dlmTree2Exp).resize(ctr->nTrees);

  // Fixed exposures for hdlmm
  ctr->dlmTree1Exp = as<Eigen::VectorXd>(model["comb1"]);
  ctr->dlmTree2Exp = as<Eigen::VectorXd>(model["comb2"]);

  (ctr->expCount).resize((ctr->expProb).size()); // Count how many are currently selected
  (ctr->expInf).resize((ctr->expProb).size());   // Sum of shrinkage prior for math later

  // [Tree pair] ===================================================
  // Rcout << "Creating tree pairs \n";
  int t; // Tree index
  
  std::vector<Node*> modTrees;
  std::vector<Node*> dlmTrees1;    // Tree1 Ensemble vector
  std::vector<Node*> dlmTrees2;    // Tree2 Ensemble vector

  // For loop iterating through trees to start the "roots"
  for (t = 0; t < ctr->nTrees; t++) {
    // Modifier tree initialize
    modTrees.push_back(new Node(0, 1));         // Root node initialize
    modTrees[t]->nodestruct = modNS->clone();   // Each modifier tree shares the same property as modNS
    Mod->updateNodeVals(modTrees[t]);           // update modifier root node

    // // dlmtree pair initialize -> No swap step for HDLMM, fix the exposures
    // ctr->dlmTree1Exp(t) = (ctr->comb1)[t]; // Each dlmtree1 is randomly assigned an exposure
    // ctr->dlmTree2Exp(t) = (ctr->comb2)[t]; // Each dlmtree2 is randomly assigned an exposure

    dlmTrees1.push_back(new Node(0, 1));        // dlmtree1 Ensemble Root node initialize
    dlmTrees2.push_back(new Node(0, 1));        // dlmtree2 Ensemble Root node initialize
    dlmTrees1[t]->nodestruct = expNS->clone();  // Each dlmtree1 shares the same property as expNS
    dlmTrees2[t]->nodestruct = expNS->clone();  // Each dlmtree2 shares the same property as expNS
    Exp[ctr->dlmTree1Exp(t)]->updateNodeVals(dlmTrees1[t]); // Update dlmtree1 root node
    Exp[ctr->dlmTree2Exp(t)]->updateNodeVals(dlmTrees2[t]); // Update dlmtree2 root node
  }
  // ================================================================
  delete modNS; // delete the object used for cloning a modifier tree
  delete expNS; // delete the object used for cloning a exposure dlmtree

  // **** Logs ****
  // Rcout << "Logs \n";
  dlmtreeLog *dgn = new dlmtreeLog;   // Initialize log for the output
  // dlmTree log with ctr->nRec
  // Rcout << "dlmTree log \n";
  (dgn->gamma).resize(ctr->pZ, ctr->nRec);    (dgn->gamma).setZero();   // fixed effect coefficient
  (dgn->sigma2).resize(ctr->nRec);            (dgn->sigma2).setZero();  // sigma2 variance (single value)
  (dgn->nu).resize(ctr->nRec);                (dgn->nu).setZero();      // Global shrinkage (single value)
  (dgn->tau).resize(ctr->nTrees, ctr->nRec);  (dgn->tau).setZero();     // Tree shrinkage (vector)
  (dgn->fhat).resize(ctr->n);                 (dgn->fhat).setZero();    // partial fitted
  (dgn->totTerm).resize(ctr->nRec);           (dgn->totTerm).setZero(); // P (total term #)

  // Exposure-specific shrinkage - SI 
  // Rcout << "Exposure-specific shrinkage \n";
  (dgn->muExp).resize(ctr->nExp, ctr->nRec);        (dgn->muExp).setZero();   // Exposure-specific variance
  if (ctr->interaction > 0) { // For TDLMMns / TDLMMall
    // A vector of all possible combinations of exposures
    (dgn->muMix).resize(ctr->nMix, ctr->nRec);      (dgn->muMix).setZero();   
    (dgn->mixInf).resize(ctr->nMix, ctr->nRec);     (dgn->mixInf).setZero(); 
    (dgn->mixCount).resize(ctr->nMix, ctr->nRec);   (dgn->mixCount).setZero();
  } else { // For TDLMMadd -> No interaction terms -> nMix is 0 so just set them as constants
    (dgn->muMix).resize(1, 1);                      (dgn->muMix).setZero();
    (dgn->mixInf).resize(1, 1);                     (dgn->mixInf).setZero(); 
    (dgn->mixCount).resize(1, 1);                   (dgn->mixCount).setZero();
  }

  // Modifier tree log
  // Rcout << "Modtree log \n";
  (dgn->modProb).resize(ctr->pM, ctr->nRec);  (dgn->modProb).setZero();   // modifier probability
  (dgn->modCount).resize(ctr->pM, ctr->nRec); (dgn->modCount).setZero();  // how many current modifiers
  (dgn->modInf).resize(ctr->pM, ctr->nRec);   (dgn->modInf).setZero();    // Sum of shrinkage prior for math
  (dgn->modKappa).resize(ctr->nRec);          (dgn->modKappa).setZero();  // (sparsity for categorical)
  (dgn->termNodesMod).resize(ctr->nTrees, ctr->nRec);   (dgn->termNodesMod).setZero();  // Number of terminal node for each modifier tree

  // Exposure log
  // Rcout << "Exposure information log \n";
  (dgn->expProb).resize(ctr->nExp, ctr->nRec);          (dgn->expProb).setZero();       // Exposure selection probability
  (dgn->expCount).resize(ctr->nExp, ctr->nRec);         (dgn->expCount).setZero();      // Number of current selected exposures
  (dgn->expInf).resize(ctr->nExp, ctr->nRec);           (dgn->expInf).setZero();        // Sum of prior shrinkage for math
  (dgn->dlmTree1Exp).resize(ctr->nTrees, ctr->nRec);    (dgn->dlmTree1Exp).setZero();      // Exposure vector of tree vector1
  (dgn->dlmTree2Exp).resize(ctr->nTrees, ctr->nRec);    (dgn->dlmTree2Exp).setZero();      // Exposure vector of tree vector2
  (dgn->termNodesDLM1).resize(ctr->nTrees, ctr->nRec);  (dgn->termNodesDLM1).setZero(); // Number of terminal node for each dlmtree1
  (dgn->termNodesDLM2).resize(ctr->nTrees, ctr->nRec);  (dgn->termNodesDLM2).setZero(); // Number of terminal node for each dlmtree2
  (dgn->mixKappa).resize(ctr->nRec);                    (dgn->mixKappa).setZero();      // This is essentially a vector of one

  // **** Initial draws ****
  // Rcout << "Partial residual setup \n";
  (ctr->fhat).resize(ctr->n);      (ctr->fhat).setZero();    // Partial residual set to zero
  ctr->R = ctr->Y;                  // Partial residual is initially Y as fhat is 0
  (ctr->gamma).resize(ctr->pZ);     // Coefficient for each fixed effect covariate

  // **** Exposure specific parameters ****
  // Rcout << "Exposure terminal node counting \n";
  (ctr->totTermExp).resize(ctr->nExp);      (ctr->totTermExp).setZero();   // total number of terminal nodes per exposure
  (ctr->sumTermT2Exp).resize(ctr->nExp);    (ctr->sumTermT2Exp).setZero(); // (delta_a)^2 per exposure for MHR math
  (ctr->muExp).resize(ctr->nExp);           (ctr->muExp).setOnes();        // Exposure-specific variance

  // If there is an interaction, define more parameters
  // Rcout << "Updating interaction parameters \n";
  // Initiate interaction effect muMix as a matrix of 1 -> start mixVar = 1
  if (ctr->interaction) { 
    (ctr->totTermMix).resize(ctr->nExp, ctr->nExp);   (ctr->totTermMix).setZero();    // total number of terminal nodes for exposure combinations
    (ctr->sumTermT2Mix).resize(ctr->nExp, ctr->nExp); (ctr->sumTermT2Mix).setZero();  // (delta_a)^2 per exposure combination for MHR math
    (ctr->muMix).resize(ctr->nExp, ctr->nExp);        (ctr->muMix).setOnes();         // interaction-specific variance
    (ctr->mixInf).resize(ctr->nExp, ctr->nExp);       (ctr->mixInf).setZero();        // interaction shrinkage prior
    (ctr->mixCount).resize(ctr->nExp, ctr->nExp);     (ctr->mixCount).setZero();      // Count of how many mixtures there are
  }

  // No MCMC run yet so zero terminal nodes
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  ctr->nu = 1.0;
  ctr->sigma2 = 1.0;

  // Rcout << "Initial model fitting \n";
  tdlmModelEst(ctr); // Initial Gaussian sampling for gamma / fhat / R

  // Shrinkage hyperparameters
  rHalfCauchyFC(&(ctr->nu), ctr->nTrees, 0.0);            // Global shrinkage: nu
  (ctr->tau).resize(ctr->nTrees);  (ctr->tau).setOnes();  // Tree shrinkage: tau(t)
  if (ctr->shrinkage > 1){
    for (t = 0; t < ctr->nTrees; t++) {
      rHalfCauchyFC(&((ctr->tau)(t)), 0.0, 0.0);
    }
  }

  // **** Initializing trees ****
  // DLM tree pairs
  // Rcout << "dlmTree pair terminal nodes # \n";
  ctr->nTermDLM.resize(ctr->nTrees);     (ctr->nTermDLM).array() = 2;    // Every tree pair starts with 2 (root) terminal nodes
  ctr->nTermDLM1.resize(ctr->nTrees);    (ctr->nTermDLM1).array() = 1;   // Every tree 1 starts with 1 (root) terminal node
  ctr->nTermDLM2.resize(ctr->nTrees);    (ctr->nTermDLM2).array() = 1;   // Every tree 2 starts with 1 (root) terminal node

  // Modifier trees
  // Rcout << "modTree pair terminal nodes # \n";
  ctr->nTermMod.resize(ctr->nTrees);   (ctr->nTermMod).array() = 1; // Every modifier tree starts with 1 (root) terminal node
  ctr->modCount.resize(ctr->pM);       (ctr->modCount).setZero();   // Modifier covariate count vector
  ctr->modInf.resize(ctr->pM);         (ctr->modInf).setZero();     // Sum starts from zero

  // Partial residual matrix
  (ctr->Rmat).resize(ctr->n, ctr->nTrees);  (ctr->Rmat).setZero();  // Each ensemble has a partial residual
  
  // **** MCMC ****
  // Rcout << "Initializing MCMC iteration \n";
  // Create a progress meter
  progressMeter* prog = new progressMeter(ctr);

  // Thinning and burn-in process
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }

    // **** Update trees and parameters ****
    // Rcout << "Update trees \n";
    ctr->R += (ctr->Rmat).col(0); 

    // Reset DLM tree parameters
    (ctr->fhat).setZero();
    ctr->totTerm = 0.0;               
    ctr->sumTermT2 = 0.0;

    // Reset exposure parameters
    (ctr->totTermExp).setZero();
    (ctr->sumTermT2Exp).setZero();
    (ctr->expCount).setZero();
    (ctr->expInf).setZero();

    // Reset interaction parameters
    (ctr->mixCount).setZero();
    (ctr->mixInf).setZero();

    if (ctr->interaction > 0) {
      (ctr->totTermMix).setZero(); 
      (ctr->sumTermT2Mix).setZero();
    }

    // Reset modifier tree parameters
    (ctr->modCount).setZero();
    (ctr->modInf).setZero();

    // For each dlmTree pair & Modifier tree, perform one MCMC iteration
    // Rcout << "Iterating through the tree ensemble \n";
    for (t = 0; t < ctr->nTrees; t++) {
      dlmtreeTDLMMGaussian_TreeMCMC(t,                    // Tree index
                                    modTrees[t],          // Modtree
                                    dlmTrees1[t],         // Tree 1 of t-th pair
                                    dlmTrees2[t],         // Tree 2 of t-th pair
                                    ctr, dgn, Mod, Exp);  // ctr, record, mod, exposure
      
      //[Bayesian backfitting] Update the fitted values after each iteration of t 
      ctr->fhat += (ctr->Rmat).col(t);

      // For each tree pair, update the partial residual
      // Note: ctr->R is used as t is iterated
      if (t < ctr->nTrees - 1)
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
    }

    // Rcout << "Tree ensemble update finished \n";
    // Update the parameters
    ctr->R = ctr->Y - ctr->fhat;
    ctr->sumTermT2 = (ctr->sumTermT2Exp).sum();
    ctr->totTerm = (ctr->totTermExp).sum();
    if(ctr->interaction > 0) {
      ctr->sumTermT2 += (ctr->sumTermT2Mix).sum();
      ctr->totTerm += (ctr->totTermMix).sum();
    }

    // Update gamma
    // Rcout << "Updating Gaussian(gamma) \n";
    tdlmModelEst(ctr);

    // Exposure Shrinkage update
    // Shrinkage
    // 2 = trees (tau(t))
    // 1 = exposures/interactions (mu(s)) - (default)
    // 0 = none
    // Rcout << "Updating the hyperparameters for Gaussian model \n";
    // Sample xi ~ Inv-Gamma(1, 1/(1 + 1/nu)) = xiInv ~ Gamma(1, 1/(1 + 1/nu))
    // Update nu full conditional

    // Global shrinkage: nu
    double xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->nu)));
    ctr->nu = 1.0 / R::rgamma(0.5 * ctr->totTerm + 0.5,
                              1.0 / (0.5 * ctr->sumTermT2 / (ctr->sigma2) + xiInv));
                              
    // Exposure shrinkage
    double sigmanu = ctr->sigma2 * ctr->nu;
    if (ctr->shrinkage == 1 || ctr->shrinkage == 3) { 
      for (int i = 0; i < ctr->nExp; i++) {
        // Eq(16) & Eq(17) from TDLMM
        xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->muExp(i))));
        ctr->muExp(i) = 1.0 / R::rgamma(0.5 * ctr->totTermExp(i) + 0.5,
                        1.0 / (0.5 * ctr->sumTermT2Exp(i) / sigmanu + xiInv));

        // Similarly update interaction-specific terms as well
        if (ctr->interaction) {
          for (int j = i; j < ctr->nExp; j++) {
            if ((j > i) || (ctr->interaction == 2)){
              rHalfCauchyFC(&(ctr->muMix(j, i)), ctr->totTermMix(j, i),
                            ctr->sumTermT2Mix(j, i) / sigmanu);
              // Rcout << "muMix: e";
              // Rcout << i;
              // Rcout << "-e";
              // Rcout << j;
              // Rcout << ": ";
              // Rcout << ctr->muMix(j, i);
              // Rcout << "\n";
            }
          } // end for loop updating interaction variances
        } // end if interactions
      } // end for loop updating exposure variances
    } // end if shrinkage == 1 or 3

    // **** Update modifier selection ****
    // Rcout << "Begin Modtree selection update after burn-in \n";
    if ((ctr->b > 1000) || (ctr->b > (0.5 * ctr->burn))) {
      // Modifier selection start
      double beta = R::rbeta(ctr->modZeta, 1.0);
      double modKappaNew = beta * ctr->pM / (1 - beta);
      // Dirichlet MHR update
      double mhrDir = logDirichletDensity(Mod->modProb, (ctr->modCount.array() + modKappaNew / ctr->pM).matrix()) -
                      logDirichletDensity(Mod->modProb, (ctr->modCount.array() + ctr->modKappa / ctr->pM).matrix());
      if (log(R::runif(0, 1) < mhrDir) && (mhrDir == mhrDir)) {
        ctr->modKappa = modKappaNew;
      }

      Mod->modProb = rDirichlet((ctr->modCount.array() + ctr->modKappa / ctr->pM).matrix());
      // end modifier selection

      // HDLMM: Exposure selection posterior calculation
      if (updateKappa) {
        double mixKappaNew = R::rgamma(1.0, ctr->nTrees/4.0);
        double mhrDirMix =
          logDirichletDensity(ctr->expProb, ((ctr->expCount).array() + mixKappaNew).matrix()) -
          logDirichletDensity(ctr->expProb, ((ctr->expCount).array() + ctr->mixKappa).matrix());

        if (log(R::runif(0, 1)) < mhrDirMix && (mhrDirMix == mhrDirMix)){
          ctr->mixKappa = mixKappaNew;
        }
      }
      
      // Update the exposure selection probability
      ctr->expProb = rDirichlet(((ctr->expCount).array() + ctr->mixKappa).matrix());
    } // HDLMM exposure selection end

    // **** Record every iteration ****
    // Rcout << "Record: Updating dgn \n";
    if (ctr->record > 0) {
      // Tree pair parameters record (prior & exposures)
      (dgn->gamma).col(ctr->record - 1) = ctr->gamma;
      (dgn->sigma2)(ctr->record - 1) = ctr->sigma2;
      (dgn->nu)(ctr->record - 1) = ctr->nu;
      (dgn->totTerm)(ctr->record - 1) = ctr->totTerm;
      (dgn->tau).col(ctr->record - 1) = ctr->tau;

      (dgn->termNodesDLM1).col(ctr->record - 1) = ctr->nTermDLM1;
      (dgn->termNodesDLM2).col(ctr->record - 1) = ctr->nTermDLM2;
      (dgn->dlmTree1Exp).col(ctr->record - 1) = ctr->dlmTree1Exp;
      (dgn->dlmTree2Exp).col(ctr->record - 1) = ctr->dlmTree2Exp;
      (dgn->expProb).col(ctr->record - 1) = ctr->expProb;
      (dgn->expCount).col(ctr->record - 1) = ctr->expCount;
      (dgn->expInf).col(ctr->record - 1) = ctr->expInf;

      // Mixtures record
      (dgn->muExp).col(ctr->record - 1) = ctr->muExp;
      (dgn->mixKappa)(ctr->record - 1) = ctr->mixKappa;

      // Interaction record
      if (ctr->interaction > 0) {
        int k = 0;
        for (int i = 0; i < ctr->nExp; i++) {
          for (int j = i; j < ctr->nExp; j++) {
            if ((j > i) || (ctr->interaction == 2)) {
              dgn->muMix(k, ctr->record - 1) = ctr->muMix(j, i);
              dgn->mixInf(k, ctr->record - 1) = ctr->mixInf(j, i);
              dgn->mixCount(k, ctr->record - 1) = ctr->mixCount(j, i);
              k++;
            }
          }
        }
      }

      // Modifier tree record
      (dgn->termNodesMod).col(ctr->record - 1) = ctr->nTermMod;
      (dgn->modProb).col(ctr->record - 1) = Mod->modProb;
      (dgn->modCount).col(ctr->record - 1) = ctr->modCount;
      (dgn->modInf).col(ctr->record - 1) = ctr->modInf / ctr->modInf.maxCoeff();
      (dgn->modKappa)(ctr->record - 1) = ctr->modKappa;

      dgn->fhat += ctr->fhat;

    } // end record

    // if (ctr->diagnostics){
    //   Rcout << "MCMC Iteration : ";
    //   Rcout << ctr->b;
    //   Rcout << "\n";

    //   Rcout << "Sigma2 : ";
    //   Rcout << ctr->sigma2;
    //   Rcout << "\n";
    //   Rcout << "nu : ";
    //   Rcout << ctr->nu;
    //   Rcout << "\n";

    //   Rcout << "tau : ";
    //   for(int p = 0; p < ctr->nTrees; p++){
    //     Rcout << ctr->tau[p];
    //     Rcout << "|";
    //   }
    //   Rcout << "\n";

    //   Rcout << "muExp(0) : ";
    //   Rcout << ctr->muExp(0);
    //   Rcout << "\n";
    //   Rcout << "muExp(1) : ";
    //   Rcout << ctr->muExp(1);
    //   Rcout << "\n";
    //   Rcout << "muExp(2) : ";
    //   Rcout << ctr->muExp(2);
    //   Rcout << "\n";

    //   if(ctr->interaction > 0){
    //     Rcout << "muMix(0-1) : ";
    //     Rcout << ctr->muMix(1, 0);
    //     Rcout << "\n";
    //     Rcout << "muMix(1-2) : ";
    //     Rcout << ctr->muMix(2, 1);
    //     Rcout << "\n";
    //     Rcout << "muMix(0-2) : ";
    //     Rcout << ctr->muMix(2, 0);
    //     Rcout << "\n";
    //   }
    
    //   Rcout << "totTerm (denoted P) : ";
    //   Rcout << ctr->totTerm;
    //   Rcout << "\n";

    //   Rcout << "sumTermT2 (denoted D) : ";
    //   Rcout << ctr->sumTermT2;
    //   Rcout << "\n";
    //   Rcout << "----------------- \n";
    // }
    
    // Progress mark
    prog->printMark();
  } // end MCMC
  // Rcout << "MCMC complete \n";


  // **** Prepare outout ****
  // Rcout << "Preparing output \n";
  Eigen::MatrixXd TreeStructs((dgn->DLMexp).size(), 10);    // TreeStructs is the same as DLM
  Rcpp::StringVector termRule(dgn->termRule.size());        // Store rules for DLM
  Rcpp::StringVector termRuleMIX(dgn->termRuleMIX.size());  // Store rules for MIX

  // Store rec() vectors
  std::size_t s;                                // std::size_t is similar to (int)
  for (s = 0; s < (dgn->DLMexp).size(); s++){    // rbind all rec()'s to build TreeStructs
    TreeStructs.row(s) = dgn->DLMexp[s];
  }

  // Modifier rule per main & interaction effect
  termRule = dgn->termRule;
  termRuleMIX = dgn->termRuleMIX;

  // Transfer variables from the dgn object
  // DLM tree pair
  // Rcout << "DLM tree pair \n";
  Eigen::VectorXd sigma2 = dgn->sigma2;
  Eigen::VectorXd nu = dgn->nu;
  Eigen::MatrixXd tau = (dgn->tau).transpose();
  Eigen::VectorXd totTerm = dgn->totTerm;
  Eigen::VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  Eigen::MatrixXd gamma = (dgn->gamma).transpose();

  // dlmTree information
  // Rcout << "Exposures \n";
  Eigen::MatrixXd termNodesDLM1 = (dgn->termNodesDLM1).transpose();
  Eigen::MatrixXd termNodesDLM2 = (dgn->termNodesDLM2).transpose();
  Eigen::MatrixXd tree1Exp = (dgn->dlmTree1Exp).transpose();
  Eigen::MatrixXd tree2Exp = (dgn->dlmTree2Exp).transpose();
  Eigen::MatrixXd expProb = (dgn->expProb).transpose();
  Eigen::MatrixXd expCount = (dgn->expCount).transpose();
  Eigen::MatrixXd expInf = (dgn->expInf).transpose();
  Eigen::MatrixXd muExp = (dgn->muExp).transpose();

  // Mixture
  Eigen::MatrixXd mixInf = (dgn->mixInf).transpose();
  Eigen::MatrixXd mixCount = (dgn->mixCount).transpose();
  Eigen::MatrixXd muMix(1, 1);  muMix.setZero();
  Eigen::MatrixXd MIX(0, 10);   MIX.setZero();
  Eigen::VectorXd mixKappa = dgn->mixKappa;

  // Rcout << "Storing interaction MCMC \n";
  // If interaction, expand muMIX and MIX
  if (ctr->interaction) {
    muMix.resize((dgn->muMix).cols(), (dgn->muMix).rows());
    muMix = (dgn->muMix).transpose();
    MIX.resize((dgn->MIXexp).size(), 10);
    for (s = 0; s < (dgn->MIXexp).size(); ++s)
      MIX.row(s) = dgn->MIXexp[s];
  }

  // Modifier tree
  // Rcout << "Modifier tree record \n";
  Eigen::MatrixXd termNodesMod = (dgn->termNodesMod).transpose();
  Eigen::VectorXd kappa = dgn->modKappa;
  Eigen::MatrixXd modProb = (dgn->modProb).transpose();
  Eigen::MatrixXd modCount = (dgn->modCount).transpose();
  Eigen::MatrixXd modInf = (dgn->modInf).transpose();

  Eigen::MatrixXd modAccept((dgn->treeModAccept).size(), 9);
  Eigen::MatrixXd dlmAccept((dgn->treeDLMAccept).size(), 9);
  for (s = 0; s < (dgn->treeDLMAccept).size(); s++){
    dlmAccept.row(s) = dgn->treeDLMAccept[s];
  }

  delete prog;
  delete ctr;
  delete dgn;
  for (s = 0; s < Exp.size(); s++){        // Exposure deletion
    delete Exp[s];
  }
  delete Mod;
  for (s = 0; s < modTrees.size(); s++) { // Trees deletion
    delete modTrees[s];
    delete dlmTrees1[s];
    delete dlmTrees2[s];
  }

  return(Rcpp::List::create(Named("TreeStructs") = wrap(TreeStructs), // This is DLMexp
                            Named("MIX") = wrap(MIX),
                            Named("termRules") = wrap(termRule),
                            Named("termRuleMIX") = wrap(termRuleMIX),
                            //Named("fhat") = wrap(fhat),
                            Named("sigma2") = wrap(sigma2),
                            Named("nu") = wrap(nu),
                            Named("tau") = wrap(tau),
                            Named("totTerm") = wrap(totTerm),
                            Named("gamma") = wrap(gamma),
                            Named("termNodesMod") = wrap(termNodesMod),
                            Named("termNodesDLM1") = wrap(termNodesDLM1),
                            Named("termNodesDLM2") = wrap(termNodesDLM2),
                            Named("mixKappa") = wrap(mixKappa),
                            Named("modProb") = wrap(modProb),
                            //Named("expInf") = wrap(expInf),
                            Named("expCount") = wrap(expCount),
                            //Named("mixInf") = wrap(mixInf),
                            //Named("mixCount") = wrap(mixCount),
                            Named("muExp") = wrap(muExp),
                            //Named("muMix") = wrap(muMix),
                            Named("modCount") = wrap(modCount),
                            Named("modInf") = wrap(modInf),
                            //Named("treeModAccept") = wrap(modAccept),
                            Named("treeDLMAccept") = wrap(dlmAccept)));

} // end dlmtreeTDLMMGaussian



void dlmtreeTDLMMGaussian_TreeMCMC(int t, Node* modTree, 
                                   Node* dlmTree1, Node* dlmTree2,
                                   dlmtreeCtr* ctr, dlmtreeLog *dgn,
                                   modDat* Mod, std::vector<exposureDat*> Exp)
{ 
  // Rcout << "Reminder: The code is in the tree pair for each t";
  // Rcout << "TreeMCMC: Defining parameters... \n";
  int m1, m2, newExp;                                 // Exposure (& potentially new exposure)
  int step, step1, step2;                             // MCMC transition steps
  double m1Var, m2Var, mixVar, newExpVar, newMixVar;  // Exposure specific variance
  int success = 0;                                    // MH ratio parameter
  double RtR = -1.0;                                  // MH ratio parameter
  double RtZVgZtR = 0;                                // MH ratio parameter
  double stepMhr = 0;                                 // MH ratio parameter
  double ratio = 0;                                   // MH ratio parameter
  double treeVar = (ctr->nu) * (ctr->tau)(t);         // Shrinkage: Global x Local
  std::size_t s;                                      // modifier tree count
  std::vector<Node*> modTerm, dlmTerm1, dlmTerm2;             // Terminal nodes
  std::vector<Node*> newModTerm, newDlmTerm1, newDlmTerm2;    // Proposed terminal nodes
  Node* newTree = 0;                                  // Potential new tree
  treeMHR mhr0, mhr;                                  // mhr object

  // Pre-calculation for MH ratio update
  // Rcout << "TreeMCMC: Precalculation for MHR update... \n";
  Eigen::VectorXd ZtR = (ctr->Z).transpose() * (ctr->R);

  // -- List terminal nodes --
  // Rcout << "TreeMCMC: Listing terminal nodes... \n";
  modTerm = modTree->listTerminal();    // modifier tree terminal nodes
  dlmTerm1 = dlmTree1->listTerminal();  // tree1 terminal nodes
  dlmTerm2 = dlmTree2->listTerminal();  // tree2 terminal nodes

  // Extract exposure information from ctr
  m1 = ctr->dlmTree1Exp[t];     // Exposure m1 of tree 1 of t th tree pair 
  m2 = ctr->dlmTree2Exp[t];     // Exposure m2 of tree 2 of t th tree pair
  m1Var = ctr->muExp(m1);    // Exposure-specific variance
  m2Var = ctr->muExp(m2);    // Exposure-specific variance
  mixVar = 0;                // No interaction variance for now
  // If there is an interaction, extract mixture variance value from the muMix matrix
  if ((ctr->interaction) && ((ctr->interaction == 2) || (m1 != m2))) { 
    if (m1 <= m2)
      mixVar = ctr->muMix(m2, m1);
    else
      mixVar = ctr->muMix(m1, m2);
  }

  // Current null state MHR
  // Rcout << "[mhr0] \n";
  mhr0 = dlmtreeTDLMM_MHR(modTerm, dlmTerm1, dlmTerm2, 
                          ctr, ZtR, treeVar, 
                          m1Var, m2Var, mixVar);
                        
  // [Tree pair update]
  // **** Propose a new TDLMM tree 1 ****
  // "New" parameters are set as the current ones until being updated by switch-exposure step.
  // Rcout << "TreeMCMC: Updating tree 1 \n";
  newExp    = m1; 
  newExpVar = m1Var;
  newMixVar = mixVar;
  stepMhr = 0;
  success = 0;

  // Choose a transition and update
  // Rcout << "TreeMCMC: Selecting transition \n";
  step1 = sampleInt(ctr->stepProb, 1);
  if ((dlmTerm1.size() == 1) && (step1 < 3)){
    step1 = 0;
  }

  // This part is different from HDLM as we have multiple exposures now
  if(step1 < 3){  // Grow / Prune / Change
    stepMhr = tdlmProposeTree(dlmTree1, Exp[m1], ctr, step1); // Take the proposed step and propose a tree
    success = dlmTree1->isProposed();                         // Evaluate the proposed tree
    newDlmTerm1 = dlmTree1->listTerminal(1);                  // Save the new terminal node
  } else { // Step 4 = Propose a new exposure
    newExp = sampleInt(ctr->expProb);
    if (newExp != m1) { // To update, the proposed exposure must be different to the current one.
      // Update the information with new exposure
      success   = 1;   // exposure-switching is always successful transition                   
      newExpVar = ctr->muExp(newExp);         // Extract the variance of the new exposure
      newTree   = new Node(*dlmTree1);        // Propose a tree with a new exposure
      newTree->setUpdate(1);                  // Update
      newDlmTerm1 = newTree->listTerminal();  // Save the terminal nodes of the tree with a new exposure

      for (Node* nt : newDlmTerm1)        // For each terminal node,
        Exp[newExp]->updateNodeVals(nt);  // Update with a new exposure

      // Update the interaction using the new exposure as well for TDLMMns/all
      if ((ctr->interaction) && ((ctr->interaction == 2) || (newExp != m2))) {
        if (newExp <= m2)
          newMixVar = ctr->muMix(m2, newExp); 
        else
          newMixVar = ctr->muMix(newExp, m2); 
      } else { // TDLMMadd -> no new MixVar
        newMixVar = 0;
      }
    }
  }

  // * dlmTree 1 MHR 
  // StepMHR: Transition ratio
  // mhr0 & mhr: Likelihood of R
  // Rcout << "TreeMCMC: Compute the MHR of the current tree/terminal nodes, called mhr0 \n";

  // If successful, calculate MH ratio
  // Rcout << "TreeMCMC: Success: -> Compute the MHR of the proposed tree/terminal nodes, called mhr \n";
  if (success) {  
    // If transition is successful, we have a new dlm terminal node.
    newDlmTerm1 = dlmTree1->listTerminal(1);
    modTree->setUpdateXmat(1);

    // MH ratio with a new terminal and the exposure: newDlmTerm1, newExpVar
    // Rcout << "[Tree1: mhr] \n";
    mhr = dlmtreeTDLMM_MHR(modTerm, 
                           newDlmTerm1, dlmTerm2, ctr, ZtR, treeVar, 
                           newExpVar, m2Var, newMixVar);

    // MH ratio - dlmTree 
    if (RtR < 0) {
      RtR = (ctr->R).dot(ctr->R);
      RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);
    }

    // MH ratio
    ratio = stepMhr + 
            mhr.logVThetaChol - mhr0.logVThetaChol -
            (0.5 * (ctr->n + 1.0) *
            (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
            log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
            (0.5 * ((mhr.nTerm1 * mhr0.nModTerm * log(treeVar * newExpVar)) -
            (mhr0.nTerm1 * mhr0.nModTerm * log(treeVar * m1Var)))); 

    // Interaction (newMixVar and mixVar cancels out if no new exposure)
    if (newMixVar != 0){ // - (new exposure of proposed)
      ratio -= 0.5 * log(treeVar * newMixVar) * mhr.nTerm1 * mhr0.nTerm2 * mhr0.nModTerm;
    }

    if (mixVar != 0){ // - (- mixVar of null)
      ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2 * mhr0.nModTerm;
    }

    // Accept / Reject
    if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
      mhr0 = mhr;
      success = 2;

      // If accepted with switch-exposure transition, update a tree more
      if (step1 == 3) { 
        // switch-exposure accept
        m1      = newExp;
        m1Var   = newExpVar;
        mixVar  = newMixVar;
        dlmTree1->replaceNodeVals(newTree);
      } else {
        // Grow/Prune/Change accept
        dlmTree1->accept(); // Otherwise, accept the tree
      }
      // Update parameters with new tree
      dlmTerm1 = dlmTree1->listTerminal(); // Get the new terminal nodes
      // No need to update Xmat for dlmtree accept
      for (Node* n : modTerm) {
        n->nodevals->updateXmat = 0;
      }

    } else {  // MH-ratio rejected
      modTree->setUpdateXmat(1);  // HDLM
    }
  }
  dlmTree1->reject();         // TDLM 

  // Reset new tree value for tree2
  if (newTree != 0)
    delete newTree;
  newTree = 0;

  // * Record tree 1
  if (ctr->diagnostics) {
    Eigen::VectorXd acc(9);
    acc << ctr->record, t, 1, step1, success, m1, dlmTerm1.size(), stepMhr, ratio;
    (dgn->treeDLMAccept).push_back(acc);
  }
  
  // [Propose a new TDLMM tree 2] - SI ==============================================
  // Rcout << "TreeMCMC: Updating tree 2...\n";
  newExp    = m2; 
  newExpVar = m2Var;
  newMixVar = mixVar;
  stepMhr = 0;
  success = 0;

  // Choose a transition and update
  step2 = sampleInt(ctr->stepProb, 1);
  if ((dlmTerm2.size() == 1) && (step2 < 3)){
    step2 = 0;
  }

  if(step2 < 3){ // Grow / Prune / Change
    stepMhr = tdlmProposeTree(dlmTree2, Exp[m2], ctr, step2);
    success = dlmTree2->isProposed();
    newDlmTerm2 = dlmTree2->listTerminal(1);

  } else { // Step 4: Propose a new exposure
    newExp = sampleInt(ctr->expProb); 
    if (newExp != m2) { // If the proposed exposure is different from the current one,
      // Update the information with new exposure
      success   = 1;
      newExpVar = ctr->muExp(newExp);     // Update the exposure list
      newTree   = new Node(*dlmTree2);    // Copy over the tree
      newTree->setUpdate(1);
      newDlmTerm2   = newTree->listTerminal();

      for (Node* nt : newDlmTerm2)
        Exp[newExp]->updateNodeVals(nt);

      // Update the interaction using the new exposure as well for TDLMMns/all
      if ((ctr->interaction) && ((ctr->interaction == 2) || (newExp != m1))) {
        if (newExp <= m1)
          newMixVar = ctr->muMix(m1, newExp);
        else
          newMixVar = ctr->muMix(newExp, m1);
      } else {
        newMixVar = 0;
      }
    }
  }

  // * dlmTree 2 MHR 
  // StepMHR: Transition ratio
  // mhr0 & mhr: Likelihood of R
  // Rcout << "TreeMCMC: Calculating MHR for tree 2...\n";

  // If successful, calculate MH ratio
  if (success) {
    newDlmTerm2 = dlmTree2->listTerminal(1);
    modTree->setUpdateXmat(1);

    // MH ratio with a new terminal and the exposure: newDlmTerm1, newExpVar
    // Rcout << "[Tree2: mhr] \n";
    mhr = dlmtreeTDLMM_MHR(modTerm, 
                           dlmTerm1, newDlmTerm2, ctr, ZtR, treeVar, 
                           m1Var, newExpVar, newMixVar);

    // MH ratio - dlmTree
    if (RtR < 0) {
      RtR = (ctr->R).dot(ctr->R);
      RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);
    }

    ratio = stepMhr + 
            mhr.logVThetaChol - mhr0.logVThetaChol -
            (0.5 * (ctr->n + 1.0) *
            (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
            log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
            (0.5 * ((mhr.nTerm2 * mhr0.nModTerm * log(treeVar * newExpVar)) -
            (mhr0.nTerm2 * mhr0.nModTerm * log(treeVar * m2Var))));

    // Interaction
    if (newMixVar != 0){ 
      ratio -= 0.5 * log(treeVar * newMixVar) * mhr0.nTerm1 * mhr.nTerm2 * mhr0.nModTerm;
    }

    if (mixVar != 0){ 
      ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2 * mhr0.nModTerm;
    }

    if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
      mhr0 = mhr;
      success = 2;

      // If accepted with switch-exposure transition, update a tree more
      if (step2 == 3) { 
        // switch-exposure accept
        m2      = newExp;
        m2Var   = newExpVar;
        mixVar  = newMixVar;
        dlmTree2->replaceNodeVals(newTree);
      } else {
        // Grow/Prune/Change accept
        dlmTree2->accept(); // Otherwise, accept the tree
      }

      // Update parameters with new tree
      dlmTerm2 = dlmTree2->listTerminal(); // Get the new terminal nodes
      // No need to update Xmat for dlmtree accept
      for (Node* n : modTerm) {
        n->nodevals->updateXmat = 0;
      }

    } else {  // MH-ratio rejected
      modTree->setUpdateXmat(1);  // HDLM
    }
  } 
  dlmTree2->reject();         // TDLM 

  if (newTree != 0)
    delete newTree;
  newTree = 0;

  // * Record tree 2
  if (ctr->diagnostics) {
    Eigen::VectorXd acc(9);
    acc << ctr->record, t, 2, step2, success, m2, dlmTerm2.size(), stepMhr, ratio;
    (dgn->treeDLMAccept).push_back(acc);
  }


  // **** Propose new modifier tree ****
  // std::vector<double> PruneChange = {ctr->stepProbMod[0], ctr->stepProbMod[1]}; // BIG UPDATE (Strict Prior) HERE
  // Rcout << "TreeMCMC: Modifier tree update...\n";
  switch (modTerm.size()) {
    case 1: step = 0;   break;
    case 2: step = sampleInt(ctr->stepProbMod, 1 - ctr->stepProbMod[3]); break;
    default: step = sampleInt(ctr->stepProbMod, 1);
  }

  stepMhr = modProposeTree(modTree, Mod, ctr, step);
  success = modTree->isProposed();

  if (success && (stepMhr == stepMhr)) {
    // New modifier node has been proposed 
    // -> X needs to be updated for MHR calculation
    newModTerm = modTree->listTerminal(1);

    // Rcout << "[modTree: mhr] \n";
    mhr = dlmtreeTDLMM_MHR(newModTerm, dlmTerm1, dlmTerm2, 
                            ctr, ZtR, treeVar,
                            m1Var, m2Var, mixVar);
                            
    ratio = stepMhr + mhr.logVThetaChol - mhr0.logVThetaChol -
      (0.5 * (ctr->n + 1.0) *
        (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
         log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2)));// -

    // Modifier grow
    if (step == 0){
      ratio -= 0.5 * (mhr0.nTerm1 * log(treeVar * m1Var) + mhr0.nTerm2 * log(treeVar * m2Var));

      if (mixVar != 0){ // interaction
        ratio -= 0.5 * mhr0.nTerm1 * mhr0.nTerm2 * log(treeVar * mixVar);
      }
    }

    // Modifier prune
    if (step == 1){
      ratio += 0.5 * (mhr0.nTerm1 * log(treeVar * m1Var) + mhr0.nTerm2 * log(treeVar * m2Var));

      if (mixVar != 0){ // interaction
        ratio += 0.5 * mhr0.nTerm1 * mhr0.nTerm2 * log(treeVar * mixVar);
      }
    }

    if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
      mhr0 = mhr;
      success = 2;
      modTree->accept();
      modTerm = modTree->listTerminal();
    }
  } // end modTree proposal
  modTree->reject();

  // * Record modifier tree
  if (ctr->diagnostics) {
    Eigen::VectorXd acc(9);
    acc << ctr->record, t, 0, step, success, 0, modTerm.size(), stepMhr, ratio;
    (dgn->treeDLMAccept).push_back(acc);
  }

  // **** Update parameters ****
  // Tree shrinkage
  double tauT2 = mhr0.term1T2 / m1Var + mhr0.term2T2 / m2Var;
  if(mixVar != 0){
    tauT2 += mhr0.mixT2 / mixVar;
  }

  if (ctr->shrinkage > 1) {
    double xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->tau)(t)));
    (ctr->tau)(t) = 1.0 / R::rgamma(0.5 * mhr0.nDlmTerm * mhr0.nModTerm + 0.5,
                                    1.0 / ((0.5 * tauT2 / (ctr->sigma2 * ctr->nu)) + xiInv));
  }



  // Rcout << "TreeMCMC: Rmat \n";
  ctr->Rmat.col(t) = mhr0.fitted;


  // Update exposure information
  // Rcout << "TreeMCMC: expCount, expInf, totTermExp, nTermMod, nTerm ... \n";
  (ctr->dlmTree1Exp)(t) = m1;
  (ctr->dlmTree2Exp)(t) = m2;
  (ctr->expCount)(m1)++;                  // Increment count of m1 exposure
  (ctr->expCount)(m2)++;                  // Increment count of m2 exposure
  (ctr->expInf)(m1) += ((ctr->tau)(t));   // Add tree variance
  (ctr->expInf)(m2) += ((ctr->tau)(t));   // Add tree variance
  (ctr->totTermExp)(m1) += mhr0.nModTerm * mhr0.nTerm1;   // Updating P_m
  (ctr->totTermExp)(m2) += mhr0.nModTerm * mhr0.nTerm2;   // Updating P_m
  (ctr->sumTermT2Exp)(m1) += mhr0.term1T2 / (ctr->tau)(t);  // Updating D_m
  (ctr->sumTermT2Exp)(m2) += mhr0.term2T2 / (ctr->tau)(t);  // Updating D_m

  // modTree / dlmtree numbers information
  ctr->nTermMod(t) = mhr0.nModTerm;
  ctr->nTermDLM(t) = mhr0.nDlmTerm;
  ctr->nTermDLM1(t) = mhr0.nTerm1;
  ctr->nTermDLM2(t) = mhr0.nTerm2;  

  // Rcout << "TreeMCMC: Mixture update ... \n";
  if (mixVar != 0) {
    if (m1 <= m2) {
      (ctr->mixCount)(m2, m1)++;
      (ctr->totTermMix)(m2, m1) += mhr0.nModTerm * mhr0.nTerm1 * mhr0.nTerm2;
      (ctr->sumTermT2Mix)(m2, m1) += mhr0.mixT2 / (ctr->tau)(t);
      (ctr->mixInf)(m2, m1) += ((ctr->tau)(t));
    } else {
      (ctr->mixCount)(m1, m2)++;
      (ctr->totTermMix)(m1, m2) += mhr0.nModTerm * mhr0.nTerm1 * mhr0.nTerm2;
      (ctr->sumTermT2Mix)(m1, m2) += mhr0.mixT2 / (ctr->tau)(t);
      (ctr->mixInf)(m1, m2) += ((ctr->tau)(t));
    }
  }
  
  // Rcout << "TreeMCMC: Modifier count in modTree ... \n";
  // **** Count modifiers used in tree ****
  Eigen::VectorXd modCount = countMods(modTree, Mod);
  ctr->modCount += modCount;

  // **** Record **** 
  // Rcout << "TreeMCMC: Record ... \n";
  if (ctr->record > 0) {    
    // Rcout << "Modifier count record ... \n";
    for (int i = 0; i < modCount.size(); i++) {
      if (modCount(i) > 0){
        ctr->modInf(i) += (ctr->tau)(t);
      }
    }

    // **** Update DLM partial estimate ****
    // Rcout << "Creating a rule ... \n";
    // Rule: Takes a node and modifier data and returns a string 
    // describing current series of modifier rules leading to terminal node
    std::string rule; 

    // 1. which MCMC iteration
    // 2. which tree of ensemble
    // 3. which terminal node of a modifier tree
    // 4. which exposure
    // 5. minimum value of exposure
    // 6. maximum value of exposure
    // 7. t_min
    // 8. t_max
    // 9. delta_a full-conditional sampled for a dlmTree node
    // 10. Exposure-specific variance

    Eigen::VectorXd rec(10); 
    Eigen::VectorXd mix(10);
    rec << ctr->record, t, 0, 0, 0, 0, 0, 0, 0, 0;
    mix << ctr->record, t, 0, 0, 0, 0, 0, 0, 0, 0;

    // For each terminal node of the modifier trees `s`,
    // Rcout << "Going through each modifier tree terminals ... \n";
    for (s = 0; s < modTerm.size(); s++) {
      // rule object to setup a string rules for outputs
      // Rcout << "Setting up a rule \n";
      rule = modRuleStr(modTerm[s], Mod); 

      // Iterate through dlmTree 1 and the mixture
      // Rcout << "Tree 1 mixture record \n";
      std::size_t k = 0; // Mixture count
      for (std::size_t i = 0; i < mhr0.nTerm1; i++) {
        rec[2] = s;                                    // Modifier terminal count
        rec[3] = 1;                                    // which tree pair
        rec[4] = i;                                    // terminal node count
        rec[5] = m1;                                   // Exposure
        rec[6] = (dlmTerm1[i]->nodestruct)->get(3);    // tmin
        rec[7] = (dlmTerm1[i]->nodestruct)->get(4);    // tmax
        rec[8] = mhr0.drawAll(s * mhr0.nDlmTerm + i);  // delta_a sample
        rec[9] = (ctr->tau)(t) * m1Var;                // Exposure-variance

        (dgn->termRule).push_back(rule); // Save the modifier terminal rule (location)
        (dgn->DLMexp).push_back(rec);    // Save the tree1 information
        
        for (std::size_t j = 0; j < mhr0.nTerm2; j++) {
          if (i == 0) {     
            rec[2] = s;     // Still the same terminal node of the modifier tree
            rec[3] = 2;                                    // which tree pair
            rec[4] = j;                                    // terminal node count
            rec[5] = m2;                                   // Exposure
            rec[6] = (dlmTerm2[j]->nodestruct)->get(3);    // tmin
            rec[7] = (dlmTerm2[j]->nodestruct)->get(4);    // tmax
            rec[8] = mhr0.drawAll(s * mhr0.nDlmTerm + mhr0.nTerm1 + j);
            rec[9] = (ctr->tau)(t) * m2Var;
            
            (dgn->termRule).push_back(rule); // Save the modifier terminal rule (location)
            (dgn->DLMexp).push_back(rec);    // Save the tree2 information
          }

          if (mixVar != 0) {
            mix[2] = s; // Still the same terminal node of the modifier tree
            if (m1 <= m2) {
              mix[3] = m1;
              mix[4] = (dlmTerm1[i]->nodestruct)->get(3); // tmin
              mix[5] = (dlmTerm1[i]->nodestruct)->get(4); // tmax
              mix[6] = m2;
              mix[7] = (dlmTerm2[j]->nodestruct)->get(3); // tmin
              mix[8] = (dlmTerm2[j]->nodestruct)->get(4); // tmax
            } else {
              mix[6] = m1;
              mix[7] = (dlmTerm1[i]->nodestruct)->get(3); // tmin
              mix[8] = (dlmTerm1[i]->nodestruct)->get(4); // tmax
              mix[3] = m2;
              mix[4] = (dlmTerm2[j]->nodestruct)->get(3); // tmin
              mix[5] = (dlmTerm2[j]->nodestruct)->get(4); // tmax
            }
            mix[9] = mhr0.drawAll(s * mhr0.nDlmTerm + mhr0.nTerm1 + mhr0.nTerm2 + k);
            
            (dgn->termRuleMIX).push_back(rule); // Save the modifier terminal rule (location) just for the mixture
            (dgn->MIXexp).push_back(mix);       // Save the mixture information
            k++;
          } // mixVar != 0 end
        } // Tree 2 iteration end
      } // Tree 1 iteration end
    } // modifier tree iteration end
  } // ctr->record > 0 end
  // Rcout << "TreeMCMC: Record finished ... \n";

} // end dlmtreeTDLMMGaussian_TreeMCMC function


treeMHR dlmtreeTDLMM_MHR(std::vector<Node*> modTerm,    // Modifier terminal node -> Need this for calculating different effect size
                          std::vector<Node*> dlmTerm1,  // Just need one set of dlm Terminal nodes from tree1 as this is a shared dlmtree
                          std::vector<Node*> dlmTerm2,  // Just need one set of dlm Terminal nodes from tree2 as this is a shared dlmtree
                          dlmtreeCtr* ctr, 
                          Eigen::VectorXd ZtR, 
                          double treeVar, double m1Var, double m2Var, double mixVar)
// Calculate part of Metropolis-Hastings ratio and make draws from full conditional. 
{
  std::size_t s;    // Count integer
  treeMHR out;      // Storage object

  // exposure variance
  out.m1Var = m1Var;
  out.m2Var = m2Var;

  // modifier tree
  // Rcout << "Counting the number of modifier terminal nodes \n";
  int pXMod = modTerm.size();

  // Number of terminal nodes
  // Rcout << "Counting the number of dlmtree terminal nodes \n";
  int pXDlm1 = dlmTerm1.size();
  int pXDlm2 = dlmTerm2.size();
  int pXDlm = pXDlm1 + pXDlm2;
  int interaction = 0;
  if (mixVar != 0) { 
    pXDlm += pXDlm1 * pXDlm2; // Add (nodes x nodes) to the parameter # (the rectangle figure)
    interaction = 1;          // Interaction flag
  }

  // Rcout << "Counting the totoal number of dlmtree x modifier terminal nodes \n";
  int pXComb = pXMod * pXDlm;   // Mod * (Tree1 + Tree2 + Tree1x2)

  // Rcout << "Define parameters for calculation \n";
  Eigen::MatrixXd Xd, ZtX, VgZtX;
  Eigen::VectorXd diagVar;
  Xd.resize(ctr->n, pXDlm);         Xd.setZero(); // Same terminal for all modifiers
  ZtX.resize(ctr->pZ, pXComb);      ZtX.setZero();
  VgZtX.resize(ctr->pZ, pXComb);    VgZtX.setZero();
  diagVar.resize(pXDlm);            diagVar.setZero(); // Need a diagVar as we have multiple exposures
  
  // Building exposure matrix, Xd & ZtX
  int i, j, k;
  // dlmTree 1
  // Rcout << "TDLMM tree 1 \n";
  for (i = 0; i < pXDlm1; i++) {                // iterate through terminal nodes of tree 1
    Xd.col(i) = (dlmTerm1[i]->nodevals)->X; // Each column of Xd is a vector X associated with a node from nodes1
    diagVar(i) = 1.0 / (m1Var * treeVar);       // Exposure specific variance
    ZtX.col(i) = (dlmTerm1[i]->nodevals)->ZtX;
  }

  // dlmTree 2
  // Rcout << "TDLMM tree 2 \n";
  for (j = 0; j < pXDlm2; j++) {                // iterate through terminal nodes of tree 2
    k = pXDlm1 + j;                             // k is shifted back because pX2 needs to be after pX1.
    Xd.col(k) = (dlmTerm2[j]->nodevals)->X; 
    diagVar(k) = 1.0 / (m2Var * treeVar);       // Exposure specific variance
    ZtX.col(k) = (dlmTerm2[j]->nodevals)->ZtX;
  }

  // Mixture & interaction
  if(interaction) {
    for (i = 0; i < pXDlm1; i++) {
      for (j = 0; j < pXDlm2; j++) {
        // tree1's first terminal with all terminals of tree2, 
        // tree1's second terminal with all terminals of tree2, and so on.
        k = (pXDlm1 + pXDlm2) + i * pXDlm2 + j; 
        Xd.col(k) = (((dlmTerm1[i]->nodevals)->X).array() * ((dlmTerm2[j]->nodevals)->X).array()).matrix();
        diagVar(k) = 1.0 / (mixVar * treeVar);
        ZtX.col(k) = ctr->Zw.transpose() * Xd.col(k); 
      }
    }
  }

  // MH ratio element calculation
  if(pXMod == 1){   // Modifier has one node -> No subsetting required -> same code as TDLMM    
    // Rcout << "Modifier has only one terminal node \n";
    const Eigen::MatrixXd VgZtX = ctr->Vg * ZtX; // V_gamma * Z_t * X = (pxp)(pxn)(nx1) = (px1)
    Eigen::MatrixXd tempV(pXDlm, pXDlm);        // tempV: Xa^T * V_z * Xa (first part of V_theta, note that inverse is the same)
    Eigen::VectorXd XtVzInvR; // was ctr->n

    // V_theta = V_delta_a from eq(11) & Xd = X_a in paper
    tempV.triangularView<Eigen::Lower>() = Xd.transpose() * Xd;   // (dlm# x dlm#)
    tempV.noalias() -= ZtX.transpose() * VgZtX; 
    out.tempV = tempV;
    
    // Finalize Eq.(9)
    tempV.diagonal().noalias() += diagVar;    // tempV is now V_theta(a)_inverse
    XtVzInvR = Xd.transpose() * ctr->R - VgZtX.transpose() * ZtR; // (dlm# x 1)

    // Calculate Vtheta & Vtheta cholesky
    // inverse tempV to obtain V_theta(a)_inverse
    // Rcout << "pxMod == 1: Vtheta calculation \n";
    Eigen::MatrixXd VTheta(pXDlm, pXDlm); 
    VTheta.triangularView<Eigen::Lower>() =
      tempV.selfadjointView<Eigen::Lower>().llt().solve(Eigen::MatrixXd::Identity(pXDlm, pXDlm));
    const Eigen::MatrixXd VThetaChol = VTheta.selfadjointView<Eigen::Lower>().llt().matrixL();

    // Sampling from Supplemental Eq(11) & (12) (Theta = Delta_a)
    // Mean
    // Rcout << "pxMod == 1: Mean of Delta_a \n";
    const Eigen::VectorXd ThetaHat = VTheta.selfadjointView<Eigen::Lower>() * XtVzInvR; // Eq. (11) mean

    // Store the sampled values and also add variance
    // Rcout << "pxMod == 1: Variance of Delta_a \n";
    Eigen::VectorXd ThetaDraw = ThetaHat; // (pXDlm x 1) vector
    ThetaDraw.noalias() += VThetaChol * as<Eigen::VectorXd>(rnorm(pXDlm, 0, sqrt(ctr->sigma2)));

    // Sort out Equation 11 Full conditional draw
    out.drawAll = ThetaDraw;

    // Temporary vectors for calculation below
    // Rcout << "pxMod == 1: drawTemp \n";
    Eigen::VectorXd drawTemp(pXDlm);      // Temporary draw vector for accessing the block
    Eigen::VectorXd draw1Temp(pXDlm1);    // Can resize with a fixed value since this is a shared tree
    Eigen::VectorXd draw2Temp(pXDlm2);    // Can resize with a fixed value since this is a shared tree

    // Exposure 1 return
    // Rcout << "pxMod == 1: Exposure 1 segment \n";
    draw1Temp = ThetaDraw.head(pXDlm1);             // column (1 ~ pX1) => DLM effect of exposure 1 sampled
    out.term1T2 = (draw1Temp).dot(draw1Temp);      // draw1 dot product for MHR math
    out.nTerm1 = double(pXDlm1);                    // Number of terminal nodes of tree 1

    // Exposure 2 return
    // Rcout << "pxMod == 1: Exposure 2 segment \n";
    draw2Temp = ThetaDraw.segment(pXDlm1, pXDlm2);  // column (pX1 ~ pX2) => DLM effect of exposure 2 sampled
    out.term2T2 = (draw2Temp).dot(draw2Temp);      // draw2 dot product for MHR math
    out.nTerm2 = double(pXDlm2);                    // Number of terminal nodes of tree 2

    // Mixture return
    if (interaction) {    
      // Rcout << "pxMod == 1: Interaction/Mixture segment \n";
      out.drawMix = ThetaDraw.tail(pXDlm - pXDlm1 - pXDlm2);  // Extract last (pX1 x pX2) element for mixture
      out.mixT2 = (out.drawMix).dot(out.drawMix);             // Mixture dot product for MHR math
      out.nTermMix = double(pXDlm - pXDlm1 - pXDlm2);
    }

    // Calculate and store parameters for MHR math 
    // Rcout << "pxMod == 1: MH ratio math calculation \n";
    out.fitted = Xd * out.drawAll;  // exposure matrix x delta

    out.beta = ThetaHat.dot(XtVzInvR);  // The big chunk in Eq(8)
    out.logVThetaChol = VThetaChol.diagonal().array().log().sum(); // |V_theta_a|
    out.termT2 = (out.drawAll).dot(out.drawAll);
    
    out.nDlmTerm = pXDlm * 1.0; 
    out.nModTerm = pXMod * 1.0;    // = 1

    return(out);
  } // End of modifier tree with only one node

  // **** Multiple Modifier nodes ****
  // Rcout << "Multiple modifier nodes -> Build blocks by subetting with the modifier tree \n";
  // Xtemp, Ztemp, Rtemp will be updated throughout the modifier group
  Eigen::MatrixXd Xtemp, Ztemp;
  Eigen::VectorXd Rtemp;
  Eigen::MatrixXd XtXblock(pXComb, pXComb);    XtXblock.setZero(); // BC x BC
  Eigen::VectorXd XtR(pXComb);                 XtR.setZero();       // BC
  Eigen::MatrixXd LInv(pXDlm, pXDlm);          LInv.setZero();      // C x C
  LInv.diagonal().array() += diagVar.array();

  // Create block matrices corresponding to modifier nodes
  // Rcout << "Subsetting data and exposures according to modifier terminal nodes \n";
  int start = 0;
  for (Node* n : modTerm) { // For each terminal node of a modifier tree,
    //if (n->nodevals->updateXmat) { // update matrices for current node 
    if(true){
      // Retrieve indices subset by the modifier tree
      Xtemp.resize(n->nodevals->idx.size(), pXDlm);     Xtemp.setZero();  // Exposure subset
      Ztemp.resize(n->nodevals->idx.size(), ctr->pZ);   Ztemp.setZero();  // Fixed effect
      Rtemp.resize(n->nodevals->idx.size());            Rtemp.setZero();  // Partial residual
      
      // Loop through indices and (subset)update Xtemp, Ztemp, Rtemp
      j = 0;
      for (int i : n->nodevals->idx) {
        Xtemp.row(j) = Xd.row(i);   // Extract from Xd and store to Xtemp -> includes interaction here
        Ztemp.row(j) = ctr->Z.row(i);   // Extract from Z and store to Ztemp
        Rtemp(j) = ctr->R(i);           // Extract from R(vector) and store to Rtemp
        j++;
      } // end loop over node indices
       
      n->nodevals->XtX.resize(pXDlm, pXDlm);
      n->nodevals->XtX = Xtemp.transpose() * Xtemp;
      n->nodevals->ZtXmat.resize(ctr->pZ, pXDlm);
      n->nodevals->ZtXmat = Ztemp.transpose() * Xtemp;
      n->nodevals->VgZtXmat.resize(ctr->pZ, pXDlm);
      n->nodevals->VgZtXmat = ctr->Vg * n->nodevals->ZtXmat;
      n->nodevals->updateXmat = 0; // Xmat has been updated
      
      XtR.segment(start, pXDlm) = Xtemp.transpose() * Rtemp;
      
    } else { // No need to update -> go through idx and update XtR
      for (int i : n->nodevals->idx) {
        XtR.segment(start, pXDlm).noalias() += Xd.row(i).transpose() * ctr->R(i);
      } // end loop over node indices
    } // end update xblock and ztx block (.block = accessing a block in a matrix)

    XtXblock.block(start, start, pXDlm, pXDlm) = ((n->nodevals->XtX) + LInv).inverse();
    ZtX.block(0, start, ctr->pZ, pXDlm) = n->nodevals->ZtXmat;
    VgZtX.block(0, start, ctr->pZ, pXDlm) = n->nodevals->VgZtXmat;
    
    // Move to the next block
    start += pXDlm;
  } // End of subsetting and building blocks

  // Compute elements of MH ratio
  // Rcout << "pxMod != 1: VTheta / ThetaHat calculation \n";
  Eigen::MatrixXd ZtXXi = ZtX * XtXblock;
  Eigen::MatrixXd VTheta = XtXblock;
  VTheta.noalias() += ZtXXi.transpose() * 
    (ctr->VgInv - ZtXXi * ZtX.transpose()).inverse() * ZtXXi;
  Eigen::VectorXd XtVzInvR = XtR;
  XtVzInvR.noalias() -= VgZtX.transpose() * ZtR;
  Eigen::VectorXd ThetaHat = VTheta * XtVzInvR;
  Eigen::MatrixXd VThetaChol = VTheta.llt().matrixL();

  // Rcout << "Sorting thetahat \n";
  // Store the sampled values and also add variance 
  // Rcout << "pxMod != 1: Sampling ThetaDraw Mean + Variance and saving as drawAll\n";
  out.drawAll = ThetaHat;
  out.drawAll.noalias() += 
      VThetaChol * as<Eigen::VectorXd>(rnorm(pXComb, 0.0, sqrt(ctr->sigma2)));

  // Sort out Equation 11 Full conditional draw
  // drawAll = mod1-dlm1 / mod1-dlm2 / mod1-dlm1&2 / mod2-dlm1 / mod2-dlm2 / mod2-dlm1&2 / mod3 ...
  // Note that we do not have draw1/draw2/drawMix when we have a modifier tree

  // Fitted value & terminal effect draws
  // Rcout << "pxMod != 1: Setting up temporary drawTheta segment (main effect) \n";
  out.fitted.resize(ctr->n);                  // Estimation
  Eigen::VectorXd drawTemp(pXDlm);            // Temporary draw vector for accessing the block
  Eigen::VectorXd draw1Temp(pXDlm1);          // Can resize with a fixed value since this is a shared tree
  Eigen::VectorXd draw2Temp(pXDlm2);          // Can resize with a fixed value since this is a shared tree
  Eigen::VectorXd drawMixTemp(pXDlm1 * pXDlm2);
  
  if(interaction){  // Prevent resizing drawMixTemp to zero in a case of HDLMMadd
    // Rcout << "pxMod != 1: Setting up temporary drawTheta segment (interaction effect) \n";
    drawMixTemp.resize(pXDlm - pXDlm1 - pXDlm2);
  }

  // Store Term1T2/Term2T2/MixT2 & nTerm1/nTerm2 for the MHR math
  // Rcout << "pxMod != 1: Setting up MH ratio parameters \n";
  out.term1T2 = 0;
  out.term2T2 = 0;
  out.mixT2 = 0;
  
  // Rcout << "pxMod != 1: Iterating through blocks for each modifier terminal node \n";
  for (s = 0; s < modTerm.size(); s++) { // Iterate through each mod terms to calculate the corresponding dlm
    // Extract draws from the block matrix for each modifier terminal node
    drawTemp = out.drawAll.segment(s * pXDlm, pXDlm); // subset drawAll by a corresponding modifier terminal
    
    // Note: 
    // 1. The code uses fixed pXDlm1, pXDlm2 as they do not change due to shared tree setting
    // 2. But the effects are different from the different blocks

    // Tree 1 draw 
    draw1Temp = drawTemp.head(pXDlm1);
    out.term1T2 += (draw1Temp).dot(draw1Temp);

    // Tree 2 draw
    draw2Temp = drawTemp.segment(pXDlm1, pXDlm2);
    out.term2T2 += (draw2Temp).dot(draw2Temp);

    // [Mixture] Tree 1 x Tree 2 draw
    if (interaction) {
      drawMixTemp = drawTemp.tail(pXDlm - pXDlm1 - pXDlm2);
      out.mixT2 += (drawMixTemp).dot(drawMixTemp);
      out.nTermMix = drawMixTemp.size();
    }

    // Calculate fitted for each modifier terminal
    for (int i : modTerm[s]->nodevals->idx){           // all indices that belong to mod terminal node s
      out.fitted(i) = Xd.row(i) * drawTemp;       // calculate exposure x corresponding dlm
    }
  }

  // Calculation return
  // Rcout << "pxMod != 1: MH ratio math calculation \n";
  out.beta = ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
  out.termT2 = (out.drawAll).dot(out.drawAll);

  out.nTerm1 = double(pXDlm1);
  out.nTerm2 = double(pXDlm2);
  out.nDlmTerm = pXDlm * 1.0; 
  out.nModTerm = pXMod * 1.0;

  return(out);
}
