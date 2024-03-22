/**
 * @file dlmtreeHDLMMGaussian.cpp
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
void dlmtreeHDLMMGaussian_TreeMCMC(int t, NodeStruct* expNS,
                                   Node* modTree, 
                                   Node* dlmTree1, Node* dlmTree2,
                                   dlmtreeCtr* ctr, dlmtreeLog *dgn,
                                   modDat* Mod, std::vector<exposureDat*> Exp);

// MHR updated
// 1. Tree pair: Tree1 & Tree 2
// 2. Exposure-variance: m1Var, m2Var, mixVar                      
treeMHR dlmtreeHDLMM_MHR(std::vector<Node*> modTerm,
                         std::vector<Node*> dlmTerm1, 
                         std::vector<Node*> dlmTerm2,
                         dlmtreeCtr* ctr, Eigen::VectorXd ZtR, 
                         double treeVar, double m1Var, double m2Var, double mixVar);

// [[Rcpp::export]]
Rcpp::List dlmtreeHDLMMGaussian(const Rcpp::List model){ 
  // *** Set up general control variables ***
  dlmtreeCtr *ctr = new dlmtreeCtr;

  // MCMC parameters
  ctr->iter = as<int>(model["nIter"]); 
  ctr->burn = as<int>(model["nBurn"]); 
  ctr->thin = as<int>(model["nThin"]); 
  ctr->nRec = floor(ctr->iter / ctr->thin);
  ctr->nTrees = as<int>(model["nTrees"]);  

  // Data setup
  ctr->Y = as<Eigen::VectorXd>(model["Y"]);
  ctr->n = (ctr->Y).size();

  // Modifier tree hyperparameter
  ctr->modZeta = as<double>(model["zeta"]); 
  ctr->modKappa = 100;        

  // Model estimation initialization
  ctr->Z = as<Eigen::MatrixXd>(model["Z"]); 
  ctr->Zw = ctr->Z; 
  ctr->pZ = (ctr->Z).cols(); 

  ctr->VgInv = (ctr->Z).transpose() * (ctr->Z); 
  ctr->VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg = ctr->VgInv.inverse();  
  ctr->VgChol = (ctr->Vg).llt().matrixL();      

  // Binomial flag              
  ctr->binomial = as<bool>(model["binomial"]);  
  ctr->zinb = as<bool>(model["zinb"]);  

  // Tree prior probabilities for the modifier tree and DLM tree
  ctr->stepProbMod = as<std::vector<double>>(model["stepProbMod"]);
  ctr->stepProb = as<std::vector<double>>(model["stepProbTDLM"]);
  ctr->treePriorMod = as<std::vector<double>>(model["treePriorMod"]); 
  ctr->treePrior = as<std::vector<double>>(model["treePriorTDLM"]);


  // Diagnostics & messages
  ctr->verbose = bool(model["verbose"]); 
  ctr->diagnostics = bool(model["diagnostics"]);

  // [Mixture & Shrinkage] 
  // Store shrinkage: 
  // 3 = all (mu(s), tau(t))
  // 2 = trees (tau(t))
  // 1 = exposures/interactions (mu(s)) - (default)
  // 0 = none
  ctr->shrinkage = as<int>(model["shrinkage"]);

  // Sparsity hyperparameter (Deprecated now as HDLMM does not perform exposure selection)
  ctr->mixKappa = as<double>(model["mixPrior"]); 
  bool updateKappa = true;  
  if (ctr->mixKappa < 0) {       
    updateKappa = true;
    ctr->mixKappa = 1;
  }

  // *** Setup modifier data ***
  modDat *Mod = new modDat(as<std::vector<int>>(model["modIsNum"]), 
                           as<Rcpp::List>(model["modSplitIdx"]), 
                           as<std::vector<int>>(model["fullIdx"]));

  NodeStruct *modNS; 
  modNS = new ModStruct(Mod); 
  ctr->pM = Mod->nMods; 

  // *** Pre-calculate single node tree pair matrices ***
  std::vector<exposureDat*> Exp;                     
  Rcpp::List exp_dat = as<Rcpp::List>(model["X"]); 
  ctr->nExp = exp_dat.size();                  
  for (int i = 0; i < ctr->nExp; i++) {
    Exp.push_back(
      new exposureDat(
        as<Eigen::MatrixXd>(
          as<Rcpp::List>(exp_dat[i])["Tcalc"]), ctr->Z, ctr->Vg));
  }

  ctr->pX = Exp[0]->pX; 
  ctr->nSplits = 0;

  // *** Mixture / interaction management ***
  ctr->interaction = as<int>(model["interaction"]);
  ctr->nMix = 0;
  if (ctr->interaction) { 
    ctr->nMix += int (ctr->nExp * (ctr->nExp - 1.0) / 2.0); 
    if (ctr->interaction == 2) {
      ctr->nMix += ctr->nExp;
    }
  }

  // *** Create trees ***
  NodeStruct *expNS;
  expNS = new DLNMStruct(0, ctr->nSplits + 1,
                         1, int (ctr->pX), 
                         as<Eigen::VectorXd>(model["splitProb"]),
                         as<Eigen::VectorXd>(model["timeProb"]));

  // Exposure information vectors
  ctr->expProb = as<Eigen::VectorXd>(model["expProb"]);

  (ctr->dlmTree1Exp).resize(ctr->nTrees);    
  (ctr->dlmTree2Exp).resize(ctr->nTrees);

  // Fixed exposures for hdlmm
  (ctr->expCount).resize((ctr->expProb).size()); 
  (ctr->expInf).resize((ctr->expProb).size());

  // *** Tree pair ***
  int t; // Tree index
  
  std::vector<Node*> modTrees;
  std::vector<Node*> dlmTrees1;
  std::vector<Node*> dlmTrees2;

  // For loop iterating through trees to start the "roots"
  for (t = 0; t < ctr->nTrees; t++) {
    // Modifier tree initialize
    modTrees.push_back(new Node(0, 1));         
    modTrees[t]->nodestruct = modNS->clone();   
    Mod->updateNodeVals(modTrees[t]);           

    ctr->dlmTree1Exp(t) = sampleInt(ctr->expProb);
    ctr->dlmTree2Exp(t) = sampleInt(ctr->expProb);
    dlmTrees1.push_back(new Node(0, 1));        
    dlmTrees2.push_back(new Node(0, 1));        
    dlmTrees1[t]->nodestruct = expNS->clone();  
    dlmTrees2[t]->nodestruct = expNS->clone();  
    Exp[ctr->dlmTree1Exp(t)]->updateNodeVals(dlmTrees1[t]); 
    Exp[ctr->dlmTree2Exp(t)]->updateNodeVals(dlmTrees2[t]); 
  }

  delete modNS;
  // delete expNS;

  // *** Logs ***
  dlmtreeLog *dgn = new dlmtreeLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec);    (dgn->gamma).setZero();  
  (dgn->sigma2).resize(ctr->nRec);            (dgn->sigma2).setZero(); 
  (dgn->nu).resize(ctr->nRec);                (dgn->nu).setZero();     
  (dgn->tau).resize(ctr->nTrees, ctr->nRec);  (dgn->tau).setZero();    
  (dgn->fhat).resize(ctr->n);                 (dgn->fhat).setZero();   
  (dgn->totTerm).resize(ctr->nRec);           (dgn->totTerm).setZero();

  // Exposure-specific shrinkage
  (dgn->muExp).resize(ctr->nExp, ctr->nRec);        (dgn->muExp).setZero();
  if (ctr->interaction > 0) { 
    (dgn->muMix).resize(ctr->nMix, ctr->nRec);      (dgn->muMix).setZero();   
    (dgn->mixInf).resize(ctr->nMix, ctr->nRec);     (dgn->mixInf).setZero(); 
    (dgn->mixCount).resize(ctr->nMix, ctr->nRec);   (dgn->mixCount).setZero();
  } else { 
    (dgn->muMix).resize(1, 1);                      (dgn->muMix).setZero();
    (dgn->mixInf).resize(1, 1);                     (dgn->mixInf).setZero(); 
    (dgn->mixCount).resize(1, 1);                   (dgn->mixCount).setZero();
  }

  // Modifier tree log
  (dgn->modProb).resize(ctr->pM, ctr->nRec);  (dgn->modProb).setZero();   
  (dgn->modCount).resize(ctr->pM, ctr->nRec); (dgn->modCount).setZero();  
  (dgn->modInf).resize(ctr->pM, ctr->nRec);   (dgn->modInf).setZero();  
  (dgn->modKappa).resize(ctr->nRec);          (dgn->modKappa).setZero(); 
  (dgn->termNodesMod).resize(ctr->nTrees, ctr->nRec);   (dgn->termNodesMod).setZero();

  // Exposure log
  (dgn->expProb).resize(ctr->nExp, ctr->nRec);          (dgn->expProb).setZero();       
  (dgn->expCount).resize(ctr->nExp, ctr->nRec);         (dgn->expCount).setZero();      
  (dgn->expInf).resize(ctr->nExp, ctr->nRec);           (dgn->expInf).setZero();        
  (dgn->dlmTree1Exp).resize(ctr->nTrees, ctr->nRec);    (dgn->dlmTree1Exp).setZero();   
  (dgn->dlmTree2Exp).resize(ctr->nTrees, ctr->nRec);    (dgn->dlmTree2Exp).setZero();   
  (dgn->termNodesDLM1).resize(ctr->nTrees, ctr->nRec);  (dgn->termNodesDLM1).setZero(); 
  (dgn->termNodesDLM2).resize(ctr->nTrees, ctr->nRec);  (dgn->termNodesDLM2).setZero(); 
  (dgn->mixKappa).resize(ctr->nRec);                    (dgn->mixKappa).setZero();      

  // *** Initial draws ***
  (ctr->fhat).resize(ctr->n);      
  (ctr->fhat).setZero();   
  ctr->R = ctr->Y;                  
  (ctr->gamma).resize(ctr->pZ); 

  // *** Exposure specific parameters ***
  (ctr->totTermExp).resize(ctr->nExp);      (ctr->totTermExp).setZero();   
  (ctr->sumTermT2Exp).resize(ctr->nExp);    (ctr->sumTermT2Exp).setZero(); 
  (ctr->muExp).resize(ctr->nExp);           (ctr->muExp).setOnes();        

  // If there is an interaction, define more parameters
  if (ctr->interaction) { 
    (ctr->totTermMix).resize(ctr->nExp, ctr->nExp);   (ctr->totTermMix).setZero();   
    (ctr->sumTermT2Mix).resize(ctr->nExp, ctr->nExp); (ctr->sumTermT2Mix).setZero(); 
    (ctr->muMix).resize(ctr->nExp, ctr->nExp);        (ctr->muMix).setOnes();        
    (ctr->mixInf).resize(ctr->nExp, ctr->nExp);       (ctr->mixInf).setZero();       
    (ctr->mixCount).resize(ctr->nExp, ctr->nExp);     (ctr->mixCount).setZero();     
  }

  // Initialize hyperparameters
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  ctr->nu = 1.0;
  ctr->sigma2 = 1.0;

  // Initial model fitting
  tdlmModelEst(ctr); 

  // Shrinkage hyperparameters
  rHalfCauchyFC(&(ctr->nu), ctr->nTrees, 0.0);           
  (ctr->tau).resize(ctr->nTrees);  (ctr->tau).setOnes(); 
  if (ctr->shrinkage > 1){
    for (t = 0; t < ctr->nTrees; t++) {
      rHalfCauchyFC(&((ctr->tau)(t)), 0.0, 0.0);
    }
  }

  // *** Initializing trees ***
  // DLM tree pairs
  ctr->nTermDLM.resize(ctr->nTrees);     (ctr->nTermDLM).array() = 2; 
  ctr->nTermDLM1.resize(ctr->nTrees);    (ctr->nTermDLM1).array() = 1;
  ctr->nTermDLM2.resize(ctr->nTrees);    (ctr->nTermDLM2).array() = 1;

  // Modifier trees
  ctr->nTermMod.resize(ctr->nTrees);   (ctr->nTermMod).array() = 1; 
  ctr->modCount.resize(ctr->pM);       (ctr->modCount).setZero();   
  ctr->modInf.resize(ctr->pM);         (ctr->modInf).setZero();     

  // Partial residual matrix
  (ctr->Rmat).resize(ctr->n, ctr->nTrees);  (ctr->Rmat).setZero();
  
  // *** MCMC ***
  // Create a progress meter
  progressMeter* prog = new progressMeter(ctr);

  // Thinning and burn-in process
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }

    // *** Update trees and parameters ***
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
    for (t = 0; t < ctr->nTrees; t++) {
      dlmtreeHDLMMGaussian_TreeMCMC(t, expNS,              
                                    modTrees[t],          
                                    dlmTrees1[t],       
                                    dlmTrees2[t], 
                                    ctr, dgn, Mod, Exp);
      
      // Update the fitted values after each iteration of t 
      ctr->fhat += (ctr->Rmat).col(t);

      // For each tree pair, update the partial residual
      if (t < ctr->nTrees - 1)
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
    }

    // Rcout << "Tree ensemble update finished \n";
    // Update the parameters
    ctr->R = ctr->Y - ctr->fhat;
    ctr->sumTermT2 = (ctr->sumTermT2Exp).sum();
    ctr->totTerm = (ctr->totTermExp).sum();
    if(ctr->interaction > 0) {
      ctr->totTerm += (ctr->totTermMix).sum();
      ctr->sumTermT2 += (ctr->sumTermT2Mix).sum();
    }

    // Update gamma
    tdlmModelEst(ctr);

    // Exposure Shrinkage update
    // Shrinkage
    // 2 = trees (tau(t))
    // 1 = exposures/interactions (mu(s)) - (default)
    // 0 = none
    // Global shrinkage: nu
    double xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->nu)));
    ctr->nu = 1.0 / R::rgamma(0.5 * ctr->totTerm + 0.5,
                              1.0 / (0.5 * ctr->sumTermT2 / (ctr->sigma2) + xiInv));
                              
    // Exposure shrinkage
    double sigmanu = ctr->sigma2 * ctr->nu;
    if (ctr->shrinkage == 1 || ctr->shrinkage == 3) { 
      for (int i = 0; i < ctr->nExp; i++) {
        xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->muExp(i))));
        ctr->muExp(i) = 1.0 / R::rgamma(0.5 * ctr->totTermExp(i) + 0.5,
                        1.0 / (0.5 * ctr->sumTermT2Exp(i) / sigmanu + xiInv));

        // Similarly update interaction-specific terms
        if (ctr->interaction) {
          for (int j = i; j < ctr->nExp; j++) {
            if ((j > i) || (ctr->interaction == 2)){
              rHalfCauchyFC(&(ctr->muMix(j, i)), ctr->totTermMix(j, i),
                            ctr->sumTermT2Mix(j, i) / sigmanu);
            }
          } // end for loop updating interaction variances
        } // end if interactions
      } // end for loop updating exposure variances
    } // end if shrinkage == 1 or 3

    // *** Update modifier selection ***
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

      Mod->modProb = rDirichlet((ctr->modCount.array() + ctr->modKappa / ctr->pM).matrix());  // (m_1 + kappa / J, ..., m_J + kappa / J)
      // end modifier selection

      // HDLMM: Exposure selection posterior calculation (Deprecated)
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

    // *** Record every iteration ***
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

    // Progress mark
    prog->printMark();
  } // end MCMC
  // Rcout << "MCMC complete \n";


  // *** Prepare outout ***
  // Rcout << "Preparing output \n";
  Eigen::MatrixXd TreeStructs((dgn->DLMexp).size(), 10);    
  Rcpp::StringVector termRule(dgn->termRule.size());        
  Rcpp::StringVector termRuleMIX(dgn->termRuleMIX.size());  

  // Store rec() vectors
  std::size_t s; 
  for (s = 0; s < (dgn->DLMexp).size(); s++){ 
    TreeStructs.row(s) = dgn->DLMexp[s];
  }

  // Modifier rule per main & interaction effect
  termRule = dgn->termRule;
  termRuleMIX = dgn->termRuleMIX;

  // Transfer variables from the dgn object
  // DLM tree pair
  Eigen::VectorXd sigma2 = dgn->sigma2;
  Eigen::VectorXd nu = dgn->nu;
  Eigen::MatrixXd tau = (dgn->tau).transpose();
  Eigen::VectorXd totTerm = dgn->totTerm;
  Eigen::VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  Eigen::MatrixXd gamma = (dgn->gamma).transpose();

  // dlmTree information
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

  // If interaction, expand muMIX and MIX
  if (ctr->interaction) {
    muMix.resize((dgn->muMix).cols(), (dgn->muMix).rows());
    muMix = (dgn->muMix).transpose();
    MIX.resize((dgn->MIXexp).size(), 10);
    for (s = 0; s < (dgn->MIXexp).size(); ++s)
      MIX.row(s) = dgn->MIXexp[s];
  }

  // Modifier tree
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
  for (s = 0; s < Exp.size(); s++){       
    delete Exp[s];
  }
  delete Mod;
  for (s = 0; s < modTrees.size(); s++) { 
    delete modTrees[s];
    delete dlmTrees1[s];
    delete dlmTrees2[s];
  }

  return(Rcpp::List::create(Named("TreeStructs") = wrap(TreeStructs), 
                            Named("MIX") = wrap(MIX),
                            Named("termRules") = wrap(termRule),
                            Named("termRuleMIX") = wrap(termRuleMIX),
                            //Named("fhat") = wrap(fhat),
                            Named("sigma2") = wrap(sigma2),
                            Named("nu") = wrap(nu),
                            Named("tau") = wrap(tau),
                            //Named("totTerm") = wrap(totTerm),
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
                            Named("muMix") = wrap(muMix),
                            Named("modCount") = wrap(modCount),
                            Named("modInf") = wrap(modInf),
                            //Named("treeModAccept") = wrap(modAccept),
                            Named("treeDLMAccept") = wrap(dlmAccept)));

} // end dlmtreeTDLMMGaussian



void dlmtreeHDLMMGaussian_TreeMCMC(int t, NodeStruct* expNS, Node* modTree, 
                                   Node* dlmTree1, Node* dlmTree2,
                                   dlmtreeCtr* ctr, dlmtreeLog *dgn,
                                   modDat* Mod, std::vector<exposureDat*> Exp)
{ 

  int m1, m2, newExp;                                 
  int step, step1, step2;                             
  double m1Var, m2Var, mixVar, newExpVar, newMixVar;  
  int success = 0;                                    
  double RtR = -1.0;                                  
  double RtZVgZtR = 0;                                
  double stepMhr = 0;                                 
  double ratio = 0;                                   
  double treeVar = (ctr->nu) * (ctr->tau)(t);         
  std::size_t s;                                      
  std::vector<Node*> modTerm, dlmTerm1, dlmTerm2;          
  std::vector<Node*> newModTerm, newDlmTerm1, newDlmTerm2; 
  Node* newTree = 0;                                  
  treeMHR mhr0, mhr;                                  

  // Pre-calculation for MH ratio update
  Eigen::VectorXd ZtR = (ctr->Z).transpose() * (ctr->R);

  // -- List terminal nodes --
  modTerm = modTree->listTerminal();   
  dlmTerm1 = dlmTree1->listTerminal(); 
  dlmTerm2 = dlmTree2->listTerminal(); 

  // Extract exposure information from ctr
  m1 = ctr->dlmTree1Exp[t]; 
  m2 = ctr->dlmTree2Exp[t]; 
  m1Var = ctr->muExp(m1);   
  m2Var = ctr->muExp(m2);   
  mixVar = 0;             

  // muMix for interaction
  if ((ctr->interaction) && ((ctr->interaction == 2) || (m1 != m2))) { 
    if (m1 <= m2)
      mixVar = ctr->muMix(m2, m1);
    else
      mixVar = ctr->muMix(m1, m2);
  }

  // Current null state MHR
  mhr0 = dlmtreeHDLMM_MHR(modTerm, dlmTerm1, dlmTerm2, 
                          ctr, ZtR, treeVar, 
                          m1Var, m2Var, mixVar);

  // // [Tree pair update]
  // // *** Propose a new TDLMM tree 1 ***
  // newExp    = m1; 
  // newExpVar = m1Var;
  // newMixVar = mixVar;
  // stepMhr = 0;
  // success = 0;

  // // *** Propose a new dlmtree 1 ***
  // // Choose a transition and update
  // step1 = sampleInt(ctr->stepProb, 1);
  // if ((dlmTerm1.size() == 1) && (step1 < 3)){
  //   step1 = 0;
  // }

  // if(step1 < 3){ 
  //   stepMhr = tdlmProposeTree(dlmTree1, Exp[m1], ctr, step1);     // Given current tree -> propose a new tree
  //   success = dlmTree1->isProposed();                             // Update to proposed
  //   newDlmTerm1 = dlmTree1->listTerminal(1);                      // If proposed, list the terminal nodes of the tree
  // } else {
  //   newExp = sampleInt(ctr->expProb);         // Sample an exposure
  //   if (newExp != m1) {                       // If new Exposure is not the same as the previous one,
  //     success   = 1;                          // It is good
  //     newExpVar = ctr->muExp(newExp);         // Find the exposure-specific variance for the new exposure
  //     newTree   = new Node(*dlmTree1);        // We need a new tree for a new exposure
  //     newTree->setUpdate(1);                  // Set an update flag + for children as well (This is used for "update NodeVals")
  //     newDlmTerm1 = newTree->listTerminal();  // List the terminal node of this new tree

  //     for (Node* nt : newDlmTerm1)            // Go through the new terminal node
  //       Exp[newExp]->updateNodeVals(nt);      // Update the math terms

  //     // Update the interaction using the new exposure as well for TDLMMns/all
  //     if ((ctr->interaction) && ((ctr->interaction == 2) || (newExp != m2))) {
  //       if (newExp <= m2)
  //         newMixVar = ctr->muMix(m2, newExp); 
  //       else
  //         newMixVar = ctr->muMix(newExp, m2); 
  //     } else {
  //       newMixVar = 0;
  //     }
  //   }
  // }

  // // * dlmTree 1 MHR 
  // // StepMHR: Transition ratio
  // // mhr0 & mhr: Likelihood of R
  // // If successful, calculate MH ratio
  // if (success) {  
  //   // If transition is successful, we have a new dlm terminal node.
  //   newDlmTerm1 = dlmTree1->listTerminal(1);
  //   modTree->setUpdateXmat(1);

  //   // MH ratio with a new terminal and the exposure: newDlmTerm1, newExpVar
  //   mhr = dlmtreeHDLMM_MHR(modTerm, 
  //                          newDlmTerm1, dlmTerm2, ctr, ZtR, treeVar, 
  //                          newExpVar, m2Var, newMixVar);

  //   // MH ratio - dlmTree 
  //   if (RtR < 0) {
  //     RtR = (ctr->R).dot(ctr->R);
  //     RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);
  //   }

  //   // MH ratio
  //   ratio = stepMhr + 
  //           mhr.logVThetaChol - mhr0.logVThetaChol -
  //           (0.5 * (ctr->n + 1.0) *
  //           (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
  //           log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
  //           (0.5 * ((mhr.nTerm1 * mhr0.nModTerm * log(treeVar * newExpVar)) -
  //           (mhr0.nTerm1 * mhr0.nModTerm * log(treeVar * m1Var)))); 

  //   // Interaction
  //   if (newMixVar != 0){ 
  //     ratio -= 0.5 * log(treeVar * newMixVar) * mhr.nTerm1 * mhr0.nTerm2 * mhr0.nModTerm;
  //   }

  //   if (mixVar != 0){ 
  //     ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2 * mhr0.nModTerm;
  //   }

  //   // Accept / Reject
  //   if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
  //     mhr0 = mhr;
  //     success = 2;

  //     if (step1 == 3) { 
  //       m1      = newExp;
  //       m1Var   = newExpVar;
  //       mixVar  = newMixVar;
  //       dlmTree1->replaceNodeVals(newTree);
  //     } else {
  //       dlmTree1->accept(); 
  //     }
  //     // Update parameters with new tree
  //     dlmTerm1 = dlmTree1->listTerminal();
  //     for (Node* n : modTerm) {
  //       n->nodevals->updateXmat = 0;
  //     }

  //   } else {  
  //     modTree->setUpdateXmat(1); 
  //   }
  // }
  // dlmTree1->reject();

  // // Reset new tree value for tree2
  // if (newTree != 0)
  //   delete newTree;
  // newTree = 0;

  // // * Record tree 1
  // if (ctr->diagnostics) {
  //   Eigen::VectorXd acc(9);
  //   acc << ctr->record, t, 1, step1, success, m1, dlmTerm1.size(), stepMhr, ratio;
  //   (dgn->treeDLMAccept).push_back(acc);
  // }

  
  // // *** Propose a new TDLMM tree 2 ***
  // newExp    = m2; 
  // newExpVar = m2Var;
  // newMixVar = mixVar;
  // stepMhr = 0;
  // success = 0;

  // // Choose a transition and update
  // step2 = sampleInt(ctr->stepProb, 1);
  // if ((dlmTerm2.size() == 1) && (step2 < 3)){
  //   step2 = 0;
  // }

  // if(step2 < 3){ 
  //   stepMhr = tdlmProposeTree(dlmTree2, Exp[m2], ctr, step2);
  //   success = dlmTree2->isProposed();
  //   newDlmTerm2 = dlmTree2->listTerminal(1);

  // } else { 
  //   newExp = sampleInt(ctr->expProb); 
  //   if (newExp != m2) { 
  //     success   = 1;
  //     newExpVar = ctr->muExp(newExp);    
  //     newTree   = new Node(*dlmTree2);   
  //     newTree->setUpdate(1);
  //     newDlmTerm2   = newTree->listTerminal();

  //     for (Node* nt : newDlmTerm2)
  //       Exp[newExp]->updateNodeVals(nt);

  //     if ((ctr->interaction) && ((ctr->interaction == 2) || (newExp != m1))) {
  //       if (newExp <= m1)
  //         newMixVar = ctr->muMix(m1, newExp);
  //       else
  //         newMixVar = ctr->muMix(newExp, m1);
  //     } else {
  //       newMixVar = 0;
  //     }
  //   }
  // }

  // // * dlmTree 2 MHR 
  // // StepMHR: Transition ratio
  // // mhr0 & mhr: Likelihood of R
  // // If successful, calculate MH ratio
  // if (success) {
  //   newDlmTerm2 = dlmTree2->listTerminal(1);
  //   modTree->setUpdateXmat(1);

  //   // MH ratio with a new terminal and the exposure: newDlmTerm1, newExpVar
  //   mhr = dlmtreeHDLMM_MHR(modTerm, 
  //                          dlmTerm1, newDlmTerm2, ctr, ZtR, treeVar, 
  //                          m1Var, newExpVar, newMixVar);

  //   // MH ratio - dlmTree
  //   if (RtR < 0) {
  //     RtR = (ctr->R).dot(ctr->R);
  //     RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);
  //   }

  //   ratio = stepMhr + 
  //           mhr.logVThetaChol - mhr0.logVThetaChol -
  //           (0.5 * (ctr->n + 1.0) *
  //           (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
  //           log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
  //           (0.5 * ((mhr.nTerm2 * mhr0.nModTerm * log(treeVar * newExpVar)) -
  //           (mhr0.nTerm2 * mhr0.nModTerm * log(treeVar * m2Var))));

  //   // Interaction
  //   if (newMixVar != 0){ 
  //     ratio -= 0.5 * log(treeVar * newMixVar) * mhr0.nTerm1 * mhr.nTerm2 * mhr0.nModTerm;
  //   }

  //   if (mixVar != 0){ 
  //     ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2 * mhr0.nModTerm;
  //   }

  //   if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
  //     mhr0 = mhr;
  //     success = 2;

  //     if (step2 == 3) { 
  //       m2      = newExp;
  //       m2Var   = newExpVar;
  //       mixVar  = newMixVar;
  //       dlmTree2->replaceNodeVals(newTree);
  //     } else {
  //       dlmTree2->accept();
  //     }

  //     dlmTerm2 = dlmTree2->listTerminal(); 
  //     // No need to update Xmat for dlmtree accept
  //     for (Node* n : modTerm) {
  //       n->nodevals->updateXmat = 0;
  //     }

  //   } else { 
  //     modTree->setUpdateXmat(1); 
  //   }
  // } 
  // dlmTree2->reject();  

  // if (newTree != 0)
  //   delete newTree;
  // newTree = 0;

  // // * Record tree 2
  // if (ctr->diagnostics) {
  //   Eigen::VectorXd acc(9);
  //   acc << ctr->record, t, 2, step2, success, m2, dlmTerm2.size(), stepMhr, ratio;
  //   (dgn->treeDLMAccept).push_back(acc);
  // }


  // [New tree proposal method] ------------------------------------------------------------------------------------------------
  // *** Propose a new TDLMM tree 1 ***
  newExp    = m1; 
  newExpVar = m1Var;
  newMixVar = mixVar;
  stepMhr = 0;
  success = 1;
  step1 = 0;

  // *** Propose a new dlmtree 1 ***
  newExp = sampleInt(ctr->expProb);       // Sample an exposure for the new tree
  newExpVar = ctr->muExp(newExp);         // Find the exposure-specific variance for the new exposure

  // Update the interaction using the new exposure as well for HDLMMns/all
  if ((ctr->interaction) && ((ctr->interaction == 2) || (newExp != m2))) {
    if (newExp <= m2)
      newMixVar = ctr->muMix(m2, newExp); 
    else
      newMixVar = ctr->muMix(newExp, m2); 
  } else {
    newMixVar = 0;
  }

  // *** Create a new tree for proposal ***
  newTree = new Node(0, 1);               // Start from the root
  newTree->nodestruct = expNS->clone();   // Construct nodestruct
  drawTree(newTree, newTree, ctr->treePrior[0], ctr->treePrior[1]); // Grow a tree structure from the root
  newTree->setUpdate(1);                  // Set a flag for updateNodeVals
  newDlmTerm1 = newTree->listTerminal();  // List the number of terminal nodes for the new tree
  for (Node* nt : newDlmTerm1) {          // Go through the new terminal node
    Exp[newExp]->updateNodeVals(nt);      // Update the calculations
  }
  
  // MH ratio
  modTree->setUpdateXmat(1);

  // MH ratio with a new terminal and the exposure: newDlmTerm1, newExpVar
  mhr = dlmtreeHDLMM_MHR(modTerm, 
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
          (0.5 * (((mhr.nTerm1 + mhr.nTerm2) * mhr0.nModTerm * log(treeVar * newExpVar)) -
          ((mhr0.nTerm1 + mhr0.nTerm2) * mhr0.nModTerm * log(treeVar * m1Var)))); 

  // Interaction
  if (newMixVar != 0){ 
    ratio -= 0.5 * log(treeVar * newMixVar) * mhr.nTerm1 * mhr.nTerm2 * mhr.nModTerm;
  }

  if (mixVar != 0){ 
    ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2 * mhr0.nModTerm;
  }

  // Accept / Reject
  if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
    // Rcout << "Accepted \n";

    mhr0 = mhr;
    success = 2;
    
    // Update exposures
    m1      = newExp;
    m1Var   = newExpVar;
    mixVar  = newMixVar;

    // Replace with new tree
    dlmTree1->replaceTree(newTree); 
    dlmTerm1 = dlmTree1->listTerminal();

    for (Node* n : modTerm) {
      n->nodevals->updateXmat = 0;
    }

  } else {  
    modTree->setUpdateXmat(1); 
    // Rcout << "Rejected \n";
    dlmTree1->reject(); 
  }
 
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


  // [New tree proposal method] 
  // *** Propose a new TDLMM tree 2 ***
  newExp    = m2; 
  newExpVar = m2Var;
  newMixVar = mixVar;
  stepMhr = 0;
  success = 1;

  // *** Create a new tree for proposal ***
  newTree = new Node(0, 1);               // Start from the root
  newTree->nodestruct = expNS->clone();   // Construct nodestruct
  drawTree(newTree, newTree, ctr->treePrior[0], ctr->treePrior[1]); // Grow a tree structure from the root
  newTree->setUpdate(1);                  // Update
  newDlmTerm2 = newTree->listTerminal();  // List the number of terminal nodes for the new tree
  for (Node* nt : newDlmTerm2) {          // Go through the new terminal node
    Exp[newExp]->updateNodeVals(nt);      // Update the calculations
  }

  // *** Propose a new dlmtree 1 ***
  step2 = 0;                              // Always propose
  newExp = sampleInt(ctr->expProb);       // Sample an exposure for the new tree
  newExpVar = ctr->muExp(newExp);         // Find the exposure-specific variance for the new exposure



  // HDLMMns
  // Update the interaction using the new exposure as well for TDLMMns/all
  if ((ctr->interaction) && ((ctr->interaction == 2) || (newExp != m1))) {
    if (newExp <= m1)
      newMixVar = ctr->muMix(m1, newExp);
    else
      newMixVar = ctr->muMix(newExp, m1);
  } else {
    newMixVar = 0;
  }
  
  // MH ratio
  modTree->setUpdateXmat(1);

  // MH ratio with a new terminal and the exposure: newDlmTerm1, newExpVar
  mhr = dlmtreeHDLMM_MHR(modTerm, 
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
          (0.5 * (((mhr.nTerm1 + mhr.nTerm2) * mhr0.nModTerm * log(treeVar * newExpVar)) -
          ((mhr0.nTerm1 + mhr0.nTerm2) * mhr0.nModTerm * log(treeVar * m2Var))));

  // Interaction
  if (newMixVar != 0){ 
    ratio -= 0.5 * log(treeVar * newMixVar) * mhr.nTerm1 * mhr.nTerm2 * mhr.nModTerm;
  }

  if (mixVar != 0){ 
    ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2 * mhr0.nModTerm;
  }

  // Accept / Reject
  if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
    mhr0 = mhr;
    success = 2;
    
    // Update exposures
    m2      = newExp;
    m2Var   = newExpVar;
    mixVar  = newMixVar;
    
    // Update parameters with new tree
    dlmTree2->replaceTree(newTree); 
    dlmTerm2 = dlmTree2->listTerminal(); 

    for (Node* n : modTerm) {
      n->nodevals->updateXmat = 0;
    }

  } else {  
    modTree->setUpdateXmat(1); 
    dlmTree2->reject(); 
  }
  
  // Reset new tree value for tree2
  if (newTree != 0)
    delete newTree;
  newTree = 0;

  // * Record tree 2
  if (ctr->diagnostics) {
    Eigen::VectorXd acc(9);
    acc << ctr->record, t, 2, step2, success, m2, dlmTerm2.size(), stepMhr, ratio;
    (dgn->treeDLMAccept).push_back(acc);
  }


  // *** Propose new modifier tree ***
  switch (modTerm.size()) {
    case 1: step = 0;   break;
    case 2: step = sampleInt(ctr->stepProbMod, 1 - ctr->stepProbMod[3]); break; // sampleInt(A vector of probability, total probability) // Fourth step is removed
    default: step = sampleInt(ctr->stepProbMod, 1);
  }

  stepMhr = modProposeTree(modTree, Mod, ctr, step);
  success = modTree->isProposed();

  if (success && (stepMhr == stepMhr)) {
    newModTerm = modTree->listTerminal(1);

    mhr = dlmtreeHDLMM_MHR(newModTerm, dlmTerm1, dlmTerm2, 
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

  // *** Update parameters ***
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

  // Update Rmat
  ctr->Rmat.col(t) = mhr0.fitted;

  // Update exposure information
  (ctr->dlmTree1Exp)(t) = m1;
  (ctr->dlmTree2Exp)(t) = m2;
  (ctr->expCount)(m1)++;                 
  (ctr->expCount)(m2)++;                 
  (ctr->expInf)(m1) += ((ctr->tau)(t)); 
  (ctr->expInf)(m2) += ((ctr->tau)(t)); 
  (ctr->totTermExp)(m1) += mhr0.nModTerm * mhr0.nTerm1;   
  (ctr->totTermExp)(m2) += mhr0.nModTerm * mhr0.nTerm2;   
  (ctr->sumTermT2Exp)(m1) += mhr0.term1T2 / (ctr->tau)(t); 
  (ctr->sumTermT2Exp)(m2) += mhr0.term2T2 / (ctr->tau)(t); 

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
  // *** Count modifiers used in tree ***
  Eigen::VectorXd modCount = countMods(modTree, Mod);
  ctr->modCount += modCount;

  // *** Record *** 
  // Rcout << "TreeMCMC: Record ... \n";
  if (ctr->record > 0) {    
    for (int i = 0; i < modCount.size(); i++) {
      if (modCount(i) > 0){
        ctr->modInf(i) += (ctr->tau)(t);
      }
    }

    // *** Update DLM partial estimate ***
    // Rcout << "Creating a rule ... \n";
    std::string rule; 
    Eigen::VectorXd rec(10); 
    Eigen::VectorXd mix(10);
    rec << ctr->record, t, 0, 0, 0, 0, 0, 0, 0, 0;
    mix << ctr->record, t, 0, 0, 0, 0, 0, 0, 0, 0;

    // For each terminal node of the modifier trees `s`,
    for (s = 0; s < modTerm.size(); s++) {
      rule = modRuleStr(modTerm[s], Mod); 

      // Iterate through dlmTree 1 and the mixture
      std::size_t k = 0; 
      for (std::size_t i = 0; i < mhr0.nTerm1; i++) {
        rec[2] = s;                                  
        rec[3] = 1;                                   
        rec[4] = i;                                   
        rec[5] = m1;                                  
        rec[6] = (dlmTerm1[i]->nodestruct)->get(3);   
        rec[7] = (dlmTerm1[i]->nodestruct)->get(4);   
        rec[8] = mhr0.drawAll[s * mhr0.nDlmTerm + i]; 
        rec[9] = (ctr->tau)(t) * m1Var;               

        (dgn->termRule).push_back(rule); 
        (dgn->DLMexp).push_back(rec);    
        
        for (std::size_t j = 0; j < mhr0.nTerm2; j++) {
          if (i == 0) {     
            rec[2] = s;   
            rec[3] = 2;                                   
            rec[4] = j;                                   
            rec[5] = m2;                                  
            rec[6] = (dlmTerm2[j]->nodestruct)->get(3);   
            rec[7] = (dlmTerm2[j]->nodestruct)->get(4);   
            rec[8] = mhr0.drawAll[s * mhr0.nDlmTerm + mhr0.nTerm1 + j];
            rec[9] = (ctr->tau)(t) * m2Var;
            
            (dgn->termRule).push_back(rule);
            (dgn->DLMexp).push_back(rec); 
          }

          if (mixVar != 0) {
            mix[2] = s; 
            if (m1 <= m2) {
              mix[3] = m1;
              mix[4] = (dlmTerm1[i]->nodestruct)->get(3); 
              mix[5] = (dlmTerm1[i]->nodestruct)->get(4); 
              mix[6] = m2;
              mix[7] = (dlmTerm2[j]->nodestruct)->get(3); 
              mix[8] = (dlmTerm2[j]->nodestruct)->get(4); 
            } else {
              mix[6] = m1;
              mix[7] = (dlmTerm1[i]->nodestruct)->get(3); 
              mix[8] = (dlmTerm1[i]->nodestruct)->get(4); 
              mix[3] = m2;
              mix[4] = (dlmTerm2[j]->nodestruct)->get(3); 
              mix[5] = (dlmTerm2[j]->nodestruct)->get(4);
            }
            mix[9] = mhr0.drawAll[s * mhr0.nDlmTerm + mhr0.nTerm1 + mhr0.nTerm2 + k];
            
            (dgn->termRuleMIX).push_back(rule); 
            (dgn->MIXexp).push_back(mix);       
            k++;
          } // mixVar != 0 end
        } // Tree 2 iteration end
      } // Tree 1 iteration end
    } // modifier tree iteration end
  } // ctr->record > 0 end
  // Rcout << "TreeMCMC: Record finished ... \n";

} // end dlmtreeHDLMMGaussian_TreeMCMC function


treeMHR dlmtreeHDLMM_MHR(std::vector<Node*> modTerm,  
                          std::vector<Node*> dlmTerm1,
                          std::vector<Node*> dlmTerm2,
                          dlmtreeCtr* ctr, 
                          Eigen::VectorXd ZtR, 
                          double treeVar, double m1Var, double m2Var, double mixVar)

{
  std::size_t s; 
  treeMHR out; 

  // exposure variance
  out.m1Var = m1Var;
  out.m2Var = m2Var;

  // modifier tree
  int pXMod = modTerm.size();

  // Number of terminal nodes
  int pXDlm1 = dlmTerm1.size();
  int pXDlm2 = dlmTerm2.size();
  int pXDlm = pXDlm1 + pXDlm2;
  int interaction = 0;
  if (mixVar != 0) { 
    pXDlm += pXDlm1 * pXDlm2; 
    interaction = 1; 
  }

  // Rcout << "Counting the totoal number of dlmtree x modifier terminal nodes \n";
  int pXComb = pXMod * pXDlm;   

  // Rcout << "Define parameters for calculation \n";
  Eigen::MatrixXd Xd, ZtX, VgZtX;
  Eigen::VectorXd diagVar;
  Xd.resize(ctr->n, pXDlm);         Xd.setZero(); 
  ZtX.resize(ctr->pZ, pXComb);      ZtX.setZero();
  VgZtX.resize(ctr->pZ, pXComb);    VgZtX.setZero();
  diagVar.resize(pXDlm);            diagVar.setZero(); 
  
  // Building exposure matrix, Xd & ZtX
  int i, j, k;

  // dlmTree 1
  for (i = 0; i < pXDlm1; i++) {               
    Xd.col(i) = (dlmTerm1[i]->nodevals)->X;
    diagVar(i) = 1.0 / (m1Var * treeVar);
    ZtX.col(i) = (dlmTerm1[i]->nodevals)->ZtX;
  }

  // dlmTree 2
  for (j = 0; j < pXDlm2; j++) {               
    k = pXDlm1 + j;                            
    Xd.col(k) = (dlmTerm2[j]->nodevals)->X; 
    diagVar(k) = 1.0 / (m2Var * treeVar);      
    ZtX.col(k) = (dlmTerm2[j]->nodevals)->ZtX;
  }

  // Mixture & interaction
  if(interaction) {
    for (i = 0; i < pXDlm1; i++) {
      for (j = 0; j < pXDlm2; j++) {
        k = (pXDlm1 + pXDlm2) + i * pXDlm2 + j; 
        Xd.col(k) = (((dlmTerm1[i]->nodevals)->X).array() * ((dlmTerm2[j]->nodevals)->X).array()).matrix();
        diagVar(k) = 1.0 / (mixVar * treeVar);
        ZtX.col(k) = ctr->Zw.transpose() * Xd.col(k); 
      }
    }
  }

  // MH ratio element calculation
  if(pXMod == 1){   
    const Eigen::MatrixXd VgZtX = ctr->Vg * ZtX; 
    Eigen::MatrixXd tempV(pXDlm, pXDlm);        
    Eigen::VectorXd XtVzInvR; 

    tempV.triangularView<Eigen::Lower>() = Xd.transpose() * Xd;  
    tempV.noalias() -= ZtX.transpose() * VgZtX; 
    out.tempV = tempV;
    
    tempV.diagonal().noalias() += diagVar; 
    XtVzInvR = Xd.transpose() * ctr->R - VgZtX.transpose() * ZtR;

    // Calculate Vtheta & Vtheta cholesky
    Eigen::MatrixXd VTheta(pXDlm, pXDlm); 
    VTheta.triangularView<Eigen::Lower>() =
      tempV.selfadjointView<Eigen::Lower>().llt().solve(Eigen::MatrixXd::Identity(pXDlm, pXDlm));
    const Eigen::MatrixXd VThetaChol = VTheta.selfadjointView<Eigen::Lower>().llt().matrixL();

    // Thetahat mean
    const Eigen::VectorXd ThetaHat = VTheta.selfadjointView<Eigen::Lower>() * XtVzInvR; 

    // Thetahat variance
    Eigen::VectorXd ThetaDraw = ThetaHat;
    ThetaDraw.noalias() += VThetaChol * as<Eigen::VectorXd>(rnorm(pXDlm, 0, sqrt(ctr->sigma2)));

    // Full conditional draw
    out.drawAll = ThetaDraw;

    // Temporary vectors for calculation below
    Eigen::VectorXd drawTemp(pXDlm);    
    Eigen::VectorXd draw1Temp(pXDlm1);  
    Eigen::VectorXd draw2Temp(pXDlm2);  

    // Exposure 1 return
    // Rcout << "pxMod == 1: Exposure 1 segment \n";
    draw1Temp = ThetaDraw.head(pXDlm1);            
    out.term1T2 = (draw1Temp).dot(draw1Temp);      
    out.nTerm1 = double(pXDlm1);                   

    // Exposure 2 return
    // Rcout << "pxMod == 1: Exposure 2 segment \n";
    draw2Temp = ThetaDraw.segment(pXDlm1, pXDlm2); 
    out.term2T2 = (draw2Temp).dot(draw2Temp);      
    out.nTerm2 = double(pXDlm2);                   

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

  // *** Multiple Modifier nodes ***
  Eigen::MatrixXd Xtemp, Ztemp;
  Eigen::VectorXd Rtemp;
  Eigen::MatrixXd XtXblock(pXComb, pXComb);    XtXblock.setZero();
  Eigen::VectorXd XtR(pXComb);                 XtR.setZero();
  Eigen::MatrixXd LInv(pXDlm, pXDlm);          LInv.setZero();  
  LInv.diagonal().array() += diagVar.array();

  // Create block matrices corresponding to modifier nodes
  int start = 0;
  for (Node* n : modTerm) { 
    // Retrieve indices subset by the modifier tree
    Xtemp.resize(n->nodevals->idx.size(), pXDlm);     Xtemp.setZero();  
    Ztemp.resize(n->nodevals->idx.size(), ctr->pZ);   Ztemp.setZero();  
    Rtemp.resize(n->nodevals->idx.size());            Rtemp.setZero();  
    
    // Loop through indices and (subset) update Xtemp, Ztemp, Rtemp
    j = 0;
    for (int i : n->nodevals->idx) {
      Xtemp.row(j) = Xd.row(i);  
      Ztemp.row(j) = ctr->Z.row(i); 
      Rtemp(j) = ctr->R(i);         
      j++;
    } // end loop over node indices
      
    n->nodevals->XtX.resize(pXDlm, pXDlm);
    n->nodevals->XtX = Xtemp.transpose() * Xtemp;
    n->nodevals->ZtXmat.resize(ctr->pZ, pXDlm);
    n->nodevals->ZtXmat = Ztemp.transpose() * Xtemp;
    n->nodevals->VgZtXmat.resize(ctr->pZ, pXDlm);
    n->nodevals->VgZtXmat = ctr->Vg * n->nodevals->ZtXmat;
    n->nodevals->updateXmat = 0;
    
    XtR.segment(start, pXDlm) = Xtemp.transpose() * Rtemp;
      
    // Update blocks
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

  // Store the sampled values and also add variance 
  out.drawAll = ThetaHat;
  out.drawAll.noalias() += 
      VThetaChol * as<Eigen::VectorXd>(rnorm(pXComb, 0.0, sqrt(ctr->sigma2)));

  // drawAll = mod1-dlm1 / mod1-dlm2 / mod1-dlm1&2 / mod2-dlm1 / mod2-dlm2 / mod2-dlm1&2 / mod3 ...
  // Fitted value & terminal effect draws
  out.fitted.resize(ctr->n);                 
  Eigen::VectorXd drawTemp(pXDlm);         
  Eigen::VectorXd draw1Temp(pXDlm1);     
  Eigen::VectorXd draw2Temp(pXDlm2);     
  Eigen::VectorXd drawMixTemp(pXDlm1 * pXDlm2);
  
  if(interaction){  
    drawMixTemp.resize(pXDlm - pXDlm1 - pXDlm2);
  }

  // Store Term1T2/Term2T2/MixT2 & nTerm1/nTerm2
  out.term1T2 = 0;
  out.term2T2 = 0;
  out.mixT2 = 0;
  
  for (s = 0; s < modTerm.size(); s++) { 
    // Extract draws from the block matrix for each modifier terminal node
    drawTemp = out.drawAll.segment(s * pXDlm, pXDlm); 

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
      out.nTermMix = double(pXDlm - pXDlm1 - pXDlm2);
    }

    // Calculate fitted for each modifier terminal
    for (int i : modTerm[s]->nodevals->idx){
      out.fitted(i) = Xd.row(i) * drawTemp;
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
