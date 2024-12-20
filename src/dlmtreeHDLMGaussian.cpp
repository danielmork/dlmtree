#include "RcppEigen.h"
#include "Node.h"
#include "NodeStruct.h"
#include "modDat.h"
#include "exposureDat.h"
#include "Fncs.h"
#include "modelCtr.h"
using namespace Rcpp;


void dlmtreeHDLMGaussian_TreeMCMC(int t, NodeStruct* expNS, NodeStruct* modNS, 
                                  Node* modTree, Node* dlmTree,
                                  dlmtreeCtr* ctr, dlmtreeLog *dgn,
                                  modDat* Mod, exposureDat* Exp);

treeMHR dlmtreeTDLM_MHR(std::vector<Node*> modTerm,
                        std::vector<Node*> dlmTerm,
                        dlmtreeCtr* ctr, 
                        Eigen::VectorXd ZtR, 
                        double treevar);

double calcLogRatioHDLM(treeMHR mhr0, treeMHR mhr, double RtR, double RtZVgZtR,
                        dlmtreeCtr* ctr, double stepMhr, double treevar)
{
  if (ctr->binomial) {
    return(stepMhr + mhr.logVThetaChol - mhr0.logVThetaChol +
           0.5 * (mhr.beta - mhr0.beta) -
           log(treevar) * 0.5 * round(mhr.totTerm - mhr0.totTerm));
  } else {
    return(stepMhr + mhr.logVThetaChol - mhr0.logVThetaChol -
          (0.5 * (ctr->n + 1.0) *
            (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
            log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
          (log(treevar) * 0.5 * round(mhr.totTerm - mhr0.totTerm)));
  }
}


//' dlmtree model with shared HDLM approach
//'
//' @param model A list of parameter and data contained for the model fitting
//' @returns A list of dlmtree model fit, mainly posterior mcmc samples
//' @export
// [[Rcpp::export]]
Rcpp::List dlmtreeHDLMGaussian(const Rcpp::List model)
{
  int t;
  // ---- Set up general control variables ----
  dlmtreeCtr *ctr = new dlmtreeCtr;
  ctr->iter       = as<int>(model["nIter"]);
  ctr->burn       = as<int>(model["nBurn"]);
  ctr->thin       = as<int>(model["nThin"]);
  ctr->nRec       = floor(ctr->iter / ctr->thin);
  ctr->nTrees     = as<int>(model["nTrees"]);
  // ctr->Y      = as<VectorXd>(model["Y"]);
  ctr->Y0     = as<Eigen::VectorXd>(model["Y"]);      
  ctr->Ystar  = as<Eigen::VectorXd>(model["Y"]);  
  ctr->n          = (ctr->Y0).size();
  ctr->modZeta    = as<double>(model["zeta"]);
  ctr->modKappa   = 100;

  ctr->Z      = as<Eigen::MatrixXd>(model["Z"]);
  ctr->Zw     = ctr->Z;
  ctr->pZ     = (ctr->Z).cols();
  ctr->VgInv  = (ctr->Z).transpose() * (ctr->Z);
  ctr->VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg     = ctr->VgInv.inverse();
  ctr->VgChol = (ctr->Vg).llt().matrixL();
  
  ctr->binomial     = as<bool>(model["binomial"]);;
  ctr->zinb         = 0;
  ctr->verbose      = bool (model["verbose"]);
  ctr->diagnostics  = bool (model["diagnostics"]);
  ctr->stepProb     = as<std::vector<double> >(model["stepProbTDLM"]);
  ctr->stepProbMod  = as<std::vector<double> >(model["stepProbMod"]);
  ctr->treePrior    = as<std::vector<double> >(model["treePriorMod"]);
  ctr->treePriorMod = as<std::vector<double> >(model["treePriorTDLM"]);
  ctr->shrinkage    = as<int>(model["shrinkage"]);

  // * Set up data for logistic model
  ctr->binomialSize.resize(ctr->n);           ctr->binomialSize.setZero();
  ctr->kappa.resize(ctr->n);                  ctr->kappa.setOnes();
  ctr->Omega.resize(ctr->n);                  ctr->Omega.setOnes();
  // Rcout << ctr->binomial;
  if (ctr->binomial) {
    ctr->binomialSize = as<VectorXd>(model["binomialSize"]);
    ctr->kappa        = ctr->Y0 - 0.5 * (ctr->binomialSize);
    ctr->Ystar            = ctr->kappa;
  }

  // ---- Setup modifier data ----
  modDat *Mod = new modDat(as<std::vector<int> >(model["modIsNum"]),
                           as<Rcpp::List>(model["modSplitIdx"]),
                           as<std::vector<int> >(model["fullIdx"]));

  NodeStruct *modNS;
  modNS   = new ModStruct(Mod);
  ctr->pM = Mod->nMods;

  // ---- Pre-calculate single node tree matrices ----
  exposureDat *Exp;
  // if (as<int>(model["nSplits"]) == 0) { // DLM
    if (ctr->binomial)
      Exp = new exposureDat(as<MatrixXd>(model["Tcalc"]));
    else
      Exp = new exposureDat(as<MatrixXd>(model["Tcalc"]),
                            ctr->Z, ctr->Vg);
  // } else { // DLNM
  //   if (ctr->binomial)
  //     Exp = new exposureDat(as<MatrixXd>(model["X"]),
  //                           as<MatrixXd>(model["SE"]),
  //                           as<VectorXd>(model["Xsplits"]),
  //                           as<MatrixXd>(model["Xcalc"]),
  //                           as<MatrixXd>(model["Tcalc"]));
  //   else
  //     Exp = new exposureDat(as<MatrixXd>(model["X"]),
  //                           as<MatrixXd>(model["SE"]),
  //                           as<VectorXd>(model["Xsplits"]),
  //                           as<MatrixXd>(model["Xcalc"]),
  //                           as<MatrixXd>(model["Tcalc"]),
  //                           ctr->Z, ctr->Vg);
  // }
  ctr->pX         = Exp->pX;
  ctr->nSplits    = Exp->nSplits;

  ctr->X1         = (Exp->Tcalc).col(ctr->pX - 1);
  ctr->ZtX1       = (ctr->Z).transpose() * (ctr->X1);
  ctr->VgZtX1     = (ctr->Vg).selfadjointView<Eigen::Lower>() * (ctr->ZtX1);
  ctr->VTheta1Inv = (ctr->X1).dot(ctr->X1) - (ctr->ZtX1).dot(ctr->VgZtX1);

  // ---- Create trees ----
  NodeStruct *expNS;
  expNS = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                         as<Eigen::VectorXd>(model["splitProb"]),
                         as<Eigen::VectorXd>(model["timeProb"]));

  std::vector<Node*> modTrees;
  std::vector<Node*> dlmTrees;
  for (t = 0; t < ctr->nTrees; ++t) {
    modTrees.push_back(new Node(0, 1));
    modTrees[t]->nodestruct = modNS->clone();
    Mod->updateNodeVals(modTrees[t]);

    dlmTrees.push_back(new Node(0, 1));
    dlmTrees[t]->nodestruct = expNS->clone();
    Exp->updateNodeVals(dlmTrees[t]);
  }


  // ---- Logs ----
  dlmtreeLog *dgn = new dlmtreeLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec);    (dgn->gamma).setZero();
  (dgn->sigma2).resize(ctr->nRec);            (dgn->sigma2).setZero();
  (dgn->nu).resize(ctr->nRec);                (dgn->nu).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec);  (dgn->tau).setZero();
  (dgn->fhat).resize(ctr->n);                 (dgn->fhat).setZero();
  (dgn->totTerm).resize(ctr->nRec);           (dgn->totTerm).setZero();

  (dgn->modProb).resize(ctr->pM, ctr->nRec);  (dgn->modProb).setZero();
  (dgn->modCount).resize(ctr->pM, ctr->nRec); (dgn->modCount).setZero();
  (dgn->modInf).resize(ctr->pM, ctr->nRec);   (dgn->modInf).setZero();
  (dgn->modKappa).resize(ctr->nRec);          (dgn->modKappa).setZero();

  (dgn->termNodesMod).resize(ctr->nTrees, ctr->nRec);   (dgn->termNodesMod).setZero();
  (dgn->termNodesDLM).resize(ctr->nTrees, ctr->nRec);   (dgn->termNodesDLM).setZero();

  // ---- Initial draws ----
  (ctr->fhat).resize(ctr->n); (ctr->fhat).setZero();
  ctr->R = ctr->Ystar;
  (ctr->gamma).resize(ctr->pZ);

  // Load initial params for faster convergence in binomial model
  if (ctr->binomial) {
    // Rcout << "init";
    ctr->gamma  = as<VectorXd>(model["initParams"]);
    ctr->Omega  = rcpp_pgdraw(ctr->binomialSize, ctr->fhat + ctr->Z * ctr->gamma);
    ctr->Zw     = ctr->Omega.asDiagonal() * ctr->Z;
    ctr->VgInv  = ctr->Z.transpose() * ctr->Zw;
    ctr->VgInv.diagonal().array() += 1 / 100000.0;
    ctr->Vg     = ctr->VgInv.inverse();
    ctr->VgChol = ctr->Vg.llt().matrixL();
    // recalculate 'pseudo-Y' = kappa / omega, kappa = (y - n_b)/2
    ctr->Ystar      = ctr->kappa.array() / ctr->Omega.array();
    // Rcout << "end\n";
  }

  ctr->totTerm    = 0;
  ctr->sumTermT2  = 0;
  ctr->nu         = 1.0; // Need to define for first update of sigma2
  ctr->sigma2     = 1.0;
  tdlmModelEst(ctr);
  double xiInv    = R::rgamma(1, 0.5);
  ctr->nu = 1.0 / R::rgamma(0.5 * ctr->nTrees + 0.5, 1.0 / xiInv);

  (ctr->tau).resize(ctr->nTrees);     (ctr->tau).setOnes();
  if (ctr->shrinkage) {
    for (t = 0; t < ctr->nTrees; t++) {
      xiInv = R::rgamma(1, 0.5);
      (ctr->tau)(t) = 1.0 / R::rgamma(0.5, 1.0 / xiInv);
    }
  }
  ctr->nTerm.resize(ctr->nTrees);             ctr->nTermMod.resize(ctr->nTrees);
  (ctr->nTerm).array() = 1;                   ctr->nTermMod.array() = 1;
  (ctr->Rmat).resize(ctr->n, ctr->nTrees);    ctr->Rmat.setZero();
  ctr->modCount.resize(ctr->pM);              ctr->modCount.setZero();
  ctr->modInf.resize(ctr->pM);                ctr->modInf.setZero();
  // ctr->exDLM.resize(ctr->pX, ctr->n);
  
  // create progress meter
  progressMeter* prog = new progressMeter(ctr);

  std::size_t s;
  // ---- MCMC ----
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    Rcpp::checkUserInterrupt();
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }

    // -- Update trees --
    ctr->R += (ctr->Rmat).col(0);
    (ctr->fhat).setZero();
    ctr->totTerm = 0.0; ctr->sumTermT2 = 0.0;
    // ctr->exDLM.setZero();
    ctr->modCount.setZero();
    ctr->modInf.setZero();

    for (t = 0; t < ctr->nTrees; t++) {
      dlmtreeHDLMGaussian_TreeMCMC(t, expNS, modNS, modTrees[t], dlmTrees[t], ctr, dgn, Mod, Exp);
      ctr->fhat += (ctr->Rmat).col(t);

      if (t < ctr->nTrees - 1){
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
      }
    }

    // -- Update model --
    ctr->R = ctr->Ystar - ctr->fhat;    
    tdlmModelEst(ctr);
    xiInv   = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->nu)));
    ctr->nu = 1.0 / R::rgamma(0.5 * ctr->totTerm + 0.5,
                               1.0 / (0.5 * ctr->sumTermT2 / (ctr->sigma2) + xiInv));

    // -- Update modifier selection --
    if ((ctr->b > 1000) || (ctr->b > (0.5 * ctr->burn))) {
      double beta         = R::rbeta(ctr->modZeta, 1.0);
      double modKappaNew  = beta * ctr->pM / (1 - beta);
      double mhrDir =
        logDirichletDensity(Mod->modProb,
                            (ctr->modCount.array() + modKappaNew / ctr->pM).matrix()) -
        logDirichletDensity(Mod->modProb,
                            (ctr->modCount.array() + ctr->modKappa / ctr->pM).matrix());

      if (log(R::runif(0, 1)) < mhrDir) {
        ctr->modKappa = modKappaNew;
      }

      Mod->modProb = rDirichlet((ctr->modCount.array() + ctr->modKappa / ctr->pM).matrix());
    } // end modifier selection

    // -- Record --
    if (ctr->record > 0) {
      (dgn->gamma).col(ctr->record - 1)         = ctr->gamma;
      (dgn->sigma2)(ctr->record - 1)            = ctr->sigma2;
      (dgn->nu)(ctr->record - 1)                = ctr->nu;
      (dgn->tau).col(ctr->record - 1)           = ctr->tau;
      (dgn->termNodesDLM).col(ctr->record - 1)  = ctr->nTerm;
      (dgn->termNodesMod).col(ctr->record - 1)  = ctr->nTermMod;
      (dgn->modProb).col(ctr->record - 1)       = Mod->modProb;
      (dgn->modCount).col(ctr->record - 1)      = ctr->modCount;
      (dgn->modInf).col(ctr->record - 1)        = ctr->modInf / ctr->modInf.maxCoeff();
      (dgn->modKappa)(ctr->record - 1)          = ctr->modKappa;
      (dgn->totTerm)(ctr->record - 1)           = ctr->totTerm;
      dgn->fhat += ctr->fhat;

      // if (ctr->nSplits == 0)
      //   dlmtreeRecDLM(ctr, dgn);
    } // end record

    // -- Progress --
    prog->printMark();

    // } // end progress
  } // end MCMC


  // -- Prepare outout --
  Eigen::MatrixXd TreeStructs;
  Rcpp::StringVector termRule(dgn->termRule.size());
  TreeStructs.resize((dgn->DLMexp).size(), 9);
  for (s = 0; s < (dgn->DLMexp).size(); ++s){
    TreeStructs.row(s) = dgn->DLMexp[s];
  }
  termRule = dgn->termRule;
  
  Eigen::MatrixXd termNodesDLM  = (dgn->termNodesDLM).transpose();

  Eigen::VectorXd sigma2        = dgn->sigma2;
  Eigen::VectorXd nu            = dgn->nu;
  Eigen::MatrixXd tau           = (dgn->tau).transpose();
  Eigen::VectorXd fhat          = (dgn->fhat).array() / ctr->nRec;
  Eigen::VectorXd totTerm       = dgn->totTerm;  
  Eigen::MatrixXd gamma         = (dgn->gamma).transpose();

  Eigen::MatrixXd termNodesMod  = (dgn->termNodesMod).transpose();
  Eigen::VectorXd kappa         = dgn->modKappa;
  Eigen::MatrixXd modProb       = (dgn->modProb).transpose();
  Eigen::MatrixXd modCount      = (dgn->modCount).transpose();
  Eigen::MatrixXd modInf        = (dgn->modInf).transpose();

  Eigen::MatrixXd modAccept((dgn->treeModAccept).size(), 5);
  Eigen::MatrixXd dlmAccept((dgn->treeDLMAccept).size(), 5);

  delete prog;
  delete ctr;
  delete dgn;
  delete Exp;
  delete Mod;
  delete modNS;
  delete expNS;
  for (s = 0; s < modTrees.size(); ++s) {
    delete modTrees[s];
    delete dlmTrees[s];
  }

  return(Rcpp::List::create(// Named("DLM") = wrap(exDLM),
                            // Named("DLMse") = wrap(ex2DLM),
                            // Named("DLfun") = wrap(cumDLM),
                            // Named("DLfunse") = wrap(cum2DLM),
                            // Named("fhat") = wrap(fhat),
                            Named("TreeStructs")    = wrap(TreeStructs),
                            Named("termRules")      = wrap(termRule),
                            Named("termNodesDLM")   = wrap(termNodesDLM),
                            Named("totTerm")        = wrap(totTerm),
                            Named("sigma2")         = wrap(sigma2),
                            Named("nu")             = wrap(nu),
                            Named("tau")            = wrap(tau),
                            Named("gamma")          = wrap(gamma),
                            Named("termNodesMod")   = wrap(termNodesMod),
                            Named("kappa")          = wrap(kappa),
                            Named("modProb")        = wrap(modProb),
                            Named("modCount")       = wrap(modCount),
                            Named("modInf")         = wrap(modInf),
                            Named("treeModAccept")  = wrap(modAccept),
                            Named("treeDLMAccept")  = wrap(dlmAccept)));

} // end dlmtreeHDLMGaussian


void dlmtreeHDLMGaussian_TreeMCMC(int t, NodeStruct* expNS, NodeStruct* modNS, 
                                  Node* modTree, Node* dlmTree,
                                  dlmtreeCtr* ctr, dlmtreeLog *dgn,
                                  modDat* Mod, exposureDat* Exp)
{
  int step;
  int success     = 0;
  double stepMhr  = 0.0;
  double ratio    = 0.0;
  double treevar  = (ctr->nu) * (ctr->tau)(t);

  std::size_t s;
  std::vector<Node*> modTerm, dlmTerm, newDlmTerm, newModTerm;
  VectorXd ZtR      = ctr->Zw.transpose() * ctr->R;
  double RtR        = 0.0;
  double RtZVgZtR   = 0.0;
  if (!(ctr->binomial)) {
    RtR       = ctr->R.dot(ctr->R);
    RtZVgZtR  = ZtR.dot(ctr->Vg * ZtR);
  }

  treeMHR mhr0, mhr;

  // -- List terminal nodes --
  modTerm = modTree->listTerminal();
  dlmTerm = dlmTree->listTerminal();
  mhr0    = dlmtreeTDLM_MHR(modTerm, dlmTerm, ctr, ZtR, treevar);

  // -- Propose new TDLM tree --
  stepMhr = 0;
  success = 1;

  // *** Create a new tree for proposal ***
  Node* newTree = new Node(0, 1);               // Start from the root
  newTree->nodestruct = expNS->clone();   // Construct nodestruct
  drawTree(newTree, newTree, ctr->treePrior[0], ctr->treePrior[1]); // Grow a tree structure from the root
  newTree->setUpdate(1);                  // Set a flag for updateNodeVals
  newDlmTerm = newTree->listTerminal();  // List the number of terminal nodes for the new tree
  for (Node* nt : newDlmTerm) {          // Go through the new terminal node
    Exp->updateNodeVals(nt);      // Update the calculations
  }
  

  // if (success) {
  modTree->setUpdateXmat(1);
  mhr = dlmtreeTDLM_MHR(modTerm, newDlmTerm, ctr, ZtR, treevar);
  
  ratio = calcLogRatioHDLM(mhr0, mhr, RtR, RtZVgZtR, ctr, stepMhr, treevar);

  if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
    mhr0    = mhr;
    success = 2;

    // replace with new tree
    dlmTree->replaceTree(newTree);
    dlmTree->accept();
    dlmTerm = dlmTree->listTerminal();
    for (Node* n : modTerm) {
      n->nodevals->updateXmat = 0;
    }

  } else {
    modTree->setUpdateXmat(1);
  }
  // } // end dlmTree proposal
  // Reset new tree value 
  if (newTree != 0){
    delete newTree;
  }
  dlmTree->reject();
  newTree = 0;


  // -- Propose new modifier tree --
  // switch (modTerm.size()) {
  //   case 1: step  = 0; break;
  //   case 2: step  = sampleInt(ctr->stepProbMod, 1 - ctr->stepProbMod[3]); break;
  //   default: step = sampleInt(ctr->stepProbMod, 1);
  // }
  // stepMhr = modProposeTree(modTree, Mod, ctr, step);
  // success = modTree->isProposed();

  // if (success) {
  //   newModTerm = modTree->listTerminal(1);
  //   mhr   = dlmtreeTDLM_MHR(newModTerm, dlmTerm, ctr, ZtR, treevar);
  //   ratio = calcLogRatioHDLM(mhr0, mhr, RtR, RtZVgZtR, ctr, stepMhr, treevar);

  //   if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
  //     mhr0    = mhr;
  //     success = 2;
  //     modTree->accept();
  //     modTerm = modTree->listTerminal();
  //   } 
  // } // end modTree proposal
  // modTree->reject();


  // *** Create a new modifier tree for proposal ***
  stepMhr = 0;
  success = 1;
  newTree = new Node(0, 1);
  newTree->nodestruct = modNS->clone();   // Construct nodestruct
  drawTree(newTree, newTree, ctr->treePrior[0], ctr->treePrior[1]); // Grow a tree structure from the root
  newTree->setUpdate(1);                  // Set a flag for updateNodeVals
  newTree->setUpdateXmat(1);                  // Set a flag for updateNodeVals

  newModTerm = newTree->listTerminal(1);
  for (Node* tn : newModTerm) {
    Mod->updateNodeVals(tn);
    if (tn->nodevals->idx.size() == 0) {
      success = 0;
    } // end reject if empty
  }

  if (success) {
    mhr   = dlmtreeTDLM_MHR(newModTerm, dlmTerm, ctr, ZtR, treevar);
    ratio = calcLogRatioHDLM(mhr0, mhr, RtR, RtZVgZtR, ctr, stepMhr, treevar);

    if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
      mhr0    = mhr;
      success = 2;

      // replace with new tree
      modTree->replaceTree(newTree);
      modTree->accept();
      modTerm = modTree->listTerminal();
    } 
  }
  // Reset new tree value 
  if (newTree != 0)
    delete newTree;
  newTree = 0;


  // -- Update variance and residuals --
  if (ctr->shrinkage) {
    double xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->tau)(t)));
    (ctr->tau)(t) = 1.0 / R::rgamma(0.5 * mhr0.nDlmTerm * mhr0.nModTerm + 0.5,
                                    1.0 / ((0.5 * mhr0.termT2 / (ctr->sigma2 * ctr->nu)) + xiInv));
  }
  ctr->Rmat.col(t)  = mhr0.fitted;
  ctr->sumTermT2 += mhr0.termT2 / (ctr->tau(t));
  ctr->totTerm += mhr0.nDlmTerm * mhr0.nModTerm;
  ctr->nTermMod(t)  = mhr0.nModTerm;
  ctr->nTerm(t)     = mhr0.nDlmTerm;  
  
  // -- Count modifiers used in tree --
  Eigen::VectorXd modCount = countMods(modTree, Mod);
  ctr->modCount += modCount;

  // -- Record --
  if (ctr->record > 0) {    
    for (int i = 0; i < modCount.size(); ++i) {
      if (modCount(i) > 0){
        ctr->modInf(i) += ctr->tau(t);
      }
    }

    // -- Update DLM partial estimate --
    std::string rule;
    Eigen::VectorXd rec(9);
    // Eigen::VectorXd draw(ctr->pX);
    for (s = 0; s < modTerm.size(); ++s) {
      rule = modRuleStr(modTerm[s], Mod);
      for (std::size_t s2 = 0; s2 < dlmTerm.size(); ++s2) {
        rec << ctr->record, t, s, s2, 
          (dlmTerm[s2]->nodestruct)->get(1), // Exposure minimum value
          (dlmTerm[s2]->nodestruct)->get(2), // Exposure maximum value
          (dlmTerm[s2]->nodestruct)->get(3), // tmin
          (dlmTerm[s2]->nodestruct)->get(4), // tmax
          mhr0.draw[s * mhr0.nDlmTerm + s2];
        dgn->termRule.push_back(rule);
        dgn->DLMexp.push_back(rec);
        
        // if (ctr->nSplits == 0) {
        //   for (int t2 = dlmTerm[s2]->nodestruct->get(3) - 1; t2 < dlmTerm[s2]->nodestruct->get(4); ++t2){
        //     draw(t2) = mhr0.draw(s * mhr0.nDlmTerm + s2);
        //   }
        // }
      } // end loop over tdlm nodes
    } // end loop over modifier nodes to update partial DLM est
  } // end record


} // end dlmtreeHDLMGaussian_TreeMCMC function


treeMHR dlmtreeTDLM_MHR(std::vector<Node*> modTerm, 
                        std::vector<Node*> dlmTerm,
                        dlmtreeCtr* ctr, 
                        Eigen::VectorXd ZtR, 
                        double treevar)
// Calculate part of Metropolis-Hastings ratio and make draws from full
// conditional. 
{
  std::size_t s;
  treeMHR out;
  int pXDlm   = int(dlmTerm.size());
  int pXMod   = int(modTerm.size());
  int pXComb  = pXMod * pXDlm;

  Eigen::MatrixXd X, ZtX, VgZtX;  
  X.resize(ctr->n, pXDlm);          X.setZero();
  ZtX.resize(ctr->pZ, pXComb);      ZtX.setZero();
  VgZtX.resize(ctr->pZ, pXComb);    VgZtX.setZero();

  // Single TDLM node
  if (pXDlm == 1) {
    // Single Modifier node
    if (pXMod == 1) {
      double XtVzInvR = 0.0;
      double VTheta = 0.0;
      double tempV = 0.0;
      if (ctr->binomial) {
        const Eigen::VectorXd Xdw = (ctr->Omega).asDiagonal() * ctr->X1;  
        ZtX             = ctr->Zw.transpose() * ctr->X1;
        VgZtX           = ctr->Vg * ZtX;
        tempV           = ctr->X1.dot(Xdw);
        tempV          -= ZtX.col(0).dot(VgZtX.col(0));
        XtVzInvR        = Xdw.dot(ctr->R);
        
      } else {
        ZtX             = ctr->ZtX1;
        VgZtX           = ctr->VgZtX1;
        tempV           = ctr->X1.dot(ctr->X1);
        tempV          -= ctr->ZtX1.dot(ctr->VgZtX1);
        XtVzInvR        = ctr->X1.dot(ctr->R);
      }
      XtVzInvR         -= (VgZtX.col(0)).dot(ZtR);
      tempV            += 1.0 / treevar;
      VTheta            = 1.0 / tempV;
      double VThetaChol = sqrt(VTheta);
      double ThetaHat   = VTheta * XtVzInvR;
      out.draw.resize(1);
      out.draw(0)       = VThetaChol * R::rnorm(0, sqrt(ctr->sigma2)) + ThetaHat;
      out.fitted        = ctr->X1 * out.draw(0);
      out.beta          = ThetaHat * XtVzInvR;
      out.logVThetaChol = log(VThetaChol);
      out.termT2        = (out.draw).dot(out.draw);
      out.nDlmTerm      = 1.0; 
      out.nModTerm      = 1.0;
      out.totTerm       = 1.0;

      return(out);
    } // return single TDLM and single modifier node

  // Multiple TDLM nodes
  } else {

    for (s = 0; s < dlmTerm.size(); ++s) { // loop over dlmTerm
      X.col(s) = dlmTerm[s]->nodevals->X;
      if (pXMod == 1) {
        if (ctr->binomial) {
          ZtX.col(s)    = ctr->Zw.transpose() * X.col(s);
          VgZtX.col(s)  = ctr->Vg * ZtX.col(s);
          
        } else {
          ZtX.col(s)    = dlmTerm[s]->nodevals->ZtX;
          VgZtX.col(s)  = dlmTerm[s]->nodevals->VgZtX;
        }
      }
    } // end loop over dlmTerm

    // Single Modifier node
    if (pXMod == 1) {
      Eigen::MatrixXd tempV, VTheta, VThetaChol;
      Eigen::VectorXd XtVzInvR, ThetaHat, ThetaDraw;

      if (ctr->binomial) {
        const Eigen::MatrixXd Xdw = (ctr->Omega).asDiagonal() * X;  
        tempV           = X.transpose() *  Xdw;
        XtVzInvR        = Xdw.transpose() * ctr->R;
        
      } else {
        tempV           = X.transpose() * X;
        XtVzInvR        = X.transpose() * ctr->R;
      }

      tempV.noalias() -= ZtX.transpose() * VgZtX;
      tempV.diagonal().array() += 1.0 / treevar;
      VTheta            = tempV.inverse();
      XtVzInvR         -= VgZtX.transpose() * ZtR;
      VThetaChol        = VTheta.llt().matrixL();
      ThetaHat          = VTheta * XtVzInvR;
      out.draw          = ThetaHat;
      out.draw.noalias() += VThetaChol * as<Eigen::VectorXd>(rnorm(pXDlm, 0, sqrt(ctr->sigma2)));
      out.fitted        = X * out.draw;
      out.beta          = ThetaHat.dot(XtVzInvR);
      out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
      out.termT2        = (out.draw).dot(out.draw);
      out.nDlmTerm      = double(pXDlm); 
      out.nModTerm      = 1.0;
      out.totTerm       = out.nDlmTerm;

      return(out);
    }
  } // End if 1 modifier node

  // Multiple Modifier nodes
  Eigen::MatrixXd Xtemp, Ztemp;
  Eigen::VectorXd Rtemp, Otemp;
  Eigen::MatrixXd XXiblock(pXComb, pXComb); XXiblock.setZero();
  Eigen::VectorXd XtR(pXComb); XtR.setZero();
  Eigen::MatrixXd LInv(pXDlm, pXDlm); LInv.setZero();
  LInv.diagonal().array() += 1.0 / treevar;

  // Create block matrices corresponding to modifier nodes
  int j;
  int start = 0;
  for (Node* n : modTerm) {
    if (n->nodevals->updateXmat || ctr->binomial) { // update matrices for current node                                  
      Xtemp.resize(n->nodevals->idx.size(), pXDlm); Xtemp.setZero();
      Ztemp.resize(n->nodevals->idx.size(), ctr->pZ); Ztemp.setZero();
      Rtemp.resize(n->nodevals->idx.size()); Rtemp.setZero();
      Otemp.resize(n->nodevals->idx.size()); Otemp.setZero();

      j = 0;
      for (int i : n->nodevals->idx) {
        Xtemp.row(j) = X.row(i);
        Ztemp.row(j) = ctr->Zw.row(i);
        Rtemp(j) = ctr->R(i);
        if (ctr->binomial)
          Otemp(j) = ctr->Omega(i);
        ++j;
      } // end loop over node indices
      
      n->nodevals->XtX.resize(pXDlm, pXDlm); n->nodevals->XtX.setZero();
      n->nodevals->ZtXmat.resize(ctr->pZ, pXDlm); n->nodevals->ZtXmat.setZero();
      n->nodevals->VgZtXmat.resize(ctr->pZ, pXDlm); n->nodevals->VgZtXmat.setZero();
      if (ctr->binomial) {
        MatrixXd Xwtemp = Otemp.asDiagonal() * Xtemp;
        n->nodevals->XtX       = Xtemp.transpose() * Xwtemp;
        n->nodevals->ZtXmat    = Ztemp.transpose() * Xtemp;
        n->nodevals->VgZtXmat  = ctr->Vg * n->nodevals->ZtXmat;
        XtR.segment(start, pXDlm) = Xwtemp.transpose() * Rtemp;
      } else {
        n->nodevals->XtX       = Xtemp.transpose() * Xtemp;
        n->nodevals->ZtXmat    = Ztemp.transpose() * Xtemp;
        n->nodevals->VgZtXmat  = ctr->Vg * n->nodevals->ZtXmat;
        XtR.segment(start, pXDlm) = Xtemp.transpose() * Rtemp;
        n->nodevals->updateXmat = 0;
      }
      
      
    } else { // reuse precalculated matrices      
      for (int i : n->nodevals->idx) {
        XtR.segment(start, pXDlm).noalias() += X.row(i).transpose() * ctr->R(i);
      } // end loop over node indices
      
    } // end update xblock and ztx block
     
    if ((ctr->pZ < pXComb) && !ctr->binomial) {
      XXiblock.block(start, start, pXDlm, pXDlm)  = (n->nodevals->XtX + LInv).inverse();
    } else {
      XXiblock.block(start, start, pXDlm, pXDlm) = n->nodevals->XtX + LInv;
    }
    ZtX.block(0, start, ctr->pZ, pXDlm)   = n->nodevals->ZtXmat;
    VgZtX.block(0, start, ctr->pZ, pXDlm) = n->nodevals->VgZtXmat;
    
    start += pXDlm;
  } // end loop over modifier tree terminal nodes

  Eigen::MatrixXd VTheta(pXComb, pXComb); VTheta.setZero();
  if ((ctr->pZ < pXComb) && !ctr->binomial) {
    Eigen::MatrixXd ZtXXi = ZtX * XXiblock;
    VTheta.triangularView<Lower>() = XXiblock + ZtXXi.transpose() *
      (ctr->VgInv - ZtXXi * ZtX.transpose()).inverse() * ZtXXi;
  } else {
    VTheta.triangularView<Lower>() = (XXiblock - ZtX.transpose() * ctr->Vg * ZtX).inverse();
  }
  
  Eigen::VectorXd XtVzInvR    = XtR - VgZtX.transpose() * ZtR;
  Eigen::VectorXd ThetaHat    = VTheta.selfadjointView<Lower>() * XtVzInvR;
  Eigen::MatrixXd VThetaChol  = VTheta.selfadjointView<Lower>().llt().matrixL();

  // Calculate fitted values 
  // -> This is the same as out.Xd * out.draw in TDLMM code
  // -> Instead of calculating it outside of 'out' object, calculate and return "fitted" property 
  out.draw = ThetaHat;
  out.draw.noalias() += VThetaChol * as<Eigen::VectorXd>(rnorm(pXComb, 0.0, sqrt(ctr->sigma2)));
  out.fitted.resize(ctr->n);

  Eigen::VectorXd drawTemp(pXDlm);
  for (s = 0; s < modTerm.size(); ++s) {
    drawTemp = out.draw.segment(s * pXDlm, pXDlm); // pxDlm elements starting s*pXDlm
    for (int i : modTerm[s]->nodevals->idx){
      out.fitted(i) = X.row(i) * drawTemp;
    }
  }

  out.beta          = ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
  out.termT2        = (out.draw).dot(out.draw);
  out.nDlmTerm      = pXDlm * 1.0; 
  out.nModTerm      = pXMod * 1.0;
  out.totTerm       = out.nModTerm * out.nDlmTerm;
  
  return(out);
}