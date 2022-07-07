#include "RcppEigen.h"
#include "Node.h"
#include "NodeStruct.h"
#include "modDat.h"
#include "exposureDat.h"
#include "Fncs.h"
#include "modelCtr.h"
using namespace Rcpp;


void dlmtreeGP_Gaussian_TreeMCMC(int t, Node* modTree,
                                   dlmtreeCtr* ctr, dlmtreeLog *dgn,
                                   modDat* Mod);

treeMHR dlmtreeGP_MHR(std::vector<Node*> modTerm,
                   dlmtreeCtr* ctr, 
                   Eigen::VectorXd ZtR,
                   double treevar);


// [[Rcpp::export]]
Rcpp::List dlmtreeGPGaussian(const Rcpp::List model)
{
  int t;
  // ---- Set up general control variables ----
  dlmtreeCtr *ctr = new dlmtreeCtr;
  ctr->iter = as<int>(model["nIter"]);
  ctr->burn = as<int>(model["nBurn"]);
  ctr->thin = as<int>(model["nThin"]);
  ctr->nRec = floor(ctr->iter / ctr->thin);
  ctr->nTrees = as<int>(model["nTrees"]);
  ctr->Y = as<Eigen::VectorXd>(model["Y"]);
  ctr->n = (ctr->Y).size();
  ctr->modZeta = as<double>(model["zeta"]);
  ctr->modKappa = 100;

  ctr->Z = as<Eigen::MatrixXd>(model["Z"]);
  ctr->Zw = ctr->Z;
  ctr->pZ = (ctr->Z).cols();
  ctr->VgInv = (ctr->Z).transpose() * (ctr->Z);
  ctr->VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg = ctr->VgInv.inverse();
  ctr->VgChol = (ctr->Vg).llt().matrixL();
  
  ctr->binomial = 0;
  ctr->verbose = bool (model["verbose"]);
  ctr->diagnostics = bool (model["diagnostics"]);
  ctr->stepProbMod = as<std::vector<double> >(model["stepProbMod"]);
  ctr->treePriorMod = as<std::vector<double> >(model["treePriorTDLM"]);
  ctr->shrinkage = as<int>(model["shrinkage"]);

  // ---- Setup modifier data ----
  modDat *Mod = new modDat(as<std::vector<int> >(model["modIsNum"]),
                           as<Rcpp::List>(model["modSplitIdx"]),
                           as<std::vector<int> >(model["fullIdx"]));

  NodeStruct *modNS;
  modNS = new ModStruct(Mod);
  ctr->pM = Mod->nMods;


  ctr->X = as<Eigen::MatrixXd>(model["X"]);
  ctr->pX = ctr->X.cols();
  ctr->XtXall = ctr->X.transpose() * ctr->X;
  ctr->ZtXall = ctr->Z.transpose() * ctr->X;
  ctr->VgZtXall = ctr->Vg * ctr->ZtXall;
  ctr->VThetaInvall = ctr->XtXall - ctr->ZtXall.transpose() * ctr->VgZtXall;
  
  // ---- Setup covariance matrix ----
  ctr->covarType = as<int>(model["covarianceType"]);
  ctr->DistMat = as<Eigen::MatrixXd>(model["DistMat"]);
  Eigen::MatrixXd LambdaInvChol(ctr->pX, ctr->pX);
  ctr->LambdaInv.resize(ctr->pX, ctr->pX); ctr->LambdaInv.setZero();
  ctr->LambdaInvNew.resize(ctr->pX, ctr->pX); ctr->LambdaInvNew.setZero();
  ctr->logLambdaDet = 0; ctr->logLambdaDetNew = 0;
  ctr->phi = 1;
  ctr->phiNew = 1;
  double logphiLow = log(-log(.95));
  double logphiHigh = log(-log(.05));
  double logphi = log(ctr->phi);
  double logphiNew = logphi;
  if (ctr->covarType == 1) {
    ctr->LambdaInv = ctr->DistMat.array().exp().matrix().inverse();
    ctr->LambdaInvNew = ctr->LambdaInv;
    LambdaInvChol = ctr->LambdaInv.llt().matrixL();
    ctr->logLambdaDet = -2.0 * LambdaInvChol.diagonal().array().log().sum();
    ctr->logLambdaDetNew = ctr->logLambdaDet;
  } else {
    ctr->LambdaInv.diagonal().array() += 1;
    ctr->logLambdaDet = 0;
  }

  // ---- Create trees ----
  std::vector<Node*> modTrees;
  for (t = 0; t < ctr->nTrees; ++t) {
    modTrees.push_back(new Node(0, 1));
    modTrees[t]->nodestruct = modNS->clone();
    Mod->updateNodeVals(modTrees[t]);
    updateGPMats(modTrees[t], ctr);
  } 
  delete modNS;


  // ---- Logs ----
  dlmtreeLog *dgn = new dlmtreeLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec); (dgn->gamma).setZero();
  (dgn->sigma2).resize(ctr->nRec); (dgn->sigma2).setZero();
  (dgn->nu).resize(ctr->nRec); (dgn->nu).setZero();
  (dgn->phi).resize(ctr->nRec); (dgn->phi).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec); (dgn->tau).setZero();
  (dgn->fhat).resize(ctr->n); (dgn->fhat).setZero();
  (dgn->modProb).resize(ctr->pM, ctr->nRec); (dgn->modProb).setZero();
  (dgn->modCount).resize(ctr->pM, ctr->nRec); (dgn->modCount).setZero();
  (dgn->modInf).resize(ctr->pM, ctr->nRec); (dgn->modInf).setZero();
  (dgn->modKappa).resize(ctr->nRec); (dgn->modKappa).setZero();

  (dgn->termNodesMod).resize(ctr->nTrees, ctr->nRec);
    (dgn->termNodesMod).setZero();

  // ---- DLM estimates ----
  // dgn->exDLM.resize(ctr->pX, ctr->n); dgn->exDLM.setZero();
  // dgn->ex2DLM.resize(ctr->pX, ctr->n); dgn->ex2DLM.setZero();
  // dgn->cumDLM.resize(ctr->n); dgn->cumDLM.setZero();
  // dgn->cum2DLM.resize(ctr->n); dgn->cum2DLM.setZero();

  // ---- Initial draws ----
  (ctr->fhat).resize(ctr->n); (ctr->fhat).setZero();
  ctr->R = ctr->Y;
  (ctr->gamma).resize(ctr->pZ);
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  ctr->nu = 1.0; // Need to define for first update of sigma2
  ctr->sigma2 = 1.0;
  tdlmModelEst(ctr);
  double xiInv = R::rgamma(1, 0.5);
  ctr->nu = 1.0 / R::rgamma(0.5 * ctr->nTrees + 0.5, 1.0 / xiInv);
  (ctr->tau).resize(ctr->nTrees); (ctr->tau).setOnes();
  if (ctr->shrinkage) {
    for (t = 0; t < ctr->nTrees; t++) {
      xiInv = R::rgamma(1, 0.5);
      (ctr->tau)(t) = 1.0 / R::rgamma(0.5, 1.0 / xiInv);
    }
  }
  ctr->nTermMod.resize(ctr->nTrees);
  ctr->nTermMod.array() = 1;
  (ctr->Rmat).resize(ctr->n, ctr->nTrees); (ctr->Rmat).setZero();
  ctr->modCount.resize(ctr->pM); ctr->modCount.setZero();
  ctr->modInf.resize(ctr->pM); ctr->modInf.setZero();
  // ctr->exDLM.resize(ctr->pX, ctr->n); ctr->exDLM.setZero();
  
  // create progress meter
  progressMeter* prog = new progressMeter(ctr);

  std::size_t s;
  // ---- MCMC ----
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }

    // -- Update trees --
    ctr->R += (ctr->Rmat).col(0);
    ctr->fhat.setZero();
    ctr->totTerm = 0.0; ctr->sumTermT2 = 0.0;
    // ctr->exDLM.setZero();
    ctr->modCount.setZero();
    ctr->modInf.setZero();
    ctr->phiMH = 0; ctr->phiMHNew = 0;
    for (t = 0; t < ctr->nTrees; t++) {
      dlmtreeGP_Gaussian_TreeMCMC(t, modTrees[t], ctr, dgn, Mod);
      ctr->fhat += (ctr->Rmat).col(t);

      if (t < ctr->nTrees - 1)
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
    }

    // -- Update model --
    ctr->R = ctr->Y - ctr->fhat;
    tdlmModelEst(ctr);
    xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->nu)));
    ctr->nu = 1.0 / R::rgamma(0.5 * ctr->totTerm + 0.5,
                               1.0 / (0.5 * ctr->sumTermT2 / (ctr->sigma2) +
                                      xiInv));

    // -- Update modifier selection --
    if ((ctr->b > 1000) || (ctr->b > (0.5 * ctr->burn))) {
      double beta = R::rbeta(ctr->modZeta, 1.0);
      double modKappaNew = beta * ctr->pM / (1 - beta);
      double mhrDir =
        logDirichletDensity(Mod->modProb,
                            (ctr->modCount.array() +
                             modKappaNew / ctr->pM).matrix()) -
        logDirichletDensity(Mod->modProb,
                            (ctr->modCount.array() +
                             ctr->modKappa / ctr->pM).matrix());
      if (log(R::runif(0, 1)) < mhrDir) {
        ctr->modKappa = modKappaNew;
      }

      Mod->modProb = rDirichlet((ctr->modCount.array() +
                                 ctr->modKappa / ctr->pM).matrix());
    } // end modifier selection
    
    
    // -- Update phi for GP --
    if (ctr->covarType == 1) {
      double phiMHRatio = (ctr->logLambdaDet - ctr->logLambdaDetNew) * 
        (ctr->totTerm / (2.0 * ctr->pX)) +
        (ctr->phiMHNew - ctr->phiMH) / (2.0 * ctr->nu * ctr->sigma2) +
        (R::dgamma(ctr->phiNew, 0.5, 2.0, 1) - 
        R::dgamma(ctr->phi, 0.5, 2.0, 1));
      if (log(R::runif(0, 1) < phiMHRatio)) {
        ctr->phi = ctr->phiNew;
        logphi = logphiNew;
        ctr->LambdaInv = ctr->LambdaInvNew;
        ctr->logLambdaDet = ctr->logLambdaDetNew;
      }
      // propose new phi
      logphiNew = logphi + R::rnorm(0, 0.3);
      if (logphiNew < logphiLow)
        logphiNew = logphiLow + abs(logphiNew - logphiLow);
      if (logphiNew > logphiHigh)
        logphiNew = logphiLow + abs(logphiNew - logphiHigh);
      ctr->phiNew = exp(logphiNew);
      ctr->LambdaInvNew = 
        (ctr->phiNew * ctr->DistMat).array().exp().matrix().inverse();
      LambdaInvChol = ctr->LambdaInvNew.llt().matrixL();
      ctr->logLambdaDetNew = -2.0*LambdaInvChol.diagonal().array().log().sum();
    }
    



    // -- Record --
    if (ctr->record > 0) {
      (dgn->gamma).col(ctr->record - 1) = ctr->gamma;
      (dgn->sigma2)(ctr->record - 1) = ctr->sigma2;
      (dgn->nu)(ctr->record - 1) = ctr->nu;
      (dgn->tau).col(ctr->record - 1) = ctr->tau;
      (dgn->termNodesMod).col(ctr->record - 1) = ctr->nTermMod;
      (dgn->modProb).col(ctr->record - 1) = Mod->modProb;
      (dgn->modCount).col(ctr->record - 1) = ctr->modCount;
      (dgn->modInf).col(ctr->record - 1) = ctr->modInf / ctr->modInf.maxCoeff();
      (dgn->modKappa)(ctr->record - 1) = ctr->modKappa;
      (dgn->phi)(ctr->record - 1) = ctr->phi;
      dgn->fhat += ctr->fhat;
      // dlmtreeRecDLM(ctr, dgn);
    } // end record



    // -- Progress --
    prog->printMark();
  } // end MCMC


  // -- Prepare outout --
  // Eigen::MatrixXd exDLM = dgn->exDLM.transpose();
  // Eigen::MatrixXd ex2DLM = dgn->ex2DLM.transpose();
  // Eigen::VectorXd cumDLM = dgn->cumDLM;
  // Eigen::VectorXd cum2DLM = dgn->cum2DLM;
  
  Rcpp::StringVector termRule(dgn->termRule.size());
  Eigen::MatrixXd TreeStructs((dgn->DLMexp).size(), 3 + ctr->pX);
  for (s = 0; s < (dgn->DLMexp).size(); ++s)
    TreeStructs.row(s) = dgn->DLMexp[s];
  termRule = dgn->termRule;

  Eigen::VectorXd sigma2 = dgn->sigma2;
  Eigen::VectorXd nu = dgn->nu;
  Eigen::MatrixXd tau = (dgn->tau).transpose();
  Eigen::VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  Eigen::MatrixXd gamma = (dgn->gamma).transpose();
  Eigen::VectorXd phi = dgn->phi;

  Eigen::MatrixXd termNodesMod = (dgn->termNodesMod).transpose();
  Eigen::VectorXd kappa = dgn->modKappa;
  Eigen::MatrixXd modProb = (dgn->modProb).transpose();
  Eigen::MatrixXd modCount = (dgn->modCount).transpose();
  Eigen::MatrixXd modInf = (dgn->modInf).transpose();

  Eigen::MatrixXd modAccept((dgn->treeModAccept).size(), 5);

  delete prog;
  delete ctr;
  delete dgn;
  delete Mod;
  for (s = 0; s < modTrees.size(); ++s) {
    delete modTrees[s];
  }

  return(Rcpp::List::create(// Named("DLM") = wrap(exDLM),
                            // Named("DLMse") = wrap(ex2DLM),
                            // Named("DLfun") = wrap(cumDLM),
                            // Named("DLfunse") = wrap(cum2DLM),
                            Named("TreeStructs") = wrap(TreeStructs),
                            Named("termRules") = wrap(termRule),
                            Named("fhat") = wrap(fhat),
                            Named("sigma2") = wrap(sigma2),
                            Named("nu") = wrap(nu),
                            Named("tau") = wrap(tau),
                            Named("gamma") = wrap(gamma),
                            Named("phi") = wrap(phi),
                            Named("termNodesMod") = wrap(termNodesMod),
                            Named("kappa") = wrap(kappa),
                            Named("modProb") = wrap(modProb),
                            Named("modCount") = wrap(modCount),
                            Named("modInf") = wrap(modInf),
                            Named("treeModAccept") = wrap(modAccept)));

} // end dlmtreeGPGaussian



void dlmtreeGP_Gaussian_TreeMCMC(int t, Node* modTree,
                                   dlmtreeCtr* ctr, dlmtreeLog *dgn,
                                   modDat* Mod)
{
  int step;
  int success = 0;
  double stepMhr = 0;
  double ratio = 0;
  double treevar = (ctr->nu) * (ctr->tau)(t);
  std::size_t s;
  std::vector<Node*> modTerm, newModTerm;
  Eigen::VectorXd ZtR = (ctr->Z).transpose() * (ctr->R);
  double RtR = (ctr->R).dot(ctr->R);
  double RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);
  treeMHR mhr0, mhr;

  // -- modifier tree proposal --
  modTerm = modTree->listTerminal(); 
  mhr0 = dlmtreeGP_MHR(modTerm, ctr, ZtR, treevar);
  switch (modTerm.size()) {
    case 1: step = 0; break;
    case 2: step = sampleInt(ctr->stepProbMod, 1 - ctr->stepProbMod[3]); break;
    default: step = sampleInt(ctr->stepProbMod, 1);
  }  
  stepMhr = modProposeTree(modTree, Mod, ctr, step);
  success = modTree->isProposed();

  if (success && (stepMhr == stepMhr)) {
    newModTerm = modTree->listTerminal(1);
    mhr = dlmtreeGP_MHR(newModTerm, ctr, ZtR, treevar);
    ratio = stepMhr +
      mhr.logVThetaChol - mhr0.logVThetaChol -
      (0.5 * (ctr->n + 1.0) *
        (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
         log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2)));
    if (step == 0)
      ratio -= 0.5 * (log(treevar) * ctr->pX + ctr->logLambdaDet);
    if (step == 1)
      ratio += 0.5 * (log(treevar) * ctr->pX + ctr->logLambdaDet);
    
    if (log(R::runif(0, 1)) < ratio) {
      mhr0 = mhr;
      success = 2;
      modTree->accept();
      modTerm = modTree->listTerminal();
    }
  } // end modTree proposal
  modTree->reject();


  // -- Update variance and residuals --
  if (ctr->shrinkage) {
    double xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->tau)(t)));
    (ctr->tau)(t) = 1.0 / R::rgamma(0.5 * mhr0.draw.size() + 0.5,
                                    1.0 / ((0.5 * mhr0.termT2 /
                                            (ctr->sigma2 * ctr->nu)) + xiInv));
  }
  ctr->Rmat.col(t) = mhr0.fitted;
  ctr->sumTermT2 += mhr0.termT2 / (ctr->tau(t));
  ctr->totTerm += static_cast<double>(mhr0.draw.size());
  ctr->nTermMod(t) = static_cast<double>(modTerm.size());
  
  // -- calculate full conditionals for phi update --
  std::string rule;
  Eigen::VectorXd rec(3 + ctr->pX);
  Eigen::VectorXd draw(ctr->pX);
  for (s = 0; s < modTerm.size(); ++s) {
    draw = mhr0.draw.segment(s * ctr->pX, ctr->pX);
    if (ctr->covarType == 1) {
      ctr->phiMH -= draw.dot(ctr->LambdaInv * draw) / ctr->tau(t);
      ctr->phiMHNew -= draw.dot(ctr->LambdaInvNew * draw) / ctr->tau(t);
    }
    
    // -- Update DLM partial estimate --
    if (ctr->record > 0) {
      rule = modRuleStr(modTerm[s], Mod);
      rec << ctr->record, t, s, draw;
      dgn->termRule.push_back(rule);
      dgn->DLMexp.push_back(rec);
      // for (int i : modTerm[s]->nodevals->idx)
      //   ctr->exDLM.col(i) += draw;
    } // end record
  } // end loop over modTerm
  
  // -- Count modifiers used in tree --
  Eigen::VectorXd modCount = countMods(modTree, Mod);
  ctr->modCount += modCount;

  // -- Record --
  if (ctr->record > 0) {    
    for (int i = 0; i < modCount.size(); ++i) {
      if (modCount(i) > 0)
        ctr->modInf(i) += ctr->tau(t);
    }
  } // end record
} // end dlmtreeGP_Gaussian_TreeMCMC function



// function to calculate part of MH ratio
treeMHR dlmtreeGP_MHR(std::vector<Node*> modTerm,
                   dlmtreeCtr* ctr, 
                   Eigen::VectorXd ZtR,
                   double treevar)
{
  std::size_t s;
  treeMHR out;
  int pX = ctr->pX * modTerm.size();
  Eigen::MatrixXd Linv = ctr->LambdaInv / treevar;

  // Multiple Modifier nodes
  Eigen::MatrixXd Xtemp, Ztemp;
  Eigen::VectorXd Rtemp;
  Eigen::MatrixXd XXiblock(pX, pX); XXiblock.setZero();
  Eigen::MatrixXd ZtX(ctr->pZ, pX); ZtX.setZero();
  Eigen::MatrixXd VgZtX(ctr->pZ, pX); ZtX.setZero();
  Eigen::VectorXd XtR(pX); XtR.setZero();

  // Create block matrices corresponding to modifier nodes
  int start = 0;
  for (Node* n : modTerm) {
    if (n->nodevals->updateXmat)
      updateGPMats(n, ctr);
      
    XXiblock.block(start, start, ctr->pX, ctr->pX) = 
      (n->nodevals->XtX + Linv).inverse();
    ZtX.block(0, start, ctr->pZ, ctr->pX) = n->nodevals->ZtXmat;
    VgZtX.block(0, start, ctr->pZ, ctr->pX) = n->nodevals->VgZtXmat;
    
    for (int i : n->nodevals->idx) {
      XtR.segment(start, ctr->pX).noalias() += 
        ctr->X.row(i).transpose() * ctr->R(i);
    }
    
    start += ctr->pX;
  } // end loop over modTerm
  
  Eigen::MatrixXd ZtXXi = ZtX * XXiblock;
  Eigen::MatrixXd VTheta = XXiblock;
  VTheta.noalias() += ZtXXi.transpose() *
    (ctr->VgInv - ZtXXi * ZtX.transpose()).inverse() * ZtXXi;
  Eigen::VectorXd XtVzInvR = XtR;
  XtVzInvR.noalias() -= VgZtX.transpose() * ZtR;
  Eigen::VectorXd ThetaHat = VTheta * XtVzInvR;
  Eigen::MatrixXd VThetaChol = VTheta.llt().matrixL();

  // Calculate fitted values
  out.draw = ThetaHat;
  out.draw.noalias() +=
    VThetaChol * as<Eigen::VectorXd>(rnorm(pX, 0.0, sqrt(ctr->sigma2)));
  out.fitted.resize(ctr->n);
  out.termT2 = 0.0;
  Eigen::VectorXd drawTemp(ctr->pX);
  for (s = 0; s < modTerm.size(); ++s) {
    drawTemp = out.draw.segment(s * ctr->pX, ctr->pX);
    out.termT2 += drawTemp.dot(ctr->LambdaInv * drawTemp);
    for (int i : modTerm[s]->nodevals->idx)
      out.fitted(i) = ctr->X.row(i) * drawTemp;
  }
  
  out.beta = ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
  return(out);
} // end dlmtreeGP_MHR function
