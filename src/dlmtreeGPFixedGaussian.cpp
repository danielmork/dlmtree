#include "RcppEigen.h"
#include "Node.h"
#include "NodeStruct.h"
#include "modDat.h"
#include "exposureDat.h"
#include "Fncs.h"
#include "modelCtr.h"
using namespace Rcpp;

void dlmtreeGPFixed_Gaussian_TreeMCMC(int t, std::vector<Node*> fixedNodes,
                                      dlmtreeCtr* ctr, dlmtreeLog *dgn);

treeMHR dlmtreeFixedMHR(std::vector<Node*> fixedNodes,
                        dlmtreeCtr* ctr, 
                        Eigen::VectorXd ZtR,
                        double treevar);


//' dlmtree model with fixed Gaussian process approach
//'
//' @param model A list of parameter and data contained for the model fitting
//' @returns A list of dlmtree model fit, mainly posterior mcmc samples
//' @export
// [[Rcpp::export]]
Rcpp::List dlmtreeGPFixedGaussian(const Rcpp::List model)
{
  int t;
  // ---- Set up general control variables ----
  dlmtreeCtr *ctr = new dlmtreeCtr;
  ctr->iter     = as<int>(model["nIter"]);
  ctr->burn     = as<int>(model["nBurn"]);
  ctr->thin     = as<int>(model["nThin"]);
  ctr->nRec     = floor(ctr->iter / ctr->thin);
  ctr->nTrees   = as<int>(model["nTrees"]);
  ctr->Y        = as<Eigen::VectorXd>(model["Y"]);
  ctr->n        = (ctr->Y).size();
  ctr->modZeta  = as<double>(model["zeta"]);
  ctr->modKappa = 100;
  ctr->binomial = 0;
  ctr->zinb     = 0;

  ctr->Z      = as<Eigen::MatrixXd>(model["Z"]);
  ctr->Zw     = ctr->Z;
  ctr->pZ     = (ctr->Z).cols();
  ctr->VgInv  = (ctr->Z).transpose() * (ctr->Z);
  ctr->VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg     = ctr->VgInv.inverse();
  ctr->VgChol = (ctr->Vg).llt().matrixL();

  ctr->X            = as<Eigen::MatrixXd>(model["X"]);
  ctr->pX           = ctr->X.cols();
  ctr->DistMat      = as<Eigen::MatrixXd>(model["DistMat"]);
  ctr->LambdaInv    = ctr->DistMat.array().exp().matrix().inverse();
  ctr->LambdaInvNew = ctr->LambdaInv;

  Eigen::MatrixXd LambdaInvChol = ctr->LambdaInv.llt().matrixL();
  ctr->logLambdaDet             = -2.0 * LambdaInvChol.diagonal().array().log().sum();
  ctr->logLambdaDetNew          = ctr->logLambdaDet;
  ctr->phi    = 1;
  ctr->phiNew = 1;
  double logphiLow  = log(-log(.95));
  double logphiHigh = log(-log(.05));
  double logphi     = log(ctr->phi);
  double logphiNew  = logphi;

  // ---- Create trees ----
  std::vector<Node*> fixedNodes;
  Rcpp::List fixedIdx = as<Rcpp::List>(model["fixedIdx"]);
  Eigen::MatrixXd Xtemp, Ztemp;
  for (t = 0; t < fixedIdx.size(); ++t) {
    fixedNodes.push_back(new Node(0, 1));
    fixedNodes[t]->nodevals       = new NodeVals(ctr->n);
    fixedNodes[t]->nodevals->idx  = as<std::vector<int> >(fixedIdx[t]);
    Xtemp.resize(fixedNodes[t]->nodevals->idx.size(), ctr->pX);
    Ztemp.resize(fixedNodes[t]->nodevals->idx.size(), ctr->pZ);

    for (std::size_t i = 0; i < fixedNodes[t]->nodevals->idx.size(); ++i) {
      Xtemp.row(i) = ctr->X.row(fixedNodes[t]->nodevals->idx[i]);
      Ztemp.row(i) = ctr->Z.row(fixedNodes[t]->nodevals->idx[i]);
    }

    fixedNodes[t]->nodevals->XtX        = Xtemp.transpose() * Xtemp;
    fixedNodes[t]->nodevals->ZtXmat     = Ztemp.transpose() * Xtemp;
    fixedNodes[t]->nodevals->VgZtXmat   = ctr->Vg * fixedNodes[t]->nodevals->ZtXmat;
    fixedNodes[t]->nodevals->updateXmat = 0;
  }


  // ---- Logs ----
  dlmtreeLog *dgn = new dlmtreeLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec);    (dgn->gamma).setZero();
  (dgn->sigma2).resize(ctr->nRec);            (dgn->sigma2).setZero();
  (dgn->nu).resize(ctr->nRec);                (dgn->nu).setZero();
  (dgn->phi).resize(ctr->nRec);               (dgn->phi).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec);  (dgn->tau).setZero();
  (dgn->fhat).resize(ctr->n);                 (dgn->fhat).setZero();

  // ---- DLM estimates ----
  // dgn->exDLM.resize(ctr->pX, ctr->n); dgn->exDLM.setZero();
  // dgn->ex2DLM.resize(ctr->pX, ctr->n); dgn->ex2DLM.setZero();
  // dgn->cumDLM.resize(ctr->n); dgn->cumDLM.setZero();
  // dgn->cum2DLM.resize(ctr->n); dgn->cum2DLM.setZero();

  // ---- Initial draws ----
  (ctr->fhat).resize(ctr->n);                 (ctr->fhat).setZero();
  ctr->R = ctr->Y;
  (ctr->gamma).resize(ctr->pZ);

  // Load initial params for faster convergence in binomial model
  if (ctr->binomial) {
    ctr->gamma  = as<VectorXd>(model["initParams"]);
    ctr->Omega  = rcpp_pgdraw(ctr->binomialSize, ctr->fhat + ctr->Z * ctr->gamma);
    ctr->Zw     = ctr->Omega.asDiagonal() * ctr->Z;
    ctr->VgInv  = ctr->Z.transpose() * ctr->Zw;
    ctr->VgInv.diagonal().array() += 1 / 100000.0;
    ctr->Vg     = ctr->VgInv.inverse();
    ctr->VgChol = ctr->Vg.llt().matrixL();
    // recalculate 'pseudo-Y' = kappa / omega, kappa = (y - n_b)/2
    ctr->Y      = ctr->kappa.array() / ctr->Omega.array();
  }
  ctr->totTerm    = 0;
  ctr->sumTermT2  = 0;
  ctr->nu         = 1.0; // Need to define for first update of sigma2
  ctr->sigma2     = 1.0;
  tdlmModelEst(ctr);
  double xiInv    = R::rgamma(1, 0.5);
  ctr->nu         = 1.0 / R::rgamma(0.5 * ctr->nTrees + 0.5, 1.0 / xiInv);
  (ctr->tau).resize(ctr->nTrees);           (ctr->tau).setOnes();
  for (t = 0; t < ctr->nTrees; t++) {
    xiInv         = R::rgamma(1, 0.5);
    (ctr->tau)(t) = 1.0 / R::rgamma(0.5, 1.0 / xiInv);
  }
  // ctr->exDLM.resize(ctr->pX, ctr->n);
  (ctr->Rmat).resize(ctr->n, ctr->nTrees);  (ctr->Rmat).setZero();
  
  // Create progress meter
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
    ctr->fhat.setZero();
    ctr->totTerm    = 0.0; 
    ctr->sumTermT2  = 0.0;
    // ctr->exDLM.setZero();
    ctr->phiMH      = 0; 
    ctr->phiMHNew   = 0;
    for (t = 0; t < ctr->nTrees; t++) {
      dlmtreeGPFixed_Gaussian_TreeMCMC(t, fixedNodes, ctr, dgn);
      ctr->fhat += (ctr->Rmat).col(t);

      if (t < ctr->nTrees - 1){
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
      }
    }

    // -- Update model --
    ctr->R    = ctr->Y - ctr->fhat;
    tdlmModelEst(ctr);
    xiInv     = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->nu)));
    ctr->nu   = 1.0 / R::rgamma(0.5 * ctr->totTerm + 0.5,
                               1.0 / (0.5 * ctr->sumTermT2 / (ctr->sigma2) + xiInv));

    if ((ctr->sigma2 != ctr->sigma2) || (ctr->nu != ctr->nu)) {
      // Rcout << "\n" << ctr->sigma2 << "\n" <<
      //   ctr->nu << "\n" << ctr->totTerm << "\n" << ctr->sumTermT2 << "\n" <<
      //   ctr->xiInvSigma2 << "\n" << xiInv << "\n" <<
      //   ctr->tau;
      stop("\nNaN values occured during model run, rerun model.\n");
    }
    
    
    // -- Update phi for GP --
    double phiMHRatio = (ctr->logLambdaDet - ctr->logLambdaDetNew) * 
      (ctr->totTerm / (2.0 * ctr->pX)) +
      (ctr->phiMHNew - ctr->phiMH) / (2.0 * ctr->nu * ctr->sigma2) +
      (R::dgamma(ctr->phiNew, 0.5, 2.0, 1) - 
       R::dgamma(ctr->phi, 0.5, 2.0, 1));
    if (log(R::runif(0, 1) < phiMHRatio)) {
      ctr->phi          = ctr->phiNew;
      logphi            = logphiNew;
      ctr->LambdaInv    = ctr->LambdaInvNew;
      ctr->logLambdaDet = ctr->logLambdaDetNew;
    }
    // propose new phi
    logphiNew = logphi + R::rnorm(0, 0.3);
    if (logphiNew < logphiLow){
      logphiNew = logphiLow + abs(logphiNew - logphiLow);
    }

    if (logphiNew > logphiHigh){
      logphiNew = logphiLow + abs(logphiNew - logphiHigh);
    }

    ctr->phiNew           = exp(logphiNew);
    ctr->LambdaInvNew     = (ctr->phiNew * ctr->DistMat).array().exp().matrix().inverse();
    LambdaInvChol         = ctr->LambdaInvNew.llt().matrixL();
    ctr->logLambdaDetNew  = -2.0 * LambdaInvChol.diagonal().array().log().sum();
    
    // -- Record --
    if (ctr->record > 0) {
      (dgn->gamma).col(ctr->record - 1) = ctr->gamma;
      (dgn->sigma2)(ctr->record - 1)    = ctr->sigma2;
      (dgn->nu)(ctr->record - 1)        = ctr->nu;
      (dgn->tau).col(ctr->record - 1)   = ctr->tau;
      (dgn->phi)(ctr->record - 1)       = ctr->phi;
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
  
  Eigen::MatrixXd TreeStructs((dgn->DLMexp).size(), 3 + ctr->pX);
  for (s = 0; s < (dgn->DLMexp).size(); ++s){
    TreeStructs.row(s) = dgn->DLMexp[s];
  }

  Eigen::VectorXd sigma2  = dgn->sigma2;
  Eigen::VectorXd nu      = dgn->nu;
  Eigen::MatrixXd tau     = (dgn->tau).transpose();
  Eigen::VectorXd fhat    = (dgn->fhat).array() / ctr->nRec;
  Eigen::MatrixXd gamma   = (dgn->gamma).transpose();
  Eigen::VectorXd phi     = dgn->phi;

  delete prog;
  delete ctr;
  delete dgn;
  for (s = 0; s < fixedNodes.size(); ++s) {
    delete fixedNodes[s];
  }

  return(Rcpp::List::create(// Named("DLM") = wrap(exDLM),
                            // Named("DLMse") = wrap(ex2DLM),
                            // Named("DLfun") = wrap(cumDLM),
                            // Named("DLfunse") = wrap(cum2DLM),
                            Named("TreeStructs")  = wrap(TreeStructs),
                            Named("fhat")         = wrap(fhat),
                            Named("sigma2")       = wrap(sigma2),
                            Named("nu")           = wrap(nu),
                            Named("tau")          = wrap(tau),
                            Named("gamma")        = wrap(gamma),
                            Named("phi")          = wrap(phi)));

} // end dlmtreeGPGaussian



void dlmtreeGPFixed_Gaussian_TreeMCMC(int t, std::vector<Node*> fixedNodes,
                                      dlmtreeCtr* ctr, dlmtreeLog *dgn)
{
  double treevar = (ctr->nu) * (ctr->tau)(t);
  std::size_t s;
  Eigen::VectorXd ZtR = (ctr->Z).transpose() * (ctr->R);
  treeMHR mhr0        = dlmtreeFixedMHR(fixedNodes, ctr, ZtR, treevar);
  
  // -- Update variance and residuals --
  double xiInv      = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->tau)(t)));
  (ctr->tau)(t)     = 1.0 / R::rgamma(0.5 * mhr0.draw.size() + 0.5,
                                  1.0 / ((0.5 * mhr0.termT2 / (ctr->sigma2 * ctr->nu)) + xiInv));
  ctr->Rmat.col(t)  = mhr0.fitted;
  ctr->sumTermT2 += mhr0.termT2 / (ctr->tau(t));
  ctr->totTerm += static_cast<double>(mhr0.draw.size());
  
  // -- calculate full conditionals for phi update --
  Eigen::VectorXd draw(ctr->pX);
  Eigen::VectorXd rec(3 + ctr->pX);
  for (s = 0; s < fixedNodes.size(); ++s) {
    draw = mhr0.draw.segment(s * ctr->pX, ctr->pX);
    ctr->phiMH -= draw.dot(ctr->LambdaInv * draw) / ctr->tau(t);
    ctr->phiMHNew -= draw.dot(ctr->LambdaInvNew * draw) / ctr->tau(t);
    
    // -- Update DLM partial estimate --
    if (ctr->record > 0) {
      rec << ctr->record, t, s, draw;
      dgn->DLMexp.push_back(rec);
      // for (int i : fixedNodes[s]->nodevals->idx)
      //   ctr->exDLM.col(i) += draw;
    } // end record
  } // end loop over fixedNodes
} // end dlmtreeGP_Gaussian_TreeMCMC function



// function to calculate part of MH ratio
treeMHR dlmtreeFixedMHR(std::vector<Node*> fixedNodes,
                   dlmtreeCtr* ctr, 
                   Eigen::VectorXd ZtR,
                   double treevar)
{
  std::size_t s;
  treeMHR out;
  int pX = ctr->pX * fixedNodes.size();
  Eigen::MatrixXd Linv = ctr->LambdaInv / treevar;

  // Multiple Modifier nodes
  Eigen::MatrixXd Xtemp, Ztemp;
  Eigen::MatrixXd XXiblock(pX, pX); XXiblock.setZero();
  Eigen::MatrixXd ZtX(ctr->pZ, pX); ZtX.setZero();
  Eigen::MatrixXd VgZtX(ctr->pZ, pX); ZtX.setZero();
  Eigen::VectorXd XtR(pX); XtR.setZero();

  // Create block matrices corresponding to modifier nodes
  int start = 0;
  for (Node* n : fixedNodes) {      
    XXiblock.block(start, start, ctr->pX, ctr->pX)  = (n->nodevals->XtX + Linv).inverse();
    ZtX.block(0, start, ctr->pZ, ctr->pX)           = n->nodevals->ZtXmat;
    VgZtX.block(0, start, ctr->pZ, ctr->pX)         = n->nodevals->VgZtXmat;
    
    for (int i : n->nodevals->idx) {
      XtR.segment(start, ctr->pX).noalias() += ctr->X.row(i).transpose() * ctr->R(i);
    }
    
    start += ctr->pX;
  } // end loop over fixedNodes
  
  Eigen::MatrixXd ZtXXi       = ZtX * XXiblock;
  Eigen::MatrixXd VTheta      = XXiblock;
  VTheta.noalias() += ZtXXi.transpose() * (ctr->VgInv - ZtXXi * ZtX.transpose()).inverse() * ZtXXi;
  Eigen::VectorXd XtVzInvR  = XtR;
  XtVzInvR.noalias() -= VgZtX.transpose() * ZtR;
  Eigen::VectorXd ThetaHat    = VTheta * XtVzInvR;
  Eigen::MatrixXd VThetaChol  = VTheta.llt().matrixL();

  // Calculate fitted values
  out.draw    = ThetaHat;
  out.draw.noalias() += VThetaChol * as<Eigen::VectorXd>(rnorm(pX, 0.0, sqrt(ctr->sigma2)));
  out.fitted.resize(ctr->n);
  out.termT2  = 0.0;
  Eigen::VectorXd drawTemp(ctr->pX);
  for (s = 0; s < fixedNodes.size(); ++s) {
    drawTemp = out.draw.segment(s * ctr->pX, ctr->pX);
    out.termT2 += drawTemp.dot(ctr->LambdaInv * drawTemp);
    for (int i : fixedNodes[s]->nodevals->idx)
      out.fitted(i) = ctr->X.row(i) * drawTemp;
  }
  
  out.beta          = ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
  
  return(out);
} // end dlmtreeMHR function
