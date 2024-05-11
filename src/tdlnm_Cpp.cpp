/** @brief Method for Treed DLNM and Treed DLM
 * @version 1.0
 * @date 2021-02-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <RcppEigen.h>
#include "modelCtr.h"
#include "exposureDat.h"
#include "Node.h"
#include "NodeStruct.h"
#include "Fncs.h"
#include <random>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Lower;

// #ifdef _OPENMP
//   #include <omp.h>
//   // [[Rcpp::plugins(openmp)]]
// #else
//   #define omp_get_num_threads()  1
//   #define omp_get_thread_num()   0
//   #define omp_get_max_threads()  1
//   #define omp_get_thread_limit() 1
//   #define omp_get_num_procs()    1
//   #define omp_set_nested(a)   // empty statement to remove the call
//   #define omp_get_wtime()        0
// #endif

/**
 * @brief calculate half of metropolis-hastings ratio for given tree
 * 
 * @param nodes vector of Node pointers
 * @param ctr control data for model
 * @param ZtR Z^T * R
 * @param var nu*tau
 * @param tree pointer to top of tree
 * @param newTree if true, recalculate node-specific values
 * @return treeMHR 
 */
treeMHR dlnmMHR(std::vector<Node*> nodes, tdlmCtr *ctr,
                VectorXd ZtR, double var, Node* tree, bool newTree)
{
  treeMHR out;
  int pX = int(nodes.size());
  
  if ((pX == 1) && (!ctr->binomial) && (!ctr->zinb)) { // single terminal node, cont. response
    double VTheta = var / (var * ctr->VTheta1Inv + 1.0);
    double XtVzInvR = (ctr->X1).dot(ctr->R) - (ctr->VgZtX1).dot(ZtR);
    double ThetaHat = VTheta * XtVzInvR;
    double VThetaChol = sqrt(VTheta);

    out.Xd.resize(ctr->n, 1);
    out.Xd.col(0) = ctr->X1;
    out.draw.resize(1);
    out.draw(0) = VThetaChol * R::rnorm(0, sqrt(ctr->sigma2)) + ThetaHat;
    out.beta = ThetaHat * XtVzInvR;
    out.logVThetaChol = log(VThetaChol);

  } else { // 2+ terminal nodes or binomial response

    MatrixXd ZtX(ctr->pZ, pX);      ZtX.setZero();
    MatrixXd VgZtX(ctr->pZ, pX);    VgZtX.setZero();

    // * Create design Xd, Z^tX, and VgZ^tX matrices
    out.Xd.resize(ctr->n, pX);
    for (std::size_t s = 0; s < nodes.size(); ++s) {
      out.Xd.col(s) = (nodes[s]->nodevals)->X;
      
      if (ctr->binomial || ctr->zinb) {
        ZtX.col(s) = ctr->Zw.transpose() * (nodes[s]->nodevals)->X;
        VgZtX.col(s) = ctr->Vg * ZtX.col(s);
        
      } else {
        ZtX.col(s) = (nodes[s]->nodevals)->ZtX;
        VgZtX.col(s) = (nodes[s]->nodevals)->VgZtX;
      }
    }

    // * Calculate covariance matrix V_theta
    MatrixXd tempV(pX, pX);
    VectorXd XtVzInvR(pX);
    
    if (ctr->binomial) {
      const MatrixXd Xdw = (ctr->Omega).asDiagonal() * out.Xd;
      tempV = Xdw.transpose() * out.Xd;
      tempV.noalias() -= ZtX.transpose() * VgZtX;
      XtVzInvR = Xdw.transpose() * ctr->R;
      
    } else if (ctr->zinb){ // ZINB (subsetting at-risk observations)
      Eigen::MatrixXd outXdstar = selectIndM(out.Xd, ctr->NBidx); // subsetting
      const Eigen::MatrixXd Xdw = (selectInd(ctr->omega2, ctr->NBidx)).asDiagonal() * outXdstar;
      tempV = Xdw.transpose() * outXdstar;                                 
      tempV.noalias() -= ZtX.transpose() * VgZtX;
      XtVzInvR = Xdw.transpose() * selectInd(ctr->R, ctr->NBidx);

    } else {
      if (newTree) {
        tempV = out.Xd.transpose() * out.Xd;
        tempV.noalias() -= ZtX.transpose() * VgZtX;
        out.tempV = tempV;
      } else {
        tempV = tree->nodevals->tempV;
      }
      XtVzInvR = out.Xd.transpose() * ctr->R;
    }

    XtVzInvR.noalias() -= VgZtX.transpose() * ZtR;
    tempV.diagonal().array() += 1.0 / var;
    const MatrixXd VTheta = tempV.inverse();
    const MatrixXd VThetaChol = VTheta.llt().matrixL();
    const VectorXd ThetaHat = VTheta * XtVzInvR;
    
    // const VectorXd XtVzInvR =
    //   out.Xd.transpose() * ctr->R - VgZtX.transpose() * ZtR;

    out.draw = ThetaHat;
    out.draw.noalias() += VThetaChol * as<VectorXd>(rnorm(pX, 0, sqrt(ctr->sigma2)));
    out.beta = ThetaHat.dot(XtVzInvR);
    out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
  } // end 2+ terminal nodes

  out.termT2 = (out.draw).dot(out.draw);
  out.nTerm = double(pX);

  return(out);
} // end drawMHR

/**
 * @brief tree proposal and tree parameter sampling
 * 
 * @param t tree number
 * @param tree pointer to tree
 * @param ctr pointer to model control
 * @param dgn pointer to model log
 * @param Exp pointer to exposure data
 */
void tdlnmTreeMCMC(int t, Node *tree, tdlmCtr *ctr, tdlmLog *dgn, 
                   exposureDat *Exp)
{
  int step = 0;
  int success = 0;
  double stepMhr = 0.0;
  double ratio = 0.0;
  double treevar = ctr->nu * ctr->tau(t);
  std::size_t s;
  std::vector<Node*> dlnmTerm, newDlnmTerm;
  treeMHR mhr0, mhr;

  // List current tree terminal nodes
  dlnmTerm = tree->listTerminal();
  VectorXd ZtR = (ctr->Zw).transpose() * (ctr->R);
  mhr0 = dlnmMHR(dlnmTerm, ctr, ZtR, treevar, tree, 0);

  if (dlnmTerm.size() > 1) {
    step = sampleInt(ctr->stepProb, 1);
  } else { // if single terminal node, grow is only option
    step = 0;
  }

  // propose update
  stepMhr = tdlmProposeTree(tree, Exp, ctr, step);
  success = tree->isProposed();

  if (success) {
    // calculate new tree part of MHR and draw node effects
    newDlnmTerm = tree->listTerminal(1);
    mhr = dlnmMHR(newDlnmTerm, ctr, ZtR, treevar, tree, 1);

    // combine mhr parts into log-MH ratio
    if (ctr->binomial || ctr->zinb) {
      ratio = stepMhr + (mhr.logVThetaChol - mhr0.logVThetaChol) +
        0.5 * (mhr.beta - mhr0.beta) -
        (log(treevar) * 0.5 * (mhr.nTerm - mhr0.nTerm));
        
    } else {
      double RtR, RtZVgZtR;
      RtR = (ctr->R).dot(ctr->R);
      RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Lower>() * ZtR);
      ratio = stepMhr + (mhr.logVThetaChol - mhr0.logVThetaChol) -
                  (0.5 * (ctr->n + 1.0) *
                  (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
                  log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
                  (log(treevar) * 0.5 * (mhr.nTerm - mhr0.nTerm));
    }
    
    if (log(R::runif(0, 1)) < ratio) {
      mhr0 = mhr;
      success = 2;
      tree->accept();
      dlnmTerm = tree->listTerminal();
      if (!ctr->binomial && !(ctr->zinb)) {
        tree->nodevals->tempV.resize(dlnmTerm.size(), dlnmTerm.size());
        tree->nodevals->tempV = mhr0.tempV;
      }
    }
  }
  if (success < 2)
    tree->reject();
    
  // Update variance and residuals
  if (ctr->shrinkage > 0)
    rHalfCauchyFC(&(ctr->tau(t)), mhr0.nTerm, 
                mhr0.termT2 / (ctr->sigma2 * ctr->nu));

  if ((ctr->tau)(t) != (ctr->tau)(t)) {
    // Rcout << ctr->gamma << "\n" << ctr->sigma2 << " " << ctr->tau << "\n" << ctr->Omega.mean();
    stop("\nNaN values occured during model run, rerun model.\n");
  }
  
  ctr->nTerm(t) = mhr0.nTerm;
  ctr->totTerm += mhr0.nTerm;
  ctr->sumTermT2 += mhr0.termT2 / ctr->tau(t);
  ctr->Rmat.col(t) = mhr0.Xd * mhr0.draw;

  // Record
  if (ctr->record > 0) {
    VectorXd rec(8);
    rec << ctr->record, t, (dlnmTerm[0]->nodestruct)->get(1),
    (dlnmTerm[0]->nodestruct)->get(2), (dlnmTerm[0]->nodestruct)->get(3),
    (dlnmTerm[0]->nodestruct)->get(4), mhr0.draw(0), 0;
    (dgn->DLMexp).push_back(rec);

    for (s = 1; s < dlnmTerm.size(); ++s) {
      rec[2] = (dlnmTerm[s]->nodestruct)->get(1);
      rec[3] = (dlnmTerm[s]->nodestruct)->get(2);
      rec[4] = (dlnmTerm[s]->nodestruct)->get(3);
      rec[5] = (dlnmTerm[s]->nodestruct)->get(4);
      rec[6] = mhr0.draw[s];
      (dgn->DLMexp).push_back(rec);
    }

    if (ctr->diagnostics) {
      VectorXd acc(5);
      acc << step, success, dlnmTerm.size(), stepMhr, ratio;
      (dgn->TreeAccept).push_back(acc);
    }
  }
} // end tdlnmGaussianTreeMCMC


//' dlmtree model with tdlnm approach
//'
//' @param model A list of parameter and data contained for the model fitting
//' @return A list of dlmtree model fit, mainly posterior mcmc samples
//' @export
// [[Rcpp::export]]
Rcpp::List tdlnm_Cpp(const Rcpp::List model)
{

  // #if defined(_OPENMP)
  //   if (as<int>(model["maxThreads"]) > 0) {
  //     omp_set_num_threads(as<int>(model["maxThreads"]));
  //     Eigen::setNbThreads(as<int>(model["maxThreads"]));
  //   } else {
  //     omp_set_num_threads(int(double(omp_get_max_threads()) / 2.0));
  //     Eigen::setNbThreads(int(double(omp_get_max_threads()) / 2.0));
  //   }
    
  // #endif
  
  // * Set up model control
  tdlmCtr *ctr = new tdlmCtr;
  ctr->iter = as<int>(model["nIter"]);
  ctr->burn = as<int>(model["nBurn"]);
  ctr->thin = as<int>(model["nThin"]);
  ctr->nRec = floor(ctr->iter / ctr->thin);
  ctr->nTrees = as<int>(model["nTrees"]);
  ctr->verbose = as<bool>(model["verbose"]);
  ctr->diagnostics = as<bool>(model["diagnostics"]);

  ctr->binomial = as<bool>(model["binomial"]);
  ctr->zinb = as<bool>(model["zinb"]); 
  ctr->stepProb = as<std::vector<double> >(model["stepProbTDLM"]);
  ctr->treePrior = as<std::vector<double> >(model["treePriorTDLM"]);
  ctr->shrinkage = as<int>(model["shrinkage"]);
  ctr->modKappa = 1.0;
  

  // * Set up model data
  ctr->Y0 = as<VectorXd>(model["Y"]);    
  ctr->Ystar = as<VectorXd>(model["Y"]);
  ctr->n = (ctr->Y0).size();                

  ctr->Z = as<MatrixXd>(model["Z"]);
  ctr->Zw = ctr->Z;
  ctr->pZ = (ctr->Z).cols();

  // ZI model
  ctr->Z1 = as<Eigen::MatrixXd>(model["Z.zi"]);  
  ctr->Zw1 = ctr->Z1; 
  ctr->pZ1 = (ctr->Z1).cols(); 


  // V_gamma for ZI model
  Eigen::MatrixXd VgInv1(ctr->pZ1, ctr->pZ1);    
  VgInv1 = (ctr->Z1).transpose() * (ctr->Z1);    
  VgInv1.diagonal().array() += 1.0 / 10000.0;      
  ctr->Vg1 = VgInv1.inverse();
  VgInv1.resize(0,0);    
  ctr->VgChol1 = (ctr->Vg1).llt().matrixL();

  // V_gamma for (Gaussian / Binary / NB model) 
  Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ);
  VgInv = (ctr->Z).transpose() * (ctr->Z);
  VgInv.diagonal().array() += 1.0 / 10000.0;
  ctr->Vg = VgInv.inverse();
  VgInv.resize(0,0);
  ctr->VgChol = (ctr->Vg).llt().matrixL();    


  // *** Set up parameters for logistic model ***
  ctr->binomialSize.resize(ctr->n);               ctr->binomialSize.setZero();
  ctr->kappa.resize(ctr->n);                      ctr->kappa.setOnes();
  ctr->Omega.resize(ctr->n);                      ctr->Omega.setOnes();
  if (ctr->binomial) {
    ctr->binomialSize = as<VectorXd>(model["binomialSize"]);
    ctr->kappa = ctr->Y0 - 0.5 * (ctr->binomialSize); 
    ctr->Ystar = ctr->kappa; 
    ctr->Omega.resize(ctr->n);                    ctr->Omega.setOnes();
  }


  // *** Set up parameters for ZINB model ***
  ctr->r = 5; 

  // Useful vectors for PG variable sampling 
  ctr->ones.resize(ctr->n);     ctr->ones.setOnes(); 
  ctr->rVec = (ctr->r) * (ctr->ones).array(); 

  // If ZINB, calculate the fixed values
  if(ctr->zinb){                                      
    ctr->z2 = 0.5*((ctr->Y0).array() - ctr->r).array(); 
    ctr->Ystar = ctr->z2;
  }

  // Initialize parameters
  ctr->w.resize(ctr->n);  
  ctr->b1 = as<VectorXd>(rnorm(ctr->pZ1, 0, sqrt(100))); 
  ctr->b2 = as<VectorXd>(rnorm(ctr->pZ, 0, sqrt(100)));  

  ctr->omega1.resize(ctr->n);        ctr->omega1.setOnes();
  ctr->omega2.resize(ctr->n);        ctr->omega2.setOnes();
  ctr->z1.resize(ctr->n);            ctr->z1.setZero();

  // Store the indices of y = 0
  for (int j = 0; j < ctr->n; ++j){
    if((ctr->Y0)[j] == 0){
      (ctr->yZeroIdx).push_back(j);
      (ctr->w)[j] = 0.5;
    } else {
      (ctr->NBidx).push_back(j);
      (ctr->w)[j] = 0;
    }
  }

  // NB model specific parameters
  ctr->Zstar = (ctr->Z).array().colwise() * (1 - ctr->w.array());
  ctr->yZeroN = (ctr->yZeroIdx).size(); 
  ctr->nStar = (ctr->NBidx).size();

  // * Create exposure data management
  exposureDat *Exp;
  if (as<int>(model["nSplits"]) == 0) { // DLM
    if (ctr->binomial || ctr->zinb)
      Exp = new exposureDat(as<MatrixXd>(model["Tcalc"]));
    else
      Exp = new exposureDat(as<MatrixXd>(model["Tcalc"]),
                            ctr->Z, ctr->Vg);
  } else { // DLNM
    if (ctr->binomial || ctr->zinb)
      Exp = new exposureDat(as<MatrixXd>(model["X"]),
                            as<MatrixXd>(model["SE"]),
                            as<VectorXd>(model["Xsplits"]),
                            as<MatrixXd>(model["Xcalc"]),
                            as<MatrixXd>(model["Tcalc"]),
                            as<bool>(model["lowmem"]));
    else
      Exp = new exposureDat(as<MatrixXd>(model["X"]),
                            as<MatrixXd>(model["SE"]),
                            as<VectorXd>(model["Xsplits"]),
                            as<MatrixXd>(model["Xcalc"]),
                            as<MatrixXd>(model["Tcalc"]),
                            ctr->Z, ctr->Vg,
                            as<bool>(model["lowmem"]));
  }
  ctr->pX = Exp->pX;
  ctr->nSplits = Exp->nSplits;

  // * Calculations used in special case: single-node trees
  ctr->X1 = (Exp->Tcalc).col(ctr->pX - 1);
  ctr->ZtX1 = (ctr->Z).transpose() * (ctr->X1);
  ctr->VgZtX1 = (ctr->Vg).selfadjointView<Lower>() * (ctr->ZtX1);
  ctr->VTheta1Inv = (ctr->X1).dot(ctr->X1) - (ctr->ZtX1).dot(ctr->VgZtX1);

  // * Create trees
  int t;
  std::vector<Node*> trees;
  NodeStruct *ns;
  ns = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                      as<VectorXd>(model["splitProb"]),
                      as<VectorXd>(model["timeSplits0"]));
  for (t = 0; t < ctr->nTrees; ++t) {
    trees.push_back(new Node(0, 1));
    trees[t]->nodestruct = ns->clone();
    Exp->updateNodeVals(trees[t]);
  }
  delete ns;
  ctr->nTerm.resize(ctr->nTrees);                   ctr->nTerm.setOnes();
  ctr->Rmat.resize(ctr->n, ctr->nTrees);            ctr->Rmat.setZero();

  // * Setup model logs
  tdlmLog *dgn = new tdlmLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec);          (dgn->gamma).setZero();
  (dgn->sigma2).resize(ctr->nRec);                  (dgn->sigma2).setZero();
  (dgn->nu).resize(ctr->nRec);                      (dgn->nu).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec);        (dgn->tau).setZero();
  (dgn->fhat).resize(ctr->n);                       (dgn->fhat).setZero();
  (dgn->termNodes).resize(ctr->nTrees, ctr->nRec);  (dgn->termNodes).setZero();
  dgn->timeProbs.resize(ctr->pX - 1, ctr->nRec);    (dgn->timeProbs).setZero();
  VectorXd Yhat(ctr->n);                             Yhat.setZero();

  // ZINB specific log
  (dgn->b1).resize(ctr->pZ1, ctr->nRec);             (dgn->b1).setZero(); 
  (dgn->b2).resize(ctr->pZ, ctr->nRec);              (dgn->b2).setZero(); 
  (dgn->r).resize(ctr->nRec);                        (dgn->r).setZero(); 
  (dgn->wMat).resize(ctr->n, ctr->nRec);             (dgn->wMat).setZero(); 
  
  // * Initial values and draws
  ctr->fhat.resize(ctr->n);                         (ctr->fhat).setZero();
  ctr->R = ctr->Ystar; 
  ctr->gamma.resize(ctr->pZ);
  // Load initial params for faster convergence in binomial model
  if (ctr->binomial) {
    ctr->gamma = as<VectorXd>(model["initParams"]);
    ctr->Omega = rcpp_pgdraw(ctr->binomialSize, ctr->fhat + ctr->Z * ctr->gamma);
    ctr->Zw = ctr->Omega.asDiagonal() * ctr->Z;
    ctr->VgInv = ctr->Z.transpose() * ctr->Zw;
    ctr->VgInv.diagonal().array() += 1 / 1000.0;
    ctr->Vg = ctr->VgInv.inverse();
    ctr->VgChol = ctr->Vg.llt().matrixL();
    // recalculate 'pseudo-Y' = kappa / omega, kappa = (y - n_b)/2
    ctr->Ystar = ctr->kappa.array() / ctr->Omega.array();
  }
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  ctr->nu = 1.0; 
  ctr->sigma2 = 1.0;

  tdlmModelEst(ctr);

  rHalfCauchyFC(&(ctr->nu), ctr->nTrees, 0.0);
  (ctr->tau).resize(ctr->nTrees);                   (ctr->tau).setOnes();
  if (ctr->shrinkage > 0) {
    for (t = 0; t < ctr->nTrees; ++t) {
      rHalfCauchyFC(&(ctr->tau(t)), 0.0, 0.0);
    }
  }
  
  // * Create progress meter
  progressMeter* prog = new progressMeter(ctr);

  // * Beginning of MCMC
  std::size_t s;
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    Rcpp::checkUserInterrupt();

    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }

    // * Update trees
    ctr->R += (ctr->Rmat).col(0);
    (ctr->fhat).setZero();
    ctr->totTerm = 0.0; 
    ctr->sumTermT2 = 0.0;
    for (t = 0; t < ctr->nTrees; ++t) {
      tdlnmTreeMCMC(t, trees[t], ctr, dgn, Exp);
      ctr->fhat += (ctr->Rmat).col(t);
      if (t < ctr->nTrees - 1) {
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
      }
    } // end update trees

    // * Update model
    ctr->R = ctr->Ystar - ctr->fhat; 
    tdlmModelEst(ctr);
    rHalfCauchyFC(&(ctr->nu), ctr->totTerm, ctr->sumTermT2 / ctr->sigma2);
    if ((ctr->sigma2 != ctr->sigma2) || (ctr->nu != ctr->nu)) {
      // Rcout << ctr->gamma << "\n" << ctr->sigma2 << " " << ctr->nu << "\n" << ctr->Omega.mean() << " " << ctr->Y.mean();
      stop("\nNaN values occured during model run, rerun model.\n");
    }

    // * Record
    if (ctr->record > 0) {
      (dgn->gamma).col(ctr->record - 1) = ctr->gamma;
      (dgn->sigma2)(ctr->record - 1) = ctr->sigma2;
      (dgn->nu)(ctr->record - 1) = ctr->nu;
      (dgn->tau).col(ctr->record - 1) = ctr->tau;
      (dgn->termNodes).col(ctr->record - 1) = ctr->nTerm;
      dgn->timeProbs.col(ctr->record -1) = trees[0]->nodestruct->getTimeProbs();
      dgn->fhat += ctr->fhat;
      Yhat += ctr->fhat + ctr->Z * ctr->gamma;

      // ZINB
      (dgn->b1).col(ctr->record - 1) = ctr->b1;
      (dgn->b2).col(ctr->record - 1) = ctr->b2;
      (dgn->r)(ctr->record - 1) = ctr->r;
      (dgn->wMat).col(ctr->record - 1) = ctr->w;
    }
    // * Update progress
    prog->printMark();
  } // end MCMC

  // * Setup data for return
  MatrixXd DLM((dgn->DLMexp).size(), 8);
  for (s = 0; s < (dgn->DLMexp).size(); ++s)
    DLM.row(s) = dgn->DLMexp[s];
  VectorXd sigma2 = dgn->sigma2;
  VectorXd nu = dgn->nu;
  VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  MatrixXd gamma = (dgn->gamma).transpose();
  MatrixXd tau = (dgn->tau).transpose();
  MatrixXd termNodes = (dgn->termNodes).transpose();
  MatrixXd timeProbs = (dgn->timeProbs).transpose();
  VectorXd YhatOut = Yhat / ctr->nRec;
  MatrixXd Accept((dgn->TreeAccept).size(), 5);
  for (s = 0; s < (dgn->TreeAccept).size(); ++s){
    Accept.row(s) = dgn->TreeAccept[s];
  }

  // ZINB specific return 
  Eigen::MatrixXd b1 = (dgn->b1).transpose(); 
  Eigen::MatrixXd b2 = (dgn->b2).transpose();
  Eigen::VectorXd r = dgn->r; 
  Eigen::MatrixXd wMat = dgn->wMat; 

  delete prog;
  // delete ctr;
  delete dgn;
  delete Exp;
  for (s = 0; s < trees.size(); ++s)
    delete trees[s];

  return(Rcpp::List::create(Named("TreeStructs")  = wrap(DLM),
                            Named("fhat")         = wrap(fhat),
                            Named("Yhat")         = wrap(YhatOut),
                            Named("sigma2")       = wrap(sigma2),
                            Named("nu")           = wrap(nu),
                            Named("tau")          = wrap(tau),
                            Named("timeProbs")    = wrap(timeProbs),
                            Named("termNodes")    = wrap(termNodes),
                            Named("gamma")        = wrap(gamma),
                            Named("treeAccept")   = wrap(Accept),
                            Named("b1")           = wrap(b1),
                            Named("b2")           = wrap(b2),
                            Named("r")            = wrap(r),
                            Named("wMat")         = wrap(wMat)));
} // end tdlnm_Cpp
