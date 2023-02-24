/**
 * @file tdlnmGaussian.cpp
 * @author Daniel Mork (danielmork.github.io)
 * @brief Method for Treed DLNM and Treed DLM
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
    MatrixXd ZtX(ctr->pZ, pX);    ZtX.setZero();
    MatrixXd VgZtX(ctr->pZ, pX);  VgZtX.setZero();

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
    VectorXd XtVzInvR(ctr->n);
    
    if (ctr->binomial) {
      const MatrixXd Xdw = (ctr->Omega).asDiagonal() * out.Xd;
      tempV = Xdw.transpose() * out.Xd;
      tempV.noalias() -= ZtX.transpose() * VgZtX;
      XtVzInvR = Xdw.transpose() * ctr->R;
      
    } else if (ctr->zinb){ // ZINB (subsetting at-risk observations)
    Eigen::MatrixXd outXdstar = selectIndM(out.Xd, ctr->atRiskIdx);        // At-risk observations only
    const Eigen::MatrixXd Xdw = (selectIndM(ctr->omega2, ctr->atRiskIdx)).asDiagonal() * outXdstar; // (n*xn*) x (n*x2) = (n*x2) 
    tempV = Xdw.transpose() * outXdstar;                                   // (2xn*) x (n*x2) = (2x2)
    tempV.noalias() -= ZtX.transpose() * VgZtX;                            // (2x2) - (2x5) x (5x2)
    XtVzInvR = Xdw.transpose() * selectIndM(ctr->R, ctr->atRiskIdx);       // (2xn*) x (n*x1) = (2x1)

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
    out.draw.noalias() += 
      VThetaChol * as<VectorXd>(rnorm(pX, 0, sqrt(ctr->sigma2)));
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
  int step =        0;
  int success =     0;
  double stepMhr =  0.0;
  double ratio =    0.0;
  double treevar =  ctr->nu * ctr->tau(t);
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
  if (ctr->shrinkage)
    rHalfCauchyFC(&(ctr->tau(t)), mhr0.nTerm, 
                mhr0.termT2 / (ctr->sigma2 * ctr->nu));

  if ((ctr->tau)(t) != (ctr->tau)(t)) 
    stop("\nNaN values occured during model run, rerun model.\n");

  ctr->nTerm(t) = mhr0.nTerm;
  ctr->totTerm += mhr0.nTerm;
  ctr->sumTermT2 += mhr0.termT2 / ((ctr->tau)(t));
  ctr->Rmat.col(t) = mhr0.Xd * mhr0.draw;

  // Record
  if (ctr->record > 0) {
    VectorXd rec(8);
    rec << ctr->record, t, (dlnmTerm[0]->nodestruct)->get(1),
    (dlnmTerm[0]->nodestruct)->get(2), (dlnmTerm[0]->nodestruct)->get(3),
    (dlnmTerm[0]->nodestruct)->get(4), mhr0.draw(0), 0;
    (dgn->DLMexp).push_back(rec);

    for(s = 1; s < dlnmTerm.size(); ++s) {
      rec[2] = (dlnmTerm[s]->nodestruct)->get(1);
      rec[3] = (dlnmTerm[s]->nodestruct)->get(2);
      rec[4] = (dlnmTerm[s]->nodestruct)->get(3);
      rec[5] = (dlnmTerm[s]->nodestruct)->get(4);
      rec[6] = mhr0.draw(s);
      (dgn->DLMexp).push_back(rec);
    }

    if (ctr->diagnostics) {
      VectorXd acc(5);
      acc << step, success, dlnmTerm.size(), stepMhr, ratio;
      (dgn->TreeAccept).push_back(acc);
    }
  }
} // end tdlnmGaussianTreeMCMC

/**
 * @brief tdlnm
 * 
 * @param model list of data and model control specs from R 
 * @return Rcpp::List 
 */
// [[Rcpp::export]]
Rcpp::List tdlnm_Cpp(const Rcpp::List model)
{
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
  ctr->zinb = as<bool>(model["zinb"]);          // ZINB
  ctr->stepProb = as<std::vector<double> >(model["stepProb"]);
  ctr->treePrior = as<std::vector<double> >(model["treePrior"]);
  ctr->shrinkage = as<int>(model["shrinkage"]);
  
  // * Set up model data
  ctr->Y0 = as<VectorXd>(model["Y"]); // Change to Y0 later for the consistency later 
  ctr->n = (ctr->Y0).size(); // Change to Y0 later for the consistency later

  ctr->Z = as<MatrixXd>(model["Z"]);
  ctr->Zw = ctr->Z;
  ctr->pZ = (ctr->Z).cols();

    // binary component
  ctr->Z1 = as<Eigen::MatrixXd>(model["Z_b"]);       // Fixed effect design matrix, Z (n x p)
  ctr->Zw1 = ctr->Z1;                               // Copy Z to Zw, will be updated throughout iteration (Z Omega)
  ctr->pZ1 = (ctr->Z1).cols();                      // Number of fixed effect covariates (ncol(Z))

  // Full conditional initialization for logistic model: V_gamma
  // Compute V_gamma^{-1} = Zt*Z + I/c where c = 100000 from Supplemental Eq(2) of TDLM
  Eigen::MatrixXd VgInv1(ctr->pZ1, ctr->pZ1);        // pxp var-cov matrix
  VgInv1 = (ctr->Z1).transpose() * (ctr->Z1);        // Zt*Z
  VgInv1.diagonal().array() += 1.0 / 100.0;     // Zt*Z + I/c where c = 100000
  ctr->Vg1 = VgInv1.inverse();                      // V_gamma
  VgInv1.resize(0,0);                              // Clear VgInv now that we obtained Vg
  ctr->VgChol1 = (ctr->Vg1).llt().matrixL();        // Compute the cholesky of Vg: V_gamma = L*Lt = VgChol * t(VgChol)
                                                  // (Ref: https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html)

  Eigen::MatrixXd VgInv = (ctr->Z).transpose() * (ctr->Z);
  VgInv.diagonal().array() += 1.0 / 100.0;
  ctr->Vg = VgInv.inverse();
  VgInv.resize(0,0);
  ctr->VgChol = (ctr->Vg).llt().matrixL();
  
  
  // * Set up data for logistic model
  ctr->binomialSize.resize(ctr->n);               ctr->binomialSize.setZero();
  ctr->kappa.resize(ctr->n);                      ctr->kappa.setOnes();
  ctr->Omega.resize(ctr->n);                      ctr->Omega.setOnes();
  if (ctr->binomial) {
    ctr->binomialSize = as<VectorXd>(model["binomialSize"]);
    ctr->kappa = ctr->Y0 - 0.5 * (ctr->binomialSize); // Change to Y0 later for the consistency later
    ctr->Ystar = ctr->kappa; // Change to Ystar later for the consistency later
    ctr->Omega.resize(ctr->n);                      ctr->Omega.setOnes();
  }

  // ============ Set up parameters for ZINB model ============
  // Initialize values
  ctr->r = 5;                                       // Dispersion parameter: Starting at the center of Unif(0, 10)

  // Useful vectors for PG variable sampling 
  ctr->ones.resize(ctr->n);   ctr->ones.setOnes();  // Vector of ones
  ctr->rVec = (ctr->r) * (ctr->ones).array();       // Dispersion parameter as a vector: rep(r, n)

  // If ZINB, calculate the fixed values
  if(ctr->zinb){                                      
    ctr->z2 = 0.5*((ctr->Y0).array() - ctr->r).array();   // Compute z2
    ctr->Ystar = ctr->z2;                                 // Ystar in tdlmm_Cpp.cpp, z2 in modelEst.cpp
  }

  ctr->w.resize(ctr->n);                                        // At-risk binary latent variable
  ctr->b1 = as<Eigen::VectorXd>(rnorm(ctr->pZ1, 0, sqrt(100)));   // Prior sampling of coefficients for binary component (p x 1)
  ctr->b2 = as<Eigen::VectorXd>(rnorm(ctr->pZ, 0, sqrt(100)));   // Prior sampling of coefficients for NegBin component (p x 1)

  ctr->omega1.resize(ctr->n);        ctr->omega1.setOnes();      // Initiate with omega1(binary) as an identity matrix
  ctr->omega2.resize(ctr->n);        ctr->omega2.setOnes();      // Initiate with omega2(NegBin) as an identity matrix
  
  ctr->z1.resize(ctr->n);            ctr->z1.setZero();          // z1 for binary component

  // Store the indices of y == 0 (yZeroIdx) & y != 0 (atRiskIdx)
  for(int j = 0; j < ctr->n; j++){
    if((ctr->Y0)[j] == 0){
      (ctr->yZeroIdx).push_back(j);     // y == 0: Could be structural or at-risk zero
      (ctr->w)[j] = 0.5;
    } else {
      (ctr->atRiskIdx).push_back(j);    // y != 0: Determined to be at-risk as y != 0
      (ctr->w)[j] = 1;
    }
  }

  ctr->Zstar = (ctr->Z).array().colwise() * (ctr->w).array();    // Fixed effect matrix with [w == 0] zeroed out.

  // Calculate the size: yZeroN + nStar should equal to n
  ctr->yZeroN = (ctr->yZeroIdx).size();   // Fixed as we are just counting y == 0
  ctr->nStar = (ctr->atRiskIdx).size();   // Random as some of y == 0 can be at-risk

  // ============ Spatial information ============
  // Spatial flag
  ctr->spatial = as<bool>(model["spatial"]);

  // Spatial model control parameters (Default is all zeros)
  // As a default, number of unique areas is set as the same as sample size.
  ctr->spN = ctr->n;

  // Binary assignment matrix: areaA
  (ctr->areaA).resize(ctr->n, ctr->spN);
  (ctr->areaA).setZero();

  // Phi of CAR model (spN x 1)
  ctr->spPhi.resize(ctr->spN);
  ctr->spPhi.setZero();

  // Random uniform distribution sampling to initialize spPhi
  std::default_random_engine Generator;
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);

  // If spatial random effect is called:
  if(ctr->spatial){
    // Spatial information is not null, hence store them to ctr.
    // spNodes1 & spNodes2
    ctr->spNodes1 = as<Eigen::VectorXd>(model["spNodes1"]);
    ctr->spNodes2 = as<Eigen::VectorXd>(model["spNodes2"]);

    ctr->areaD = as<Eigen::MatrixXd>(model["areaD"]);         // Diagonal neighbor # matrix
    ctr->areaW = as<Eigen::MatrixXd>(model["areaW"]);         // Adjacency matrix
    ctr->areaA = as<Eigen::MatrixXd>(model["areaA"]);         // Assignment matrix A so that A*phi is n x 1
    ctr->spN = as<int>(model["spN"]);                         // Number of unique areas

    ctr->spPhi.resize(ctr->spN);        // spN has been updated -> Update the vector size as well
    ctr->spPhi.setZero();               // Set zero.

    ctr->spTau = 1/R::rgamma(1.0, 1.0) ;  // variance: spTau prior ~ IG(1, 1)
    ctr->rho = 1;                         // Spatial correlation

    // Initialize spPhi ~ Unif(-2, 2)
    for (int spIndex = 0; spIndex < ctr->spN; spIndex++) {
      double number = distribution(Generator);
      ctr->spPhi[spIndex] = number;
    }

    // Initialize components for full conditional of spPhi, similarly to binary component of ZINB   
    ctr->zPhi = ctr->z1;
    
    // Similarly to V_gamma, parameterize as Vp for V_phi
    Eigen::MatrixXd VpInv(ctr->spN, ctr->spN);        // Vp (spN x spN)
    VpInv = ctr->areaA.transpose() * ctr->omega1.asDiagonal() * ctr->areaA;  //
    VpInv.noalias() += (1 / ctr->spTau) * (ctr->areaD - (ctr->rho) * ctr->areaW);
    ctr->Vp = VpInv.inverse();
    VpInv.resize(0,0);                                // Clear VpInv now that we obtained Vp
    ctr->VpChol = (ctr->Vp).llt().matrixL();          // Compute the cholesky of Vp: V_phi = L*Lt = VpChol * t(VpChol)
  }

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
                            as<MatrixXd>(model["Tcalc"]));
    else
      Exp = new exposureDat(as<MatrixXd>(model["X"]),
                            as<MatrixXd>(model["SE"]),
                            as<VectorXd>(model["Xsplits"]),
                            as<MatrixXd>(model["Xcalc"]),
                            as<MatrixXd>(model["Tcalc"]),
                            ctr->Z, ctr->Vg);
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
                      as<VectorXd>(model["timeProb"]));
  for (t = 0; t < ctr->nTrees; t++) {
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

  // ZINB log
  (dgn->b1).resize(ctr->pZ1, ctr->nRec);             (dgn->b1).setZero();         // Binary component coeffient (p x MCMC)
  (dgn->b2).resize(ctr->pZ, ctr->nRec);              (dgn->b2).setZero();         // NegBin component coeffient (p x MCMC)
  (dgn->r).resize(ctr->nRec);                        (dgn->r).setZero();          // Dispersion parameter
  (dgn->wMat).resize(ctr->n, ctr->nRec);             (dgn->wMat).setZero();       // At-risk Auxiliary: Each column is w at each iteration
  (dgn->spPhi).resize(ctr->spN, ctr->nRec);          (dgn->spPhi).setZero();
  (dgn->rho).resize(ctr->nRec);                      (dgn->rho).setZero();        // Spatial correlation
  (dgn->spTau).resize(ctr->nRec);                    (dgn->spTau).setZero();        // Spatial precision

  
  // * Initial values and draws
  ctr->fhat.resize(ctr->n);                         (ctr->fhat).setZero();
  ctr->R = ctr->Ystar; // Change to Y0 or Ystar later for the consistency later
  ctr->gamma.resize(ctr->pZ);
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  ctr->nu = 1.0; // ! Need to define nu and sigma2 prior to ModelEst
  ctr->sigma2 = 1.0;
  tdlmModelEst(ctr);
  rHalfCauchyFC(&(ctr->nu), ctr->nTrees, 0.0);
  (ctr->tau).resize(ctr->nTrees);                   (ctr->tau).setOnes();
  if (ctr->shrinkage) {
    for (t = 0; t < ctr->nTrees; t++) 
      rHalfCauchyFC(&(ctr->tau(t)), 0.0, 0.0);
  }
  
  // * Create progress meter
  progressMeter* prog = new progressMeter(ctr);

  // * Beginning of MCMC
  std::size_t s;
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }

    // * Update trees
    ctr->R += (ctr->Rmat).col(0);
    (ctr->fhat).setZero();
    ctr->totTerm = 0.0; ctr->sumTermT2 = 0.0;
    for (t = 0; t < ctr->nTrees; t++) {
      tdlnmTreeMCMC(t, trees[t], ctr, dgn, Exp);
      ctr->fhat += (ctr->Rmat).col(t);
      if (t < ctr->nTrees - 1)
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
    } // end update trees

    // * Update model
    ctr->R = ctr->Ystar - ctr->fhat; // Change to Y0 or Ystar later for the consistency later
    tdlmModelEst(ctr);
    rHalfCauchyFC(&(ctr->nu), ctr->totTerm, ctr->sumTermT2 / ctr->sigma2);
    if ((ctr->sigma2 != ctr->sigma2) || (ctr->nu != ctr->nu))
      stop("\nNaN values occured during model run, rerun model.\n");

    // * Record
    if (ctr->record > 0) {
      (dgn->gamma).col(ctr->record - 1) = ctr->gamma;
      (dgn->sigma2)(ctr->record - 1) = ctr->sigma2;
      (dgn->nu)(ctr->record - 1) = ctr->nu;
      (dgn->tau).col(ctr->record - 1) = ctr->tau;
      (dgn->termNodes).col(ctr->record - 1) = ctr->nTerm;
      dgn->fhat += ctr->fhat;

      // ZINB
      (dgn->b1).col(ctr->record - 1) = ctr->b1;
      (dgn->b2).col(ctr->record - 1) = ctr->b2;
      (dgn->r)(ctr->record - 1) = ctr->r;
      (dgn->wMat).col(ctr->record - 1) = ctr->w;
      (dgn->spPhi).col(ctr->record - 1) = ctr->spPhi;
      (dgn->rho)(ctr->record - 1) = ctr->rho;
      (dgn->spTau)(ctr->record - 1) = ctr->spTau;
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
  MatrixXd Accept((dgn->TreeAccept).size(), 5);
  for (s = 0; s < (dgn->TreeAccept).size(); ++s)
    Accept.row(s) = dgn->TreeAccept[s];

  // ZINB & Spatial return
  // ZINB return 
  Eigen::MatrixXd b1 = (dgn->b1).transpose(); // binary component coefficients
  Eigen::MatrixXd b2 = (dgn->b2).transpose(); // count component coefficients (This gets return for NB)
  Eigen::VectorXd r = dgn->r;                 // dispersion parameter
  Eigen::MatrixXd wMat = dgn->wMat;           // What percentage of iteration was an individual at risk?

  // Spatial return
  Eigen::VectorXd rho = dgn->rho;
  Eigen::VectorXd spTau = dgn->spTau;
  Eigen::MatrixXd spPhi = (dgn->spPhi).transpose();

  delete prog;
  delete ctr;
  delete dgn;
  delete Exp;
  for (s = 0; s < trees.size(); ++s)
    delete trees[s];

  return(Rcpp::List::create(Named("DLM")    = wrap(DLM),
                            Named("fhat")   = wrap(fhat),
                            Named("sigma2") = wrap(sigma2),
                            Named("nu")     = wrap(nu),
                            Named("tau")    = wrap(tau),
                            Named("termNodes")  = wrap(termNodes),
                            Named("gamma")  = wrap(gamma),
                            Named("treeAccept") = wrap(Accept),
                            Named("b1") = wrap(b1),
                            Named("b2") = wrap(b2),
                            Named("r") = wrap(r),
                            Named("wMat") = wrap(wMat),
                            Named("spPhi") = wrap(spPhi),
                            Named("spTau") = wrap(spTau),
                            Named("rho") = wrap(rho)));
} // end tdlnm_Cpp