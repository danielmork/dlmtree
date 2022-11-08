/**
 * @file monotdlnm_Cpp.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-03-06
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "RcppEigen.h"
#include "mvtnorm.h"
#include "modelCtr.h"
#include "exposureDat.h"
#include "Node.h"
#include "NodeStruct.h"
#include "Fncs.h"
using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Lower;


treeMHR monoDlnmMHR(std::vector<Node*> dlnmTerm, tdlmCtr* ctr, VectorXd ZtR, double treevar, Node* tree, bool updateNested)
{
  treeMHR out;
  int totTerm = 0;
  std::vector<std::vector<std::pair<int, Node*> > > nestedTerm;
  for (Node* eta : dlnmTerm) {
    // sort nodes in increasing xmin order
    std::vector<std::pair<int, Node*> > nodeList;
    for (Node* lambda : eta->nodevals->nestedTree->listTerminal(updateNested))
      nodeList.push_back(std::make_pair(lambda->nodestruct->get(1), lambda));
    sort(nodeList.begin(), nodeList.end());
    
    totTerm += nodeList.size() - 1;
    nestedTerm.push_back(nodeList);
  }
    
  if (totTerm == 0) {
    out.Xd.resize(ctr->n, 1); out.Xd.setZero();
    out.Dtrans.resize(1, 1);  out.Dtrans.setOnes();
    out.draw.resize(1);       out.draw.setZero();      
    out.cdf =                 1.0;
    out.beta =                0.0;      
    out.logVThetaChol =       0.0;
    out.termT2 =              0.0;
    out.totTerm =             0.0;
    return(out);
  }
  
  out.Xd.resize(ctr->n, totTerm);   out.Xd.setZero();
  MatrixXd Dinv(totTerm, totTerm);  Dinv.setIdentity();
  MatrixXd ZtX(ctr->pZ, totTerm);   ZtX.setZero();
  MatrixXd VgZtX(ctr->pZ, totTerm); VgZtX.setZero();
  
  // * Define transformation matrices (D and Dinv), fill U, ZtX, VgZtX
  int k = 0;  // index for D, Dinv
  int i = 0;  // index for X matrix
  int l;      // size of block D matrix
  for (std::vector<std::pair<int, Node*> > nt : nestedTerm) {
    if (nt.size() > 1) {
      l = nt.size() - 1;
      Dinv.block(k, k, l, l).triangularView<Lower>().setOnes();
      k += nt.size() - 1;
      for (std::pair<int, Node*> nt2 : nt) {
        if (nt2.first > 0) {
          out.Xd.col(i) = nt2.second->nodevals->X;
          if (ctr->binomial) {
            ZtX.col(i) = ctr->Zw.transpose() * nt2.second->nodevals->X;
            VgZtX.col(i) = ctr->Vg * ZtX.col(i);
            
          } else {
            ZtX.col(i) = nt2.second->nodevals->ZtX;
            VgZtX.col(i) = nt2.second->nodevals->VgZtX;
          }
          ++i;
        } // end if xmin > 0
      } // end loop over nested terminal nodes
    } // end if > 1 nested terminal nodes
  } // end loop over vector of Node*
  
  out.Xd =  out.Xd * Dinv;
  ZtX =     ZtX * Dinv;
  VgZtX =   VgZtX * Dinv;  
  
  // * Calculate covariance matrix V_theta
  MatrixXd tempV(totTerm, totTerm);
  VectorXd XtVzInvR(ctr->n);
  
  if (ctr->binomial) {
    const MatrixXd Xdw = ctr->Omega.asDiagonal() * out.Xd;
    tempV.triangularView<Lower>() = Xdw.transpose() * out.Xd;
    tempV.noalias() -= ZtX.transpose() * VgZtX;
    XtVzInvR = Xdw.transpose() * ctr->R;
    
  } else {
    if (updateNested) {
      tempV.triangularView<Lower>() = out.Xd.transpose() * out.Xd;
      tempV.noalias() -= ZtX.transpose() * VgZtX;
      out.tempV = tempV;
    } else {
      tempV = tree->nodevals->tempV;
    }
    XtVzInvR = out.Xd.transpose() * ctr->R;
  }
  
  XtVzInvR.noalias() -=        VgZtX.transpose() * ZtR;
  tempV.diagonal().array() +=  1.0 / treevar;
  const MatrixXd VTheta =      tempV.selfadjointView<Lower>().llt().
                                solve(MatrixXd::Identity(totTerm, totTerm));
  const MatrixXd VThetaChol =  VTheta.llt().matrixL();
  const VectorXd ThetaHat =    VTheta * XtVzInvR;
  
  out.Dtrans =        Dinv;   
  if (ctr->debug) 
    Rcout << "\n Term = " << totTerm << "\n That = \n" << ThetaHat << "\n VTheta = \n" << ctr->sigma2 * VTheta << "\n Xcounts = \n" <<
    out.Xd.colwise().sum() << "\n";
  out.draw =          rtmvnorm(ThetaHat, ctr->sigma2 * VTheta, 1);    
  if (ctr->debug)
    Rcout << "\n draw = " << out.draw;
  out.cdf =           zeroToInfNormCDF(ThetaHat, ctr->sigma2 * VTheta);  
  if (ctr->debug)
    Rcout << "\n cdf = " << out.cdf;  
  out.beta =          ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
  out.termT2 =        out.draw.dot(out.draw);
  out.totTerm =       (double) totTerm;
  return(out);
}




/**
 * @brief update tree function
 * 
 * @param t tree number
 * @param tree tree pointer
 * @param ctr modelCtr pointer
 * @param dgn log pointer
 * @param Exp exposureDat pointer
 * @param nsX pointer to exposure concentration struct
 */
void monoTDLNMTreeUpdate(int t, Node* tree, tdlmCtr* ctr, tdlmLog* dgn, exposureDat* Exp, NodeStruct* nsX)
{
  int step =        0;
  int success =     0;
  double stepMhr =  0.0;
  double ratio =    0.0;
  double treevar =  ctr->nu * ctr->tau(t);
  std::vector<Node*> dlnmTerm, newDlnmTerm, nestedTerm;
  Node* nestedTree;
  treeMHR mhr0, mhr;
  
  // List current tree terminal nodes
  dlnmTerm = tree->listTerminal();
  VectorXd ZtR = ctr->Zw.transpose() * ctr->R;
  if (dlnmTerm.size() > 1) // if single terminal node, grow is only option
    step = sampleInt(ctr->stepProb, 1);
  
  // propose update
  stepMhr = tdlmProposeTree(tree, 0, ctr, step);
  success = tree->isProposed();
  
  if (ctr->debug)
    Rcout << "\nstep = " << step << " nnodes = " << dlnmTerm.size() << " stepMhr = " << stepMhr;
  
  if (tree->nodevals->tempV.rows() > 0)
    mhr0 = monoDlnmMHR(dlnmTerm, ctr, ZtR, treevar, tree, 0);
  else
    mhr0 = monoDlnmMHR(dlnmTerm, ctr, ZtR, treevar, tree, 1);
    
  if (success && (stepMhr == stepMhr)) {
    newDlnmTerm = tree->listTerminal(1);
    
    // draw nested trees if grow/prune
    if (step < 2) {
      for (Node* eta : newDlnmTerm) {
        if (ctr->debug) {
          Rcout << "\n"; eta->nodestruct->printStruct();
        }
        if (eta->nodevals == 0)
          eta->nodevals = new NodeVals(ctr->n);
        if (eta->nodevals->nestedTree == 0) {
          drawZirt(eta, ctr, nsX);
          
          for (Node* lambda : eta->nodevals->nestedTree->listTerminal()) {
            if (ctr->debug) {
              Rcout << "\n\t"; lambda->nodestruct->printStruct();
            }
            Exp->updateNodeVals(lambda);
          }
        }
      }
    } else { // end draw nested trees if grow/prune
      for (Node* eta : newDlnmTerm) { // update time range for nested nodes
        int tmin = eta->nodestruct->get(3);
        int tmax = eta->nodestruct->get(4);
        if (eta->nodevals == 0)
          eta->nodevals = new NodeVals(ctr->n);
        if ((tmin != eta->nodevals->nestedTree->nodestruct->get(3)) ||
            (tmax != eta->nodevals->nestedTree->nodestruct->get(4))) {
          eta->nodevals->nestedTree->nodestruct->setTimeRange(tmin, tmax);
          eta->nodevals->nestedTree->updateStruct();
          eta->nodevals->nestedTree->setUpdate(1);
          if (ctr->debug) {
            Rcout << "\n"; eta->nodestruct->printStruct();
          }
          for (Node* lambda : eta->nodevals->nestedTree->listTerminal()) {
            if (ctr->debug) {
              Rcout << "\n\t"; lambda->nodestruct->printStruct();
            }
            Exp->updateNodeVals(lambda);
          }
        }
      }
    }
    
    mhr = monoDlnmMHR(newDlnmTerm, ctr, ZtR, treevar, tree, 1);
    
    // combine mhr parts into log-MH ratio
    ratio = stepMhr + (mhr.logVThetaChol - mhr0.logVThetaChol) +
      (0.5 * (mhr.beta - mhr0.beta) * (1 / ctr->sigma2)) -
      (log(4 * ctr->sigma2 * treevar) * 0.5 * (mhr.totTerm - mhr0.totTerm)) +
      log(mhr.cdf) - log(mhr0.cdf);
    if (ctr->debug)
      Rcout << " ratio = " << ratio;
    
    if (log(R::runif(0, 1)) < ratio) {
      mhr0 = mhr;
      success = 2;
      tree->accept();
      dlnmTerm = tree->listTerminal();
      if (!(ctr->binomial))
        tree->nodevals->tempV = mhr0.tempV;
    }
  }
  if (success < 2) {
    tree->reject();
    for (Node* eta : dlnmTerm) { // update time range for nested nodes
      int tmin = eta->nodestruct->get(3);
      int tmax = eta->nodestruct->get(4);
      if ((tmin != eta->nodevals->nestedTree->nodestruct->get(3)) ||
          (tmax != eta->nodevals->nestedTree->nodestruct->get(4))) {
        eta->nodevals->nestedTree->nodestruct->setTimeRange(tmin, tmax);
        eta->nodevals->nestedTree->updateStruct();
        eta->nodevals->nestedTree->setUpdate(1);
        for (Node* lambda : eta->nodevals->nestedTree->listTerminal())
          Exp->updateNodeVals(lambda);
      }
    }
  }
  
  // * Propose new nested tree at each terminal node
  for (Node* eta : dlnmTerm) {
    int tmin = eta->nodestruct->get(3);
    int tmax = eta->nodestruct->get(4);
    nestedTree = eta->nodevals->nestedTree;
    nestedTerm = nestedTree->listTerminal();
    
    step = 0;
    if (nestedTerm.size() > 1)
      step = sampleInt(ctr->stepProb, 1);
    
    if (ctr->debug)
      Rcout << "\n\t nt step = " << step;
      
    stepMhr = tdlmProposeTree(nestedTree, Exp, ctr, step, 0.0);
    success = nestedTree->isProposed();
    
    if (success && (stepMhr == stepMhr)) {
      // grow from depth 0
      if ((step == 0) && (nestedTerm.size() == 1)) {
        stepMhr = logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 0) -
          logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 1) + 
          2 * logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], 1.0, 1);
      // prune from depth 1
      } else if ((step == 1) && (nestedTerm.size() == 2)) {
        stepMhr = logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 1) -
          logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 0) -
          2 * logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], 1.0, 1);
      }
      
      // calculate MHR
      mhr = monoDlnmMHR(dlnmTerm, ctr, ZtR, treevar, tree, 1);
      ratio = stepMhr + (mhr.logVThetaChol - mhr0.logVThetaChol) +
        (0.5 * (mhr.beta - mhr0.beta) * (1 / ctr->sigma2)) -
        (log(4 * ctr->sigma2 * treevar) * 0.5 * (mhr.totTerm - mhr0.totTerm)) +
        log(mhr.cdf) - log(mhr0.cdf);
        
      if (ctr->debug)
        Rcout << "\n\tterm0 = " << mhr0.totTerm << " term1 = " << mhr.totTerm << "\n\tratio = " << ratio;
      
      if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
        mhr0 = mhr;
        success = 2;
        nestedTree->accept();
        if (!(ctr->binomial))
          tree->nodevals->tempV = mhr0.tempV;
      } // end MHR accept/reject
    }
    nestedTree->reject();
  } // end loop to update nested trees
  
  // Update variance and residuals
  if (ctr->shrinkage)
    rHalfCauchyFC(&(ctr->tau(t)), mhr0.totTerm, 
                mhr0.termT2 / (ctr->sigma2 * ctr->nu));
  
  if ((ctr->tau)(t) != (ctr->tau)(t)) 
    stop("\nNaN values occured during model run, rerun model.\n");
  
  if (ctr->debug)
    Rcout << "\n\t tau = " << ctr->tau(t);
  ctr->nTerm(t) =     mhr0.totTerm;
  ctr->totTerm +=     mhr0.totTerm;
  ctr->sumTermT2 +=   mhr0.termT2 / ctr->tau(t);
  ctr->Rmat.col(t) =  mhr0.Xd * mhr0.draw;
  mhr0.draw =         mhr0.Dtrans * mhr0.draw;
  
  if (ctr->debug)
    Rcout << "\nMax draw = " << mhr0.draw.maxCoeff() << "\n";
  
  // Record
  if (ctr->record > 0) {
    Eigen::VectorXd rec(8);
    rec << ctr->record, t, 0, 0, 0, 0, 0, 0;
    
    int k = 0;
    for(Node* eta : dlnmTerm) {
      
      // sort nodes in increasing xmin order
      std::vector<std::pair<int, Node*> > nodeList;
      for (Node* lambda : eta->nodevals->nestedTree->listTerminal())
        nodeList.push_back(std::make_pair(lambda->nodestruct->get(1), lambda));
      sort(nodeList.begin(), nodeList.end());
      
      for (std::pair<int, Node*> lambda : nodeList) {
        rec[2] = (lambda.second->nodestruct)->get(1);
        rec[3] = (lambda.second->nodestruct)->get(2);
        rec[4] = (lambda.second->nodestruct)->get(3);
        rec[5] = (lambda.second->nodestruct)->get(4);
        if (lambda.first > 0) {
          rec[6] = mhr0.draw(k);
          ++k;
        } else {
          rec[6] = 0.0;
        }
        (dgn->DLMexp).push_back(rec);
      } // end loop over nested tree term
    } // end loop over dlnmTerm
  } // end record
} // end function monoTDLNMTreeUpdate


/**
* @brief monotone tdlnm
* 
* @param model list of data and model control specs from R
* @return Rcpp::List 
*/
// [[Rcpp::export]]
Rcpp::List monotdlnm_Cpp(const Rcpp::List model)
{
   // * Set up model control
  tdlmCtr *ctr =      new tdlmCtr;
  ctr->debug =        as<bool>(model["debug"]);
  if (ctr->debug)
    Rcout << "Create data\n";
  ctr->iter =         as<int>(model["nIter"]);
  ctr->burn =         as<int>(model["nBurn"]);
  ctr->thin =         as<int>(model["nThin"]);
  ctr->nRec =         floor(ctr->iter / ctr->thin);
  ctr->nTrees =       as<int>(model["nTrees"]);
  ctr->stepProb =     as<std::vector<double> >(model["stepProb"]);
  ctr->treePrior =    as<std::vector<double> >(model["treePriorTime"]);
  ctr->treePrior2 =   as<std::vector<double> >(model["treePriorExp"]);
  ctr->zirtP0 =       as<VectorXd>(model["zirtp0"]);
  ctr->zirtP0 = (ctr->zirtP0.array() / (1.0 - ctr->zirtP0.array())).log();
  ctr->zirtAlpha =       as<double>(model["zirtAlpha"]);
  ctr->binomial =     as<bool>(model["binomial"]);
  ctr->shrinkage =    as<int>(model["shrinkage"]);
  ctr->verbose =      as<bool>(model["verbose"]);
  ctr->diagnostics =  as<bool>(model["diagnostics"]);
  ctr->modKappa = 1.0;
  
  // * Set up model data
  ctr->Y =            as<VectorXd>(model["Y"]);
  ctr->n =            ctr->Y.size();
  ctr->Z =            as<MatrixXd>(model["Z"]);
  ctr->pZ =           ctr->Z.cols();
  ctr->Zw =           ctr->Z;
  MatrixXd VgInv =    ctr->Z.transpose() * ctr->Z;
  VgInv.diagonal().array() += 1.0 / 1000.0;
  ctr->Vg =           VgInv.inverse();
  ctr->VgChol =       ctr->Vg.llt().matrixL();
  VgInv.resize(0,0);
  
  // * Set up data for logistic model
  ctr->binomialSize.resize(ctr->n);            ctr->binomialSize.setZero();
  ctr->kappa.resize(ctr->n);                   ctr->kappa.setOnes();
  ctr->Omega.resize(ctr->n);                   ctr->Omega.setOnes();
  if (ctr->binomial) {
    ctr->binomialSize =   as<VectorXd>(model["binomialSize"]);
    ctr->kappa =          ctr->Y - 0.5 * (ctr->binomialSize);
    ctr->Y =              ctr->kappa;
  }
  
  // * Create exposure data management
  if (ctr->debug)
    Rcout << "Create expsoureDat\n";
  exposureDat *Exp;
  if (ctr->binomial)
    Exp = new exposureDat(as<MatrixXd>(model["X"]), as<MatrixXd>(model["SE"]),
                          as<VectorXd>(model["Xsplits"]), 
                          as<MatrixXd>(model["Xcalc"]),
                          as<MatrixXd>(model["Tcalc"]),
                          as<bool>(model["lowmem"]));
  else
    Exp = new exposureDat(as<MatrixXd>(model["X"]), as<MatrixXd>(model["SE"]),
                          as<VectorXd>(model["Xsplits"]),
                          as<MatrixXd>(model["Xcalc"]),
                          as<MatrixXd>(model["Tcalc"]), ctr->Z, ctr->Vg,
                          as<bool>(model["lowmem"]));
  ctr->pX =         Exp->pX;
  ctr->nSplits =    Exp->nSplits;
  ctr->timeCounts.resize(ctr->pX); ctr->timeCounts.setZero();
  ctr->zirtPsi0 =   ctr->zirtP0;

  std::vector<MatrixXd> zirtCovInv;
  std::vector<double> zirtCovDet;
  for (int i = 1; i < 20; ++i) {
    MatrixXd covMat(ctr->pX, ctr->pX); covMat.setZero();
    for (int j = 1; j < ctr->pX; ++j) {
      covMat.diagonal(j).array() = j;
      covMat.diagonal(-j).array() = j;
    }
    // covMat = covMat.array().square();
    covMat *= log(i * 0.05);
    covMat = covMat.array().exp();
    covMat *= 1 / ctr->zirtAlpha;
    MatrixXd covDet = covMat.llt().matrixL();
    zirtCovInv.push_back(covMat.inverse());
    zirtCovDet.push_back(covDet.diagonal().array().log().sum());
  }
  int curCov = 9;
  ctr->zirtCov = zirtCovInv[9];
  
  // * Calculations used in special case: single-node trees
  ctr->X1 =         Exp->Tcalc.col(ctr->pX - 1);
  ctr->ZtX1 =       ctr->Z.transpose() * ctr->X1;
  ctr->VgZtX1 =     ctr->Vg.selfadjointView<Lower>() * ctr->ZtX1;
  ctr->VTheta1Inv = ctr->X1.dot(ctr->X1) - ctr->ZtX1.dot(ctr->VgZtX1);

  // * Create trees
  if (ctr->debug)
    Rcout << "Create nodeStruct\n";
  int t;
  std::vector<Node*> trees;
  NodeStruct *nsT;
  NodeStruct *nsX;
  bool updateTimeProb = as<bool>(model["updateTimeProb"]);
  VectorXd timeProbs0 = as<VectorXd>(model["timeProb"]);
  nsT = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                      0.0 * as<VectorXd>(model["splitProb"]), 
                      timeProbs0);
  nsX = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                      as<VectorXd>(model["splitProb"]), 
                      0.0 * timeProbs0);
  
  if (ctr->debug)
    Rcout << "Create trees\n";
  for (t = 0; t < ctr->nTrees; t++) {
    trees.push_back(new Node(0, 1));
    trees[t]->nodestruct = nsT->clone();
    trees[t]->nodevals = new NodeVals(ctr->n, ctr->pZ);
    drawZirt(trees[t], ctr, nsX);
    for (Node* lambda : trees[t]->nodevals->nestedTree->listTerminal())
      Exp->updateNodeVals(lambda); 
  }
  delete nsT;
  ctr->nTerm.resize(ctr->nTrees);                   ctr->nTerm.setOnes();
  ctr->Rmat.resize(ctr->n, ctr->nTrees);            ctr->Rmat.setZero();
  
  // * Setup model logs
  tdlmLog *dgn = new tdlmLog;
  dgn->gamma.resize(ctr->pZ, ctr->nRec);            dgn->gamma.setZero();
  dgn->sigma2.resize(ctr->nRec);                    dgn->sigma2.setZero();
  dgn->nu.resize(ctr->nRec);                        dgn->nu.setZero();
  dgn->tau.resize(ctr->nTrees, ctr->nRec);          dgn->tau.setZero();
  dgn->fhat.resize(ctr->n);                         dgn->fhat.setZero();
  dgn->fhat2.resize(ctr->n);                        dgn->fhat2.setZero();
  dgn->termNodes.resize(ctr->nTrees, ctr->nRec);    dgn->termNodes.setZero();
  dgn->zirtPsi0.resize(ctr->pX, ctr->nRec);         dgn->zirtPsi0.setZero();
  dgn->zirtCov.resize(ctr->nRec);                   dgn->zirtCov.setZero();
  dgn->kappa.resize(ctr->nRec);                   dgn->kappa.setZero();
  dgn->timeProbs.resize(ctr->pX - 1, ctr->nRec);   dgn->timeProbs.setZero();
  dgn->timeCounts.resize(ctr->pX, ctr->nRec);   dgn->timeCounts.setZero();
  
  // * Initial values and draws
  if (ctr->debug)
    Rcout << "Initial draws\n";
  ctr->fhat.resize(ctr->n);                         ctr->fhat.setZero();
  ctr->gamma.resize(ctr->pZ);                       ctr->gamma.setZero();
  ctr->tau.resize(ctr->nTrees);                     ctr->tau.setOnes();
  ctr->R =          ctr->Y;
  ctr->totTerm =    0;
  ctr->sumTermT2 =  0.0;
  ctr->nu =         1.0; // ! Need to define nu and sigma2 prior to ModelEst
  ctr->sigma2 =     1.0;
  tdlmModelEst(ctr);     // initial draws for gamma, sigma2, omega (binomial)
  rHalfCauchyFC(&(ctr->nu), ctr->nTrees, 0.0);
  if (ctr->shrinkage) {
    for (t = 0; t < ctr->nTrees; t++) 
      rHalfCauchyFC(&(ctr->tau(t)), 0.0, 0.0);
  }
  
  // * Create progress meter
  progressMeter* prog = new progressMeter(ctr);
  
  std::size_t s;
  // * Begin MCMC run
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    Rcpp::checkUserInterrupt();
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }
    
    // * Update trees
    ctr->R +=         ctr->Rmat.col(0);
    ctr->totTerm =    0.0; 
    ctr->sumTermT2 =  0.0;
    ctr->fhat.setZero();
    for (t = 0; t < ctr->nTrees; t++) {
      if (ctr->debug)
        Rcout << "\n" << t << ":";
      monoTDLNMTreeUpdate(t, trees[t], ctr, dgn, Exp, nsX);
      ctr->fhat += ctr->Rmat.col(t);
      if (t < ctr->nTrees - 1)
        ctr->R += ctr->Rmat.col(t + 1) - ctr->Rmat.col(t);
    } // end update trees

    // * Update model
    ctr->R = ctr->Y - ctr->fhat;
    tdlmModelEst(ctr);
    rHalfCauchyFC(&(ctr->nu), ctr->totTerm, ctr->sumTermT2 / ctr->sigma2);
    if (ctr->debug)
      Rcout << "\nVar = " << ctr->sigma2 << " " << ctr->nu;
    if ((ctr->sigma2 != ctr->sigma2) || (ctr->nu != ctr->nu))
      stop("\nNaN values (sigma2, nu) occured during model run, rerun model.");
      




    // Update zero-inflated tree split probabilities
    VectorXd timeSplits(ctr->pX - 1); timeSplits.setZero();
    ctr->timeCounts.setZero();
    if ((ctr->b > 1000) || (ctr->b > (0.5 * ctr->burn))) {

      // Count time-tree terminal nodes for logistic model estimation
      int nTerm = 0;
      for (Node* tree : trees) {
        timeSplits += countTimeSplits(tree, ctr);
        nTerm += tree->nTerminal();
      }

      // Create split outcome (cwY) and design matrix (cwX)
      VectorXd cwY(nTerm); cwY.array() = -0.5;// polya-gamma: y=y-1/2
      MatrixXd cwX(nTerm, ctr->pX); cwX.setZero();
      Eigen::Index termIt = 0;

      for (Node* tree : trees) { // loop over all trees
        for (Node* eta : tree->listTerminal(0)) { 
          int tmin = eta->nodestruct->get(3);
          int tmax = eta->nodestruct->get(4);  
          cwX.row(termIt).segment(tmin - 1, tmax - tmin + 1).array() = 1.0;

          if (eta->nodevals->nestedTree->c1 != 0) {// single node tree
            cwY(termIt) += 1.0;
            ctr->timeCounts.segment(tmin - 1, tmax - tmin + 1).array() = 1.0;
          }
          // Rcout << "\n" << cwY(termIt) << " " << tmin << " " << tmax << " " << cwX.row(termIt);

          ++termIt;
        }
      }
      
      // draw polya-gamma for CW var selection, repeat 10x to help convergence
      VectorXd cwOnes(nTerm); cwOnes.setOnes();
      VectorXd zirtProb = ctr->zirtPsi0;
      MatrixXd cwV = cwX.transpose() * cwX;
      for (int i = 0; i < 10; ++i) {
        VectorXd psi = cwX * zirtProb;
        VectorXd cwPG = rcpp_pgdraw(cwOnes, psi);
        cwV = cwX.transpose() * cwPG.asDiagonal() * cwX;
        cwV += ctr->zirtCov;
        zirtProb = cwV.inverse() * (cwX.transpose() * cwY + ctr->zirtCov * ctr->zirtP0);
        zirtProb += cwV.inverse().llt().matrixL() * as<VectorXd>(rnorm(ctr->pX, 0, 1));
      }
      ctr->zirtPsi0 = zirtProb;
      

      // update CW var selection covariance matrix
      double covMHR = 0.0;
      int newCov = curCov;
      if (curCov == 0) {
        newCov = 1;
        covMHR += log(0.5);
      } else if (curCov == 18) {
        newCov = 17;
        covMHR += log(0.5);
      } else {
        if (R::runif(0.0, 1.0) < 0.5) {
          ++newCov;
        } else {
          --newCov;
        }
      }
      covMHR += -zirtCovDet[newCov] - 
        (ctr->zirtPsi0 - ctr->zirtP0).dot(zirtCovInv[newCov] * (ctr->zirtPsi0 - ctr->zirtP0)) + 
        zirtCovDet[curCov] + (ctr->zirtPsi0 - ctr->zirtP0).dot(zirtCovInv[curCov] * (ctr->zirtPsi0 - ctr->zirtP0));
      if (log(R::runif(0.0, 1.0)) < covMHR) {
        curCov = newCov;
        ctr->zirtCov = zirtCovInv[curCov];
      }

      // if (updateTimeProb) { // update time splitting probabilities
        VectorXd timeProbs = trees[0]->nodestruct->getTimeProbs();
        double beta = R::rbeta(1.0, 1.0);
        double modKappaNew = beta * (ctr->pX - 1.0)/ (1 - beta);
        double mhrDir = logDirichletDensity(timeProbs, timeSplits + modKappaNew * timeProbs0) - 
          logDirichletDensity(timeProbs, timeSplits + ctr->modKappa * timeProbs0);
        if (log(R::runif(0, 1)) < mhrDir)
          ctr->modKappa = modKappaNew;

        VectorXd newTimeProbs = rDirichlet(timeSplits + ctr->modKappa * timeProbs0);


        // update tree time split probabilities
        for (Node* tree : trees) {
          tree->nodestruct->setTimeProbs(newTimeProbs);
          tree->updateStruct();
        }
      // }
    } // end update of split and zero-inflated probabilities
    
    if (ctr->debug)
      Rcout << "\nsigma2=" << ctr->sigma2 << " nu = " << ctr->nu;
      
    // * Record
    if (ctr->record > 0) {
      dgn->gamma.col(ctr->record - 1) =     ctr->gamma;
      dgn->sigma2(ctr->record - 1) =        ctr->sigma2;
      dgn->nu(ctr->record - 1) =            ctr->nu;
      dgn->tau.col(ctr->record - 1) =       ctr->tau;
      dgn->termNodes.col(ctr->record - 1) = ctr->nTerm;
      dgn->fhat +=                          ctr->fhat;
      dgn->fhat2 +=                         ctr->fhat.array().square().matrix();
      dgn->zirtPsi0.col(ctr->record - 1) =  ctr->zirtPsi0;
      dgn->zirtCov(ctr->record - 1) =            curCov;
      dgn->kappa(ctr->record - 1) =            ctr->modKappa;
      dgn->timeProbs.col(ctr->record -1) = trees[0]->nodestruct->getTimeProbs();
      dgn->timeCounts.col(ctr->record - 1) = ctr->timeCounts;
    }
    
    // * Update progress
    prog->printMark();
  } // end MCMC
  
  // * Setup data for return
  Eigen::MatrixXd DLM((dgn->DLMexp).size(), 8);
  for (s = 0; s < (dgn->DLMexp).size(); ++s)
    DLM.row(s) = dgn->DLMexp[s];
  Eigen::VectorXd sigma2 = dgn->sigma2;
  Eigen::VectorXd nu = dgn->nu;
  Eigen::VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  Eigen::MatrixXd gamma = (dgn->gamma).transpose();
  Eigen::MatrixXd zirtPsi0 = dgn->zirtPsi0.transpose();
  Eigen::VectorXd kappa = dgn->kappa;
  Eigen::VectorXd zirtCov = dgn->zirtCov;
  Eigen::MatrixXd tau = (dgn->tau).transpose();
  Eigen::MatrixXd termNodes = (dgn->termNodes).transpose();
  Eigen::MatrixXd timeProbs = (dgn->timeProbs).transpose();
  Eigen::MatrixXd timeCounts = (dgn->timeCounts).transpose();
  Eigen::MatrixXd Accept((dgn->TreeAccept).size(), 5);
  for (s = 0; s < (dgn->TreeAccept).size(); ++s)
    Accept.row(s) = dgn->TreeAccept[s];
  delete prog;
  delete ctr;
  delete dgn;
  delete Exp;
  for (s = 0; s < trees.size(); ++s)
    delete trees[s];

  return(Rcpp::List::create(
    Named("DLM")    = wrap(DLM),
    Named("fhat")   = wrap(fhat),
    Named("sigma2") = wrap(sigma2),
    Named("nu")     = wrap(nu),
    Named("kappa")     = wrap(kappa),
    Named("zirtCov")     = wrap(zirtCov),
    Named("tau")    = wrap(tau),
    Named("termNodes")  = wrap(termNodes),
    Named("gamma")  = wrap(gamma),
    Named("zirt") = wrap(zirtPsi0),
    Named("timeProbs") = wrap(timeProbs),
    Named("timeCounts") = wrap(timeCounts),
    Named("treeAccept") = wrap(Accept)));
} // end function monotdlnm_Cppa