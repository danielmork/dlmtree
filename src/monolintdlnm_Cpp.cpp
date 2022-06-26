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
 
#include <RcppEigen.h>
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


class exposureDatPL {
public:
  int n, nSplits, pX;
  int pZ = 0;
  double Xmin, Xmax;
  VectorXd Xsplits;
  std::vector<MatrixXd> Xcalc, X1, ZtXcalc, ZtX1, VgZtXcalc, VgZtX1; 
  // sum of X < Xmax and count matrices
  // Xcalc[Xmax](index, tmax) = sum_t_from_0_to_tmax X(i,t) if X(i,t) < Xmax
  // X1[Xmax](index, tmax) = sum_t_from_0_to_tmax 1 if X(i,t) < Xmax
  
  /**
   * @brief Construct a new exposure Dat P L object
   * 
   * @param X 
   * @param Tcalc_in 
   * @param Xsplits_in 
   */
  exposureDatPL(MatrixXd X, MatrixXd Tcalc_in, VectorXd Xsplits_in, 
                MatrixXd Z, MatrixXd Vg)
  {
    n = X.rows();
    nSplits = Xsplits_in.size();
    pX = X.cols();
    pZ = Z.cols();
    Xsplits = Xsplits_in;
    Xmin = X.minCoeff();
    Xmax = X.maxCoeff();
    
    // subset data based on X < Xsplits
    for (int s = 0; s < nSplits; ++s) { // populate X < Xsplits
      Xcalc.push_back(MatrixXd::Zero(n, pX));
      X1.push_back(MatrixXd::Zero(n, pX));
      
      for (int tmax = 1; tmax <= pX; ++tmax) {
        for (int t = 0; t < tmax; ++t) {
          for (int i = 0; i < n; ++i) {
            
            if (X(i, t) < Xsplits[s]) {
              Xcalc.back()(i, tmax - 1) += X(i, t);
              X1.back()(i, tmax - 1) += 1.0;
            } 
            
          } // loop across indices
        } // sum_t from 0 to tmax X(i,t) if X(i,t) < Xmax
      } // loop over tmax
      ZtXcalc.push_back(Z.transpose() * Xcalc.back());
      ZtX1.push_back(Z.transpose() * X1.back());
      VgZtXcalc.push_back(Vg * ZtXcalc.back());
      VgZtX1.push_back(Vg * ZtX1.back());
    } // loop over Xmin
    
    // X < Inf (all data subset)
    Xcalc.push_back(Tcalc_in);
    X1.push_back(MatrixXd::Zero(n, pX));
    for (int tmax = 0; tmax < pX; ++tmax) { // Populate count matrix
      X1.back().col(tmax).array() += tmax + 1.0;
    }
    ZtXcalc.push_back(Z.transpose() * Xcalc.back());
    VgZtXcalc.push_back(Vg * ZtXcalc.back());
    ZtX1.push_back(Z.transpose() * X1.back());
    VgZtX1.push_back(Vg * ZtX1.back());
  } // end function exposureDatPL
  
  /**
   * @brief Calculate sum_t from tmin to tmax: X(i,t) if xmin <= X(i,t) < xmax
   * 
   * @param tmin int node min time (starts at 1)
   * @param tmax int node max time
   * @param xmin int node min xsplit (0 = negInf)
   * @param xmax int node max xsplit (nSplits+1 = posInf)
   * @return VectorXd 
   */
  VectorXd returnXvec(int tmin, int tmax, int xmin, int xmax)
  {
    VectorXd Xout = Xcalc[xmax-1].col(tmax-1);
    if (xmin > 0) // [xmax(tmax)-xmin(tmax)] ! no tmin
      Xout -= Xcalc[xmin-1].col(tmax-1);
    if (tmin > 1) { // [xmax(tmax)-xmax(tmin)] ! no xmin
      Xout -= Xcalc[xmax-1].col(tmin-2);
      if (xmin > 0) // [xmax(tmax)-xmax(tmin)]-[xmin(tmax)-xmin(tmin)]
        Xout += Xcalc[xmin-1].col(tmin-2);
    }
    return (Xout);
  } // end returnXvec function
  VectorXd returnZtXvec(int tmin, int tmax, int xmin, int xmax)
  {
    VectorXd ZtXout = ZtXcalc[xmax-1].col(tmax-1);
    if (xmin > 0) // [xmax(tmax)-xmin(tmax)] ! no tmin
      ZtXout -= ZtXcalc[xmin-1].col(tmax-1);
    if (tmin > 1) { // [xmax(tmax)-xmax(tmin)] ! no xmin
      ZtXout -= ZtXcalc[xmax-1].col(tmin-2);
      if (xmin > 0) // [xmax(tmax)-xmax(tmin)]-[xmin(tmax)-xmin(tmin)]
        ZtXout += ZtXcalc[xmin-1].col(tmin-2);
    }
    return (ZtXout);
  } // end returnZtXvec function
  VectorXd returnVgZtXvec(int tmin, int tmax, int xmin, int xmax)
  {
    VectorXd VgZtXout = VgZtXcalc[xmax-1].col(tmax-1);
    if (xmin > 0) // [xmax(tmax)-xmin(tmax)] ! no tmin
      VgZtXout -= VgZtXcalc[xmin-1].col(tmax-1);
    if (tmin > 1) { // [xmax(tmax)-xmax(tmin)] ! no xmin
      VgZtXout -= VgZtXcalc[xmax-1].col(tmin-2);
      if (xmin > 0) // [xmax(tmax)-xmax(tmin)]-[xmin(tmax)-xmin(tmin)]
        VgZtXout += VgZtXcalc[xmin-1].col(tmin-2);
    }
    return (VgZtXout);
  } // end returnZtXvec function
  /**
   * @brief Similar to returnXvec, but returns counts
   * @return VectorXd 
   */
  VectorXd returnX1vec(int tmin, int tmax, int xmin, int xmax)
  {
    VectorXd Xout = X1[xmax-1].col(tmax-1);
    if (xmin > 0)
      Xout -= X1[xmin-1].col(tmax-1);
    if (tmin > 1) {
      Xout -= X1[xmax-1].col(tmin-2);
      if (xmin > 0)
        Xout += X1[xmin-1].col(tmin-2);
    }
    return (Xout);
  } // end returnX1 function
  VectorXd returnZtX1vec(int tmin, int tmax, int xmin, int xmax)
  {
    VectorXd ZtXout = ZtX1[xmax-1].col(tmax-1);
    if (xmin > 0)
      ZtXout -= ZtX1[xmin-1].col(tmax-1);
    if (tmin > 1) {
      ZtXout -= ZtX1[xmax-1].col(tmin-2);
      if (xmin > 0)
        ZtXout += ZtX1[xmin-1].col(tmin-2);
    }
    return (ZtXout);
  } // end returnZtX1 function
  VectorXd returnVgZtX1vec(int tmin, int tmax, int xmin, int xmax)
  {
    VectorXd VgZtXout = VgZtX1[xmax-1].col(tmax-1);
    if (xmin > 0)
      VgZtXout -= VgZtX1[xmin-1].col(tmax-1);
    if (tmin > 1) {
      VgZtXout -= VgZtX1[xmax-1].col(tmin-2);
      if (xmin > 0)
        VgZtXout += VgZtX1[xmin-1].col(tmin-2);
    }
    return (VgZtXout);
  } // end returnVgZtX1 function
  
  /**
   * @brief return X value given splitting index (0 = min, nSplits+1 = max)
   * 
   * @param xsplitIdx splitting index
   * @return double 
   */
  double Xval(int xsplitIdx) {
    if (xsplitIdx == 0)
      return (Xmin);
    if (xsplitIdx == nSplits + 1)
      return (Xmax);
    return (Xsplits[xsplitIdx - 1]);
  }
  
  /**
   * @brief create matrix of X values for piecewise linear weight function
   * 
   * @param nestedTree 
   * @param ntTerm 
   * @param tmin 
   * @param tmax 
   * @return MatrixXd 
   */
  void nodeVals(Node* eta, std::vector<Node*> ntTerm, 
                bool proposed, bool binary = 0) {
    if (ntTerm.size() == 1) // no effect
      return;
    
    // sort node limits smallest to largest
    std::vector<std::pair<int, int> > xminmax;
    for (Node* lambda : ntTerm)
      xminmax.push_back(std::make_pair(lambda->nodestruct->get(1), 
                                       lambda->nodestruct->get(2)));
    std::sort(xminmax.begin(), xminmax.end());
    int tmin = eta->nodestruct->get(3);
    int tmax = eta->nodestruct->get(4);
    
    
    // create matrix of calculated piecewise linear X values               
    MatrixXd Xmat = MatrixXd::Zero(n, ntTerm.size() - 1);             
    MatrixXd ZtXmat = MatrixXd::Zero(pZ, ntTerm.size() - 1);             
    MatrixXd VgZtXmat = MatrixXd::Zero(pZ, ntTerm.size() - 1);
    // Rcout << ":" << xminmax[0].first << "-" << xminmax[0].second << ":";
    Xmat.col(0) = // populate first column of matrix
      (returnXvec(tmin, tmax, xminmax[1].first, xminmax[1].second) - 
       returnX1vec(tmin, tmax, xminmax[1].first, xminmax[1].second) * 
        Xval(xminmax[1].first)) /
      (Xval(xminmax[1].second) - Xval(xminmax[1].first));
    if (!binary) {
      ZtXmat.col(0) = // populate first column of matrix
        (returnZtXvec(tmin, tmax, xminmax[1].first, xminmax[1].second) - 
        returnZtX1vec(tmin, tmax, xminmax[1].first, xminmax[1].second) * 
          Xval(xminmax[1].first)) /
        (Xval(xminmax[1].second) - Xval(xminmax[1].first));
      VgZtXmat.col(0) = // populate first column of matrix
        (returnVgZtXvec(tmin, tmax, xminmax[1].first, xminmax[1].second) - 
        returnVgZtX1vec(tmin, tmax, xminmax[1].first, xminmax[1].second) * 
          Xval(xminmax[1].first)) /
        (Xval(xminmax[1].second) - Xval(xminmax[1].first));
    }
    
    
    // for each remaining node, add to previous and current column of matrix
    if (ntTerm.size() > 2) {
      for (size_t i = 2; i < ntTerm.size(); ++i) {
        // Rcout << ":" << xminmax[i].first << "-" << xminmax[i].second << ":";
        double lower = Xval(xminmax[i].first);
        double upper = Xval(xminmax[i].second);
        VectorXd tempX1 = 
          returnX1vec(tmin, tmax, xminmax[i].first, xminmax[i].second);
        VectorXd tempXvec =
          (returnXvec(tmin, tmax, xminmax[i].first, xminmax[i].second) -
          tempX1 * lower) / (upper - lower);
        Xmat.col(i-2) += tempX1 - tempXvec;
        Xmat.col(i-1) += tempXvec;
        
        if (!binary) {
          VectorXd tempZtX1 = 
            returnZtX1vec(tmin, tmax, xminmax[i].first, xminmax[i].second);
          VectorXd tempVgZtX1 = 
            returnVgZtX1vec(tmin, tmax, xminmax[i].first, xminmax[i].second);
          VectorXd tempZtXvec =
            (returnZtXvec(tmin, tmax, xminmax[i].first, xminmax[i].second) -
            tempZtX1 * lower) / (upper - lower);
          VectorXd tempVgZtXvec =
            (returnVgZtXvec(tmin, tmax, xminmax[i].first, xminmax[i].second) -
            tempVgZtX1 * lower) / (upper - lower);
            
          ZtXmat.col(i-2) += tempZtX1 - tempZtXvec;
          VgZtXmat.col(i-2) += tempVgZtX1 - tempVgZtXvec;
          ZtXmat.col(i-1) += tempZtXvec;
          VgZtXmat.col(i-1) += tempVgZtXvec;
        }
      }
    } // end loop over remaining nodes
    
    
    // update nodevals  
    if (proposed) {
      eta->nodevals->XplProposed = Xmat;
      eta->nodevals->ZtXmatProposed = ZtXmat;
      eta->nodevals->VgZtXmatProposed = VgZtXmat;
    } else {
      eta->nodevals->Xpl = Xmat;
      eta->nodevals->ZtXmat = ZtXmat;
      eta->nodevals->VgZtXmat = VgZtXmat;
    }
    return;
  } // end nodeVals function
}; // end exposureDatPL class


treeMHR monoDlnmMHR(std::vector<Node*> dlnmTerm, tdlmCtr* ctr,  
                     exposureDatPL* Exp, VectorXd ZtR,
                     double treevar, Node* tree, bool updateNested)
{
  if (ctr->debug)
    Rcout << "\nMHR";
  treeMHR out;
  int totTerm = 0;
  std::vector<std::vector<Node* > > nestedTerm;
  for (Node* eta : dlnmTerm) {
    nestedTerm.push_back(eta->nodevals->nestedTree->listTerminal(updateNested));
    totTerm += nestedTerm.back().size() - 1;
      
    if ((nestedTerm.back().size() > 1) && eta->update)
      Exp->nodeVals(eta, nestedTerm.back(), updateNested);
    
    if (ctr->debug)
      Rcout << " nTerm = " << nestedTerm.back().size();
  }
    
  if (totTerm == 0) { // no terminal node effects
    if (ctr->debug)
      Rcout << "empty";
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
  
  // Define transformation matrices (D and Dinv), fill Xd, ZtX, VgZtX
  int k = 0;  // index for D, Dinv
  int l;      // size of block D matrix
  int tmin, tmax;
  NodeVals* nv;
  for (size_t t = 0; t < nestedTerm.size(); ++t) {
    if (nestedTerm[t].size() == 1) continue;
    
    l = nestedTerm[t].size() - 1;
    Dinv.block(k, k, l, l).triangularView<Lower>().setOnes();
    
    tmin = dlnmTerm[t]->nodestruct->get(3);
    tmax = dlnmTerm[t]->nodestruct->get(4);
    nv =   dlnmTerm[t]->nodevals;
    
    if (updateNested && dlnmTerm[t]->update) {
      if (ctr->debug)
        Rcout << "rc1." << tmin << "," << tmax << "," << ctr->n << "=" << nv->XplProposed.rows() << "," << l << "=" << nv->XplProposed.cols() <<
        " mean=" << nv->XplProposed.mean();
      out.Xd.block(0, k, ctr->n, l) = nv->XplProposed;
      
      if (ctr->binomial) {
        ZtX.block(0, k, ctr->pZ, l) = ctr->Zw.transpose() * nv->XplProposed;
        VgZtX.block(0, k, ctr->pZ, l) = ctr->Vg * ZtX.block(0, k, ctr->pZ, l);
      } else {
        ZtX.block(0, k, ctr->pZ, l) = nv->ZtXmatProposed;
        VgZtX.block(0, k, ctr->pZ, l) = nv->VgZtXmatProposed;
      }
      
    } else {
      if (ctr->debug)
        Rcout << " rc2." << tmin << "," << tmax << "," << ctr->n << "=" << nv->Xpl.rows() << "," << l << "=" << nv->Xpl.cols() <<
        " mean=" << nv->Xpl.mean();
      out.Xd.block(0, k, ctr->n, l) = nv->Xpl;
      
      if (ctr->binomial) {
        ZtX.block(0, k, ctr->pZ, l) = ctr->Zw.transpose() * nv->Xpl;
        VgZtX.block(0, k, ctr->pZ, l) = ctr->Vg * ZtX.block(0, k, ctr->pZ, l);
      } else {
        ZtX.block(0, k, ctr->pZ, l) = nv->ZtXmat;
        VgZtX.block(0, k, ctr->pZ, l) = nv->VgZtXmat;
      }
    }
    
    k += l;
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
  out.draw =          rtmvnorm(ThetaHat, ctr->sigma2 * VTheta, 2);    
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
 * @brief update tree
 * 
 * @param t 
 * @param tree 
 * @param ctr 
 * @param dgn 
 * @param Exp 
 * @param nsX 
 */
void monoTDLNMTreeUpdate(int t, Node* tree, tdlmCtr* ctr, tdlmLog* dgn, 
                       exposureDatPL* Exp, NodeStruct* nsX)
{
  int step =        0;
  int success =     0;
  double stepMhr =  0.0;
  double ratio =    0.0;
  double treevar =  ctr->nu * ctr->tau(t);
  std::vector<Node*> dlnmTerm = tree->listTerminal();
  std::vector<Node*> newDlnmTerm, nestedTerm;
  Node* nestedTree;
  treeMHR mhr0, mhr;
  
  // List current tree terminal nodes
  VectorXd ZtR = ctr->Zw.transpose() * ctr->R;
  if (dlnmTerm.size() > 1) // if single terminal node, grow is only option
    step = sampleInt(ctr->stepProb, 1);
  
  // propose update
  stepMhr = tdlmProposeTree(tree, 0, ctr, step);
  success = tree->isProposed();
  if (ctr->debug)
    Rcout << "\nstep = " << step << " nnodes = " << dlnmTerm.size() << " stepMhr = " << stepMhr;
  mhr0 = monoDlnmMHR(dlnmTerm, ctr, Exp, ZtR, treevar, tree, 0);
  
  if (success && (stepMhr == stepMhr)) {
    newDlnmTerm = tree->listTerminal(1);
    
    // draw nested trees if grow/prune
    if (step < 2) {
      for (Node* eta : newDlnmTerm) {
        if (eta->nodevals == 0)
          eta->nodevals = new NodeVals(ctr->n);
        if (eta->update) {
          drawZirt(eta, ctr, nsX);
        }
      }
    
    // change proposal for time tree
    } else {
      for (Node* eta : dlnmTerm) { // calculate MHR
        int tmin = eta->nodestruct->get(3);
        int tmax = eta->nodestruct->get(4);
        if (eta->nodevals->nestedTree->c1 == 0)
          stepMhr -= logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 1);
        else
          stepMhr -= logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 0);
      }
      
      for (Node* eta : newDlnmTerm) {
        int tmin = eta->nodestruct->get(3);
        int tmax = eta->nodestruct->get(4);
        if (eta->nodevals->nestedTree->c1 == 0)
          stepMhr += logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 1);
        else
          stepMhr += logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 0);
      }
    } // end change proposal
    
    mhr = monoDlnmMHR(newDlnmTerm, ctr, Exp, ZtR, treevar, tree, 1);
    
    // combine mhr parts into log-MH ratio
    ratio = stepMhr + (mhr.logVThetaChol - mhr0.logVThetaChol) +
      (0.5 * (mhr.beta - mhr0.beta) * (1 / ctr->sigma2)) -
      (log(4 * ctr->sigma2 * treevar) * 0.5 * (mhr.totTerm - mhr0.totTerm)) +
      log(mhr.cdf) - log(mhr0.cdf);
    if (ctr->debug)
      Rcout << " ratio = " << ratio << " treevar = " << treevar;
    
    if (log(R::runif(0, 1)) < ratio) {
      mhr0 = mhr;
      success = 2;
      tree->accept();
      dlnmTerm = tree->listTerminal();
      
      if (!(ctr->binomial))
        tree->nodevals->tempV = mhr0.tempV;
      
      for (Node* eta : dlnmTerm) {
        if (eta->update) {
          eta->update = 0;
          if (eta->parent != 0) {
            delete eta->parent->nodevals->nestedTree;
            eta->parent->nodevals->nestedTree = 0;
            eta->parent->update = 1;
          }
          eta->nodevals->Xpl = eta->nodevals->XplProposed;
          if (!(ctr->binomial)) {
            eta->nodevals->ZtXmat = eta->nodevals->ZtXmatProposed;
            eta->nodevals->VgZtXmat = eta->nodevals->VgZtXmatProposed;
          }
        }
      }
    }
  }
  
  tree->reject();
  
  
  // * Propose new nested tree at each terminal node
  for (Node* eta : dlnmTerm) {
    eta->update = 1;
    nestedTree = eta->nodevals->nestedTree;
    nestedTerm = nestedTree->listTerminal();
    int tmin = eta->nodestruct->get(3);
    int tmax = eta->nodestruct->get(4);
    
    step = 0;
    if (nestedTerm.size() > 1)
      step = sampleInt(ctr->stepProb, 1);
    
    if (ctr->debug)
      Rcout << "\n\t nt step = " << step;
      
    stepMhr = tdlmProposeTree(nestedTree, 0, ctr, step, 1.0);
    success = nestedTree->isProposed();
    
    if (success && (stepMhr == stepMhr)) {
      mhr = monoDlnmMHR(dlnmTerm, ctr, Exp, ZtR, treevar, tree, 1);
      
      if ((step == 0) && (nestedTerm.size() == 1)) {
        stepMhr = logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 0) -
          logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 1) + 
          2 * logPSplit((ctr->treePrior2)[0], (ctr->treePrior2)[1], 2.0, 1);
                   
      } else if ((step == 1) && (nestedTerm.size() == 2)) {
        stepMhr = logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 1) -
          logZIPSplit(ctr->zirtPsi0, tmin, tmax, ctr->nTrees, 0) -
          2 * logPSplit((ctr->treePrior2)[0], (ctr->treePrior2)[1], 2.0, 1);
      }
      
      if (ctr->debug)
        Rcout << "\n\tterm0 = " << mhr0.totTerm << " term1 = " << mhr.totTerm;
      
      ratio = stepMhr + (mhr.logVThetaChol - mhr0.logVThetaChol) +
        (0.5 * (mhr.beta - mhr0.beta) * (1 / ctr->sigma2)) -
        (log(4 * ctr->sigma2 * treevar) * 0.5 * (mhr.totTerm - mhr0.totTerm)) +
        log(mhr.cdf) - log(mhr0.cdf);
        
      if (ctr->debug)
        Rcout << "\n\tratio = " << ratio;
      
      if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
        mhr0 = mhr;
        success = 2;
        nestedTree->accept();
        
        if (!(ctr->binomial))
          tree->nodevals->tempV = mhr0.tempV;
        
        eta->nodevals->Xpl = eta->nodevals->XplProposed;
        if (!(ctr->binomial)) {
          eta->nodevals->ZtXmat = eta->nodevals->ZtXmatProposed;
          eta->nodevals->VgZtXmat = eta->nodevals->VgZtXmatProposed;
        }          
      } // end MHR accept/reject
    }
    
    nestedTree->reject();
    eta->update = 0;
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
      if (eta->nodevals->nestedTree->c1 != 0) {
        rec[4] = eta->nodestruct->get(3);
        rec[5] = eta->nodestruct->get(4);
        
        // sort nodes in increasing xmin order
        std::vector<std::pair<int, int> > xminmax;
        for (Node* lambda : eta->nodevals->nestedTree->listTerminal())
          xminmax.push_back(std::make_pair(lambda->nodestruct->get(1), 
                                           lambda->nodestruct->get(2)));
        std::sort(xminmax.begin(), xminmax.end());
        
        for (std::pair<int, int> x : xminmax) {
          if (x.first == 0) continue;
          // double denom = Exp->Xval(x.second) - Exp->Xval(x.first);
          rec[2] = x.first;
          rec[3] = x.second;
          rec[6] = mhr0.draw(k);// / denom;
          // if (x.first > 0)
          //   rec[6] -= mhr0.draw(k-1);// / denom;
          dgn->DLMexp.push_back(rec);
          ++k;
          if (ctr->debug)
            Rcout << rec.transpose() << "\n";
        } // end loop over nested tree terminal node limits
      } // end if nested tree at least 2 terminal nodes
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
 List monolintdlnm_Cpp(const List model)
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
  ctr->zirtAlpha =    as<VectorXd>(model["zirtAlpha"]);
  ctr->zirtP0 =       as<VectorXd>(model["p_t"]);
  ctr->binomial =     as<bool>(model["binomial"]);
  ctr->shrinkage =    as<int>(model["shrinkage"]);
  ctr->verbose =      as<bool>(model["verbose"]);
  ctr->diagnostics =  as<bool>(model["diagnostics"]);
  
  // * Set up model data
  ctr->Y =            as<VectorXd>(model["Y"]);
  ctr->n =            ctr->Y.size();
  ctr->Z =            as<MatrixXd>(model["Z"]);
  ctr->pZ =           ctr->Z.cols();
  ctr->Zw =           ctr->Z;
  MatrixXd VgInv =    ctr->Z.transpose() * ctr->Z;
  VgInv.diagonal().array() += 1.0 / 100000.0;
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
  exposureDatPL *Exp;
  if (ctr->debug)
    Rcout << "create Exp\n";
  Exp = new exposureDatPL(as<MatrixXd>(model["X"]), 
                          as<MatrixXd>(model["Tcalc"]),
                          as<VectorXd>(model["Xsplits"]),
                          ctr->Z, ctr->Vg);
  ctr->pX =         Exp->pX;
  ctr->nSplits =    Exp->nSplits;
  ctr->zirtPsi0 = ctr->zirtP0;

  // * Create trees
  if (ctr->debug)
    Rcout << "Create nodeStruct\n";
  int t;
  std::vector<Node*> trees;
  NodeStruct *nsT;
  NodeStruct *nsX;
  nsT = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                      0.0 * as<VectorXd>(model["splitProb"]), 
                      as<VectorXd>(model["timeProb"]));
  nsX = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                      as<VectorXd>(model["splitProb"]), 
                      0.0 * as<VectorXd>(model["timeProb"]));
  
  if (ctr->debug)
    Rcout << "Create trees\n";
  for (t = 0; t < ctr->nTrees; t++) {
    trees.push_back(new Node(0, 1));
    trees[t]->nodestruct = nsT->clone();
    trees[t]->nodevals = new NodeVals(ctr->n, ctr->pZ);
    drawZirt(trees[t], ctr, nsX);
    Exp->nodeVals(trees[t], trees[t]->nodevals->nestedTree->listTerminal(), 0);
    treeMHR mhr0 = monoDlnmMHR(trees[t]->listTerminal(0), 
                                ctr, Exp, ctr->Y, 1, trees[t], 0);
    trees[t]->nodevals->tempV = mhr0.tempV;
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
      
    // * Update zirt time split probabilities
    if ((ctr->b > 1000) || (ctr->b > (0.5 * ctr->burn))) {
      VectorXd sumSplits(ctr->pX); sumSplits.setZero();
      for (Node* tree : trees) { // loop over all trees
        for (Node* eta : tree->listTerminal(0)) { 
          int tmin = eta->nodestruct->get(3);
          int tmax = eta->nodestruct->get(4);       
          if (eta->nodevals->nestedTree->c1 != 0) { // single node tree
            sumSplits.array().segment(tmin - 1, tmax - tmin + 1) += 1.0;
          }
        }
      }
      // ctr->zirtPsi0 = rDirichlet(sumSplits);
      for (int i = 0; i < ctr->pX; ++i) {
        ctr->zirtPsi0(i) = R::rbeta(1.0 + sumSplits(i), 1.0 + ctr->nTrees - sumSplits(i));
      }

      if (ctr->debug)
        Rcout << "\nsumSplits" << sumSplits << "\nzirtPsi0" << ctr->zirtPsi0;
      // for (int i = 0; i < ctr->pX; ++i) {
      //   double newProb = R::rbeta(ctr->zirtAlpha(i), ctr->nTrees);
      //   // double newProb = R::runif(0.0, 1.0);
      //   if (log(R::runif(0, 1)) <
      //       zeroInflatedTreeMHR(ctr->zirtPsi0, trees, i, newProb) +
      //       R::dbeta(newProb, ctr->zirtAlpha(i), ctr->nTrees, 1) -
      //       R::dbeta(ctr->zirtPsi0(i), ctr->zirtAlpha(i), ctr->nTrees, 1))
      //     ctr->zirtPsi0(i) = newProb;
      // }
    }

    if (ctr->debug)
      Rcout << "\nsigma2=" << ctr->sigma2 << " nu = " << ctr->nu;
      
    // * Record
    if (ctr->record > 0) {
      // Rcout << "\nsigma2=" << ctr->sigma2 << " nu = " << ctr->nu;
      dgn->gamma.col(ctr->record - 1) =     ctr->gamma;
      dgn->sigma2(ctr->record - 1) =        ctr->sigma2;
      dgn->nu(ctr->record - 1) =            ctr->nu;
      dgn->tau.col(ctr->record - 1) =       ctr->tau;
      dgn->termNodes.col(ctr->record - 1) = ctr->nTerm;
      dgn->fhat +=                          ctr->fhat;
      dgn->fhat2 +=                         ctr->fhat.array().square().matrix();
      dgn->zirtPsi0.col(ctr->record - 1) =  ctr->zirtPsi0;
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
  Eigen::MatrixXd tau = (dgn->tau).transpose();
  Eigen::MatrixXd termNodes = (dgn->termNodes).transpose();
  Eigen::MatrixXd Accept((dgn->TreeAccept).size(), 5);
  for (s = 0; s < (dgn->TreeAccept).size(); ++s)
    Accept.row(s) = dgn->TreeAccept[s];
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
                            Named("zirt") = wrap(zirtPsi0),
                            Named("treeAccept") = wrap(Accept)));
} // end function monotdlnm_Cppa