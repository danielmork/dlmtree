#include <RcppEigen.h>
#include "exposureDat.h"
#include "Node.h"
#include "NodeStruct.h"
using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Lower;

#define DOUBLE_ERROR   0.000000000000000111022302462515654042363166809082031259999
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588

double Phi(double x1, double x2)
{
  return (erf(x2 * MATH_SQRT1_2) - erf(x1 * MATH_SQRT1_2)) * 0.5;
}


VectorXd nodeCount(//Node* n, 
               exposureDat* Exp,
               double xmin, double xmax, int tmin, int tmax)
{
  int i, j;
  VectorXd Xvec(Exp->n); Xvec.setZero();

  if (Exp->se) { // exposure measurement error or smoothing
    for (j = tmin - 1; j < tmax; ++j) {
      for (i = 0; i < Exp->n; ++i) {
        Xvec(i) +=
          Phi((xmin - Exp->X(i, j)) / Exp->SE(i, j),
              (xmax - Exp->X(i, j)) / Exp->SE(i, j));
      }
    }

  } else { // stepwise
    for (j = tmin - 1; j < tmax; ++j) {
      for (i = 0; i < Exp->n; ++i) {
        if ((Exp->X(i, j) >= xmin) && (Exp->X(i, j) < xmax))
          Xvec(i) += 1.0;
      }
    }
  }
  return (Xvec);
}

exposureDat::exposureDat(MatrixXd Tcalc_in) // Binomial DLM
{
  n = Tcalc_in.rows();
  pX = Tcalc_in.cols();
  pZ = 0;
  nSplits = 0;

  Tcalc = Tcalc_in;
  preset = 0;
}
exposureDat::exposureDat(MatrixXd Tcalc_in, // Gaussian DLM
                         MatrixXd Z_in,
                         MatrixXd Vg_in)
{
  n = Tcalc_in.rows();
  pX = Tcalc_in.cols();
  pZ = Z.cols();
  nSplits = 0;

  Tcalc = Tcalc_in;
  Z = Z_in;
  Vg = Vg_in;
  ZtTcalc = Z.transpose() * Tcalc;
  VgZtTcalc = Vg.selfadjointView<Lower>() * ZtTcalc;
  Z = Z_in;
  preset = 1;
}
exposureDat::exposureDat(MatrixXd X_in, // Binomial DLNM
                         MatrixXd SE_in,
                         VectorXd Xsplits_in,
                         MatrixXd Xcalc_in,
                         MatrixXd Tcalc_in,
                         bool lowmem_in)
{
  X = X_in;
  SE = SE_in;
  if (SE.rows() == X.rows())
    se = 1;
  else
    se = 0;
  n = X.rows();
  pX = X.cols();
  pZ = 0;

  Xsplits = Xsplits_in;
  nSplits = Xsplits.size();

  Xcalc = Xcalc_in;
  Tcalc = Tcalc_in;
  preset = 0;
  
  lowmem = lowmem_in;
  if (!lowmem) {
    for (int i = 0; i < nSplits; ++i) {
      MatrixXd Xtemp(n, pX); Xtemp.setZero();
      Xsave.push_back(Xtemp);
      for (int t = 0; t < pX; ++t)
        Xsave[i].col(t) = nodeCount(this, R_NegInf, Xsplits(i), 1, t + 1);
    }
    Xsave.push_back(Tcalc);
  }
}
exposureDat::exposureDat(MatrixXd X_in,  // Gaussian DLNM
                         MatrixXd SE_in,
                         VectorXd Xsplits_in,
                         MatrixXd Xcalc_in,
                         MatrixXd Tcalc_in,
                         MatrixXd Z_in,
                         MatrixXd Vg_in,
                         bool lowmem_in)
{
  X = X_in;
  Z = Z_in;
  Vg = Vg_in;
  SE = SE_in;
  if (SE.rows() == X.rows())
    se = 1;
  else
    se = 0;
  n = X.rows();
  pX = X.cols();
  pZ = Z.cols();

  Xsplits = Xsplits_in;
  nSplits = Xsplits.size();

  Xcalc = Xcalc_in;
  ZtXcalc = Z.transpose() * Xcalc;
  VgZtXcalc = Vg * ZtXcalc;

  Tcalc = Tcalc_in;
  ZtTcalc = Z.transpose() * Tcalc;
  VgZtTcalc = Vg * ZtTcalc;
  preset = 1;
  
  lowmem = lowmem_in;
  if (!lowmem) {
    for (int i = 0; i < nSplits; ++i) {
      MatrixXd Xtemp(n, pX);  Xtemp.setZero();
      Xsave.push_back(Xtemp);
      for (int t = 0; t < pX; ++t)
        Xsave[i].col(t) = nodeCount(this, R_NegInf, Xsplits(i), 1, t + 1);
      MatrixXd ZtXtemp = Z.transpose() * Xsave[i];
      ZtXsave.push_back(ZtXtemp);
      VgZtXsave.push_back(Vg * ZtXtemp);
    }
    Xsave.push_back(Tcalc);
    ZtXsave.push_back(ZtTcalc);
    VgZtXsave.push_back(VgZtTcalc);
  }
}
exposureDat::~exposureDat()
{
}

void exposureDat::updateNodeVals(Node *n)
{
  // stop if no update needed
  if (n->update == 0)
    return;

  // create NodeVals struct, if needed
  if (n->nodevals == 0)
    n->nodevals = new NodeVals(this->n, pZ);

  // update parent nodes first, if needed
  Node *sib, *parent;
  sib = 0; parent = 0;
  if ((n->depth > 0)) {
    parent = n->parent;
    sib = n->sib();
    if (sib == 0 || parent == 0)
      stop("missing node sib or parent");
    if (parent->update) // update parent
      updateNodeVals(parent);
    if (sib->nodevals == 0) // create NodeVals for sibling node
      sib->nodevals = new NodeVals(this->n, pZ);
    if ((sib->nodestruct->get(3) == 1) || // easier update
        (sib->nodestruct->get(1) == 0)) {
      n = sib;
      sib = n->sib();
    }
  }

  if (nSplits == 0) { // time splits only
    int tmin = n->nodestruct->get(3);
    int tmax = n->nodestruct->get(4);
    
    if (tmin == 1) { // if node begins at time 1
      n->nodevals->X = Tcalc.col(tmax - 1);
      if (preset) {
        n->nodevals->ZtX = ZtTcalc.col(tmax - 1);
        n->nodevals->VgZtX = VgZtTcalc.col(tmax - 1);
      }
    } else if (tmax == pX) { // if node ends at last time period
      n->nodevals->X = Tcalc.col(pX - 1) - Tcalc.col(tmin - 2);
      if (preset) {
        n->nodevals->ZtX = ZtTcalc.col(pX - 1) - ZtTcalc.col(tmin - 2);
        n->nodevals->VgZtX = VgZtTcalc.col(pX - 1) - VgZtTcalc.col(tmin - 2);
      }
    } else { // node begins/ends in middle of time periods
      n->nodevals->X = Tcalc.col(tmax - 1) - Tcalc.col(tmin - 2);
      if (preset) {
        n->nodevals->ZtX = ZtTcalc.col(tmax - 1) - ZtTcalc.col(tmin - 2);
        n->nodevals->VgZtX = VgZtTcalc.col(tmax-1) - VgZtTcalc.col(tmin-2);
      }
    }
    
  } else { // exposure and time splits in node
    int xmin = n->nodestruct->get(1);
    int xmax = n->nodestruct->get(2);
    int tmin = n->nodestruct->get(3);
    int tmax = n->nodestruct->get(4);
    
    if ((xmin == 0) && (xmax == nSplits + 1)) { // only splits in time
      int scale = tmax - tmin;
      n->nodevals->X = Tcalc.col(scale);
      if (preset) {
        n->nodevals->ZtX = ZtTcalc.col(scale);
        n->nodevals->VgZtX = VgZtTcalc.col(scale);
      }
    } else if ((tmin == 1) && (tmax == pX)) { // only splits in exposure
      if (xmin == 0) { // lower limint at min exposure
        n->nodevals->X = Xcalc.col(xmax - 1);
        if (preset) {
          n->nodevals->ZtX = ZtXcalc.col(xmax - 1);
          n->nodevals->VgZtX = VgZtXcalc.col(xmax - 1);
        }
      } else if (xmax == nSplits + 1) { // upper limit at max exposure
        n->nodevals->X = Tcalc.col(pX - 1) - Xcalc.col(xmin - 1);
        if (preset) {
          n->nodevals->ZtX = ZtTcalc.col(pX - 1) - ZtXcalc.col(xmin - 1);
          n->nodevals->VgZtX = VgZtTcalc.col(pX-1) - VgZtXcalc.col(xmin-1);
        }
      } else { // exposure limits between min/max
        n->nodevals->X = Xcalc.col(xmax - 1) - Xcalc.col(xmin - 1);
        if (preset) {
          n->nodevals->ZtX = ZtXcalc.col(xmax - 1) - ZtXcalc.col(xmin - 1);
          n->nodevals->VgZtX = VgZtXcalc.col(xmax - 1) - VgZtXcalc.col(xmin - 1);
        }
      }
      
    } else { // splits in both exposure and time, count manually
        
      if (!lowmem) { // if !lowmem use precalculated exposure counts
        n->nodevals->X = Xsave[xmax-1].col(tmax-1);
        if (preset) {
          n->nodevals->ZtX = ZtXsave[xmax-1].col(tmax-1);
          n->nodevals->VgZtX =  VgZtXsave[xmax-1].col(tmax-1);
        }
        if (xmin > 0) { // subtract xmin data
          n->nodevals->X -= Xsave[xmin-1].col(tmax-1);
          if (preset) {
            n->nodevals->ZtX -= ZtXsave[xmin-1].col(tmax-1);
            n->nodevals->VgZtX -=  VgZtXsave[xmin-1].col(tmax-1);
          }
        }
        if (tmin > 1) { // if tmin > 1, remove lower time periods from count
            n->nodevals->X -= Xsave[xmax-1].col(tmin-2);
            if (preset) {
              n->nodevals->ZtX -= ZtXsave[xmax-1].col(tmin-2);
              n->nodevals->VgZtX -= VgZtXsave[xmax-1].col(tmin-2);
            }
          if (xmin > 0) { // subtract xmin data
            n->nodevals->X += Xsave[xmin-1].col(tmin-2);
            if (preset) {
              n->nodevals->ZtX += ZtXsave[xmin-1].col(tmin-2);
              n->nodevals->VgZtX += VgZtXsave[xmin-1].col(tmin-2);
            }
          }
        } // end if tmin != 1
      } else {
        double xmin2 = R_NegInf;
        double xmax2 = R_PosInf;
        if (xmin != 0)
          xmin2 = Xsplits(xmin - 1);
        if (xmax != nSplits + 1)
          xmax2 = Xsplits(xmax - 1);

        n->nodevals->X = nodeCount(this, xmin2, xmax2, tmin, tmax);
        if (preset) {
          n->nodevals->ZtX = Z.transpose() * n->nodevals->X;
          n->nodevals->VgZtX = Vg * n->nodevals->ZtX;
        }
      } // end if !lowmem
    }
  }

  if (n->depth > 0) { // update sibling node
    if (sib->update) {
      (sib->nodevals)->X = (parent->nodevals)->X - n->nodevals->X;
      if (preset) {
        (sib->nodevals)->ZtX = (parent->nodevals)->ZtX - n->nodevals->ZtX;
        (sib->nodevals)->VgZtX = 
          (parent->nodevals)->VgZtX - n->nodevals->VgZtX;
      }
      sib->update = 0;
    }
  }
  n->update = 0;
  return;
}
