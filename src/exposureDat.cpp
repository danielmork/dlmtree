#include <RcppEigen.h>
#include "exposureDat.h"
#include "Node.h"
#include "NodeStruct.h"
using namespace Rcpp;

#define DOUBLE_ERROR   0.000000000000000111022302462515654042363166809082031259999
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
// #define R_PosInf            1569325055
// #define R_NegInf         -1569325055

void nodeCount(Node*, exposureDat*, double, double, int, int);
double Phi(double, double);

exposureDat::exposureDat(Eigen::MatrixXd Tcalc_in)
{
  n = Tcalc_in.rows();
  pX = Tcalc_in.cols();
  pZ = 0;
  nSplits = 0;

  Tcalc = Tcalc_in;
  preset = 0;
}
exposureDat::exposureDat(Eigen::MatrixXd Tcalc_in,
                         Eigen::MatrixXd Z_in,
                         Eigen::MatrixXd Vg_in)
{
  n = Tcalc_in.rows();
  pX = Tcalc_in.cols();
  pZ = Z.cols();
  nSplits = 0;

  Tcalc = Tcalc_in;
  Z = Z_in;
  Vg = Vg_in;
  ZtTcalc = Z.transpose() * Tcalc;
  VgZtTcalc = Vg.selfadjointView<Eigen::Lower>() * ZtTcalc;
  Z = Z_in;
  preset = 1;
}
exposureDat::exposureDat(Eigen::MatrixXd X_in,
                         Eigen::MatrixXd SE_in,
                         Eigen::VectorXd Xsplits_in,
                         Eigen::MatrixXd Xcalc_in,
                         Eigen::MatrixXd Tcalc_in)
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
}
exposureDat::exposureDat(Eigen::MatrixXd X_in,
                         Eigen::MatrixXd SE_in,
                         Eigen::VectorXd Xsplits_in,
                         Eigen::MatrixXd Xcalc_in,
                         Eigen::MatrixXd Tcalc_in,
                         Eigen::MatrixXd Z_in,
                         Eigen::MatrixXd Vg_in)
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
}
exposureDat::~exposureDat()
{
}

void exposureDat::updateNodeVals(Node *n)
{
  // stop if no update needed
  if (n->update == 0)
    return;

  // load NodeVals struct
  if (n->nodevals == 0) {
    // Rcout << "\nmake nodevals\n";
    n->nodevals = new NodeVals(this->n, pZ);
  }

  // update parent nodes first, if needed
  Node *sib, *parent;
  sib = 0; parent = 0;
  if ((n->depth > 0)) {
    // Rcout << "\ndepth>0\n";
    parent = n->parent;
    // Rcout << parent << "\n";
    sib = n->sib();
    // Rcout << sib << "\n";
    if (sib == 0 || parent == 0)
      stop("missing node sib or parent");

    if (parent->update) {
      // Rcout << "\nupdate parent\n";
      updateNodeVals(parent);
    }
    if (sib->nodevals == 0) {
      // Rcout << "\nmake sib nodevals\n";
      sib->nodevals = new NodeVals(this->n, pZ);
    }
  }

  // time splits only
  if (nSplits == 0) {
    // Rcout << "\n0 splits\n";
    // (n->nodestruct)->printStruct();
    int tmin = (n->nodestruct)->get(3);
    int tmax = (n->nodestruct)->get(4);

    // top of tree
    // if (n->depth == 0) {
    //   // Rcout << "\nupdate depth 0\n";
    //   (n->nodevals)->X = Tcalc.col(pX - 1);
    //   if (preset) {
    //     (n->nodevals)->ZtX = ZtTcalc.col(pX - 1);
    //     (n->nodevals)->VgZtX = VgZtTcalc.col(pX - 1);
    //   }
    //   n->update = 0;
    //   return;

    //   // starting at time 1
    // } else 
    if (tmin == 1) {
      // Rcout << "\nupdate tmin 1\n";
      (n->nodevals)->X = Tcalc.col(tmax - 1);
      if (preset) {
        (n->nodevals)->ZtX = ZtTcalc.col(tmax - 1);
        (n->nodevals)->VgZtX = VgZtTcalc.col(tmax - 1);
      }

      // ending at last time
    } else if (tmax == pX) {
      // Rcout << "\nupdate tmax px\n";
      (n->nodevals)->X = Tcalc.col(pX - 1) - Tcalc.col(tmin - 2);
      if (preset) {
        (n->nodevals)->ZtX = ZtTcalc.col(pX - 1) - ZtTcalc.col(tmin - 2);
        (n->nodevals)->VgZtX = VgZtTcalc.col(pX - 1) - VgZtTcalc.col(tmin - 2);
      }

      // start/end times in middle
    } else {
      // Rcout << "\nupdate other\n";
      (n->nodevals)->X = Tcalc.col(tmax - 1) - Tcalc.col(tmin - 2);
      if (preset) {
        (n->nodevals)->ZtX = ZtTcalc.col(tmax - 1) - ZtTcalc.col(tmin - 2);
        (n->nodevals)->VgZtX = VgZtTcalc.col(tmax - 1) - VgZtTcalc.col(tmin - 2);
      }
    }

    // exposure and time splits
  } else {
    // (n->nodestruct)->printStruct();
    int xmin = (n->nodestruct)->get(1);
    // if (xmin == 0)
    // xmin = R_NegInf;
    int xmax = (n->nodestruct)->get(2);
    // if (xmax == nSplits + 1)
    // xmax = R_PosInf;
    int tmin = (n->nodestruct)->get(3);
    int tmax = (n->nodestruct)->get(4);
    // Rcout << xmin << " " << xmax << " " << tmin << " " << tmax << "\n";

    // top of tree
    // if (n->depth == 0) {
    //   // Rcout << "update top";
    //   (n->nodevals)->X = Tcalc.col(pX - 1);
    //   if (preset) {
    //     (n->nodevals)->ZtX = ZtTcalc.col(pX - 1);
    //     (n->nodevals)->VgZtX = VgZtTcalc.col(pX - 1);
    //   }
    //   n->update = 0;
    //   return;

    //   // first split on time
    // } else 
    if ((xmin == 0) && (xmax == nSplits + 1)) {
      // Rcout << "update time";
      int scale = tmax - tmin;
      (n->nodevals)->X = Tcalc.col(scale);
      if (preset) {
        (n->nodevals)->ZtX = ZtTcalc.col(scale);
        (n->nodevals)->VgZtX = VgZtTcalc.col(scale);
      }

      // first split on exposure
    } else if ((tmin == 1) && (tmax == pX)) {
      // Rcout << "update exposure ";
      if (xmin == 0) {
        // Rcout << "xmin inf";
        (n->nodevals)->X = Xcalc.col(xmax - 1);
        if (preset) {
          (n->nodevals)->ZtX = ZtXcalc.col(xmax - 1);
          (n->nodevals)->VgZtX = VgZtXcalc.col(xmax - 1);
        }
      } else if (xmax == nSplits + 1) {
        // Rcout << "xmax inf";
        (n->nodevals)->X = Tcalc.col(pX - 1) - Xcalc.col(xmin - 1);
        if (preset) {
          (n->nodevals)->ZtX = ZtTcalc.col(pX - 1) - ZtXcalc.col(xmin - 1);
          (n->nodevals)->VgZtX = VgZtTcalc.col(pX - 1) - VgZtXcalc.col(xmin - 1);
        }
      } else {
        // Rcout << "other";
        (n->nodevals)->X = Xcalc.col(xmax - 1) - Xcalc.col(xmin - 1);
        if (preset) {
          (n->nodevals)->ZtX = ZtXcalc.col(xmax - 1) - ZtXcalc.col(xmin - 1);
          (n->nodevals)->VgZtX = VgZtXcalc.col(xmax - 1) - VgZtXcalc.col(xmin - 1);
        }
      }

      // count manually
    } else {
      // Rcout << "update manual";
      double xmin2 = R_NegInf;
      double xmax2 = R_PosInf;
      if (xmin != 0)
        xmin2 = Xsplits(xmin - 1);

      if (xmax != nSplits + 1)
        xmax2 = Xsplits(xmax - 1);

      ((n->nodevals)->X).resize(this->n); ((n->nodevals)->X).setZero();
      nodeCount(n, this, xmin2, xmax2, tmin, tmax);
      // (n->nodevals)->X /= (sqrt(this->n) + double (pX));
      if (preset) {
        // Rcout << "ZtX, VgZtX";
        (n->nodevals)->ZtX = Z.transpose() * (n->nodevals)->X;
        (n->nodevals)->VgZtX = Vg * (n->nodevals)->ZtX;
      }
    }
  }

  if (parent == 0) {
    n->update = 0;
    return;
  }

  // update sibling node
  if (sib->update) {
    (sib->nodevals)->X = (parent->nodevals)->X - (n->nodevals)->X;
    if (preset) {
      (sib->nodevals)->ZtX = (parent->nodevals)->ZtX - (n->nodevals)->ZtX;
      (sib->nodevals)->VgZtX = (parent->nodevals)->VgZtX - (n->nodevals)->VgZtX;
    }
    sib->update = 0;
  }
  n->update = 0;
  // Rcout << " done\n";
  return;
}

double Phi(double x1, double x2)
{
  return (erf(x2 * MATH_SQRT1_2) - erf(x1 * MATH_SQRT1_2)) * 0.5;
}

void nodeCount(Node* n, exposureDat* Exp,
               double xmin, double xmax, int tmin, int tmax)
{
  int i, j;
  // double inc = 1.0 / (sqrt(Exp->n)* double(Exp->pX));//

  if (Exp->se) {      
    for (j = tmin - 1; j < tmax; ++j) {
      for (i = 0; i < Exp->n; ++i) {
        (n->nodevals)->X(i) +=
          Phi((xmin - Exp->X(i, j)) / Exp->SE(i, j),
              (xmax - Exp->X(i, j)) / Exp->SE(i, j));// * inc;
      }
    }

  } else {
    // auto f = [&xmin, &xmax](double x) {
    //   return double ((x >= xmin) && (x < xmax));
    // };
    // n->nodevals->X = Exp->X.block(0, tmin - 1, Exp->n, tmax - tmin + 1).unaryExpr(f).array().rowwise().sum();
    for (j = tmin - 1; j < tmax; ++j) {
      for (i = 0; i < Exp->n; ++i) {
        if ((Exp->X(i, j) >= xmin) && (Exp->X(i, j) < xmax))
          ((n->nodevals)->X(i)) += 1.0;//inc;
      }
    }
  }
}
