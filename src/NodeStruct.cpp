#include <RcppEigen.h>
#include <Rcpp.h>
#include "NodeStruct.h"
#include "Fncs.h"
using namespace Rcpp;

NodeStruct::NodeStruct() {}
NodeStruct::~NodeStruct() {}
NodeStruct::NodeStruct(const NodeStruct&) {}
NodeStruct* NodeStruct::clone() {return(0);}
NodeStruct* NodeStruct::subStruct(bool b) {return(0);}
bool NodeStruct::valid() {return(0);}
void NodeStruct::updateStruct(NodeStruct* a, bool b) {};
double NodeStruct::logPRule() {return(0);}
bool NodeStruct::proposeSplit() {return(0);}
void NodeStruct::dropSplit() {}
void NodeStruct::printStruct() {}
int NodeStruct::get(int a) {return(0);}

DLNMStruct::DLNMStruct(int xmin_in, int xmax_in, int tmin_in, int tmax_in,
                       Eigen::VectorXd Xp_in, Eigen::VectorXd Tp_in)
{
  xmin = xmin_in; xmax = xmax_in; tmin = tmin_in; tmax = tmax_in;
  Xp = Xp_in; Tp = Tp_in;
  xsplit = 0; tsplit = 0;
  totXp = Xp.segment(xmin, xmax - xmin - 1).array().sum();
  totTp = Tp.segment(tmin - 1, tmax - tmin).array().sum();
}

DLNMStruct::~DLNMStruct() {}

DLNMStruct::DLNMStruct(const DLNMStruct& ns)
{
  xmin = ns.xmin; xmax = ns.xmax; tmin = ns.tmin; tmax = ns.tmax;
  Xp = ns.Xp; Tp = ns.Tp; xsplit = ns.xsplit; tsplit = ns.tsplit;
  totXp = Xp.segment(xmin, xmax - xmin - 1).array().sum();
  totTp = Tp.segment(tmin - 1, tmax - tmin).array().sum();
}

NodeStruct* DLNMStruct::clone() { return new DLNMStruct(*this); }


void DLNMStruct::updateStruct(NodeStruct* parStruct, bool left)
{
  int xs = parStruct->get(5);
  int ts = parStruct->get(6);
  if (left) {
    if (xs != 0) {
      xmin = parStruct->get(1);
      xmax = xs;
      tmin = parStruct->get(3);
      tmax = parStruct->get(4);
    } else if (ts != 0) {
      xmin = parStruct->get(1);
      xmax = parStruct->get(2);
      tmin = parStruct->get(3);
      tmax = ts;
    }
  } else {
    if (xs != 0) {
      xmin = xs;
      xmax = parStruct->get(2);
      tmin = parStruct->get(3);
      tmax = parStruct->get(4);
    } else if (ts != 0) {
      xmin = parStruct->get(1);
      xmax = parStruct->get(2);
      tmin = ts + 1;
      tmax = parStruct->get(4);
    }
  }
  totXp = Xp.segment(xmin, xmax - xmin - 1).array().sum();
  totTp = Tp.segment(tmin - 1, tmax - tmin).array().sum();
}

bool DLNMStruct::valid()
{
  if ((xmin >= xmax) || (tmin > tmax)) {
    return(0);
  }
  // if (xsplit != 0) {
  //   if ((xsplit >= xmax) || (xsplit <= xmin))
  //     return(0);
  // }
  // if (tsplit != 0) {
  //   if ((tsplit >= tmax) || (tsplit < tmin))
  //     return(0);
  // }
  return(1);
}

double DLNMStruct::logPRule()
{
  // totXp = Xp.segment(xmin, xmax - xmin - 1).array().sum();
  // totTp = Tp.segment(tmin - 1, tmax - tmin).array().sum();
  if (xsplit != 0) {
    return(log(Xp(xsplit - 1) / (totXp + totTp)));
  } else if (tsplit != 0) {
    return(log(Tp(tsplit - 1) / (totXp + totTp)));
  } else {
    return(0);
  }
}

bool DLNMStruct::proposeSplit()
{
  if ((xmin >= xmax - 1) && (tmin >= tmax)) {
    return(0);
  }

  // totXp = Xp.segment(xmin, xmax - xmin - 1).sum();
  // totTp = Tp.segment(tmin - 1, tmax - tmin).sum();

  int split;
  if ((xmin >= xmax - 1) || (totXp == 0)) { // sample T
    tsplit = sampleInt(Tp.segment(tmin - 1, tmax - tmin)) + tmin;
    xsplit = 0;
    // if ((tsplit < tmin) || (tsplit >= tmax)) {
    //   return(0);
    // }

  } else if ((tmin >= tmax) || (totTp == 0)) { // sample X
    xsplit = sampleInt(Xp.segment(xmin, xmax - xmin - 1)) + xmin + 1;
    tsplit = 0;
    // if ((xsplit <= xmin) || (xsplit >= xmax)) {
    //   return(0);
    // }

  } else { // sample X or T
    if (R::runif(0, 1) < (totXp / (totXp + totTp))) {
      xsplit = sampleInt(Xp.segment(xmin, xmax - xmin - 1)) + xmin + 1;
      tsplit = 0;
      // if ((xsplit <= xmin) || (xsplit >= xmax)) {
      //   return(0);
      // }
    } else {
      tsplit = sampleInt(Tp.segment(tmin - 1, tmax - tmin)) + tmin;
      xsplit = 0;
      // if ((tsplit < tmin) || (tsplit >= tmax)) {
      //   return(0);
      // }
    }
  }

  return(1);
}

NodeStruct* DLNMStruct::subStruct(bool left)
{
  NodeStruct *ns;
  if (left) {
    if (xsplit > 0) {
      ns = new DLNMStruct(xmin, xsplit, tmin, tmax, Xp, Tp);

    } else {
      ns = new DLNMStruct(xmin, xmax, tmin, tsplit, Xp, Tp);
    }

  } else {
    if (xsplit > 0) {
      ns = new DLNMStruct(xsplit, xmax, tmin, tmax, Xp, Tp);

    } else {
      ns = new DLNMStruct(xmin, xmax, tsplit + 1, tmax, Xp, Tp);
    }
  }
  // if (n0 != 0) {
  //   ns->setSplit(n0);
  //   // ns->xsplit = n0->get(5);
  //   // ns->tsplit = n0->get(6);
  // }

  return(ns);
}

void DLNMStruct::dropSplit()
{
  xsplit = 0;
  tsplit = 0;
}

void DLNMStruct::printStruct()
{
  Rcout << "Struct:" << "x in [" << xmin << "," << xmax << "] split at " <<
    xsplit << ", t in [" <<
      tmin << "," << tmax << "] split at " << tsplit << ", logPRule = " <<
        logPRule() << "\n";
}

int DLNMStruct::get(int a)
{
  switch(a) {
  case 1: return(xmin);
  case 2: return(xmax);
  case 3: return(tmin);
  case 4: return(tmax);
  case 5: return(xsplit);
  case 6: return(tsplit);
  }
  return(0);
}
