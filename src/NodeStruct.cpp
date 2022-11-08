#include "RcppEigen.h"
#include "NodeStruct.h"
#include "modDat.h"
#include "Fncs.h"
using namespace Rcpp;

// Generic function placeholders
NodeStruct::NodeStruct() {}
NodeStruct::~NodeStruct() {}
NodeStruct::NodeStruct(const NodeStruct&) {}
NodeStruct* NodeStruct::clone() {return(0);}
NodeStruct* NodeStruct::subStruct(bool b) {return(0);}
bool NodeStruct::valid() {return(0);}
void NodeStruct::updateStruct(NodeStruct* a, bool b) {}
double NodeStruct::logPRule() {return(0);}
bool NodeStruct::proposeSplit() {return(0);}
void NodeStruct::dropSplit() {}
void NodeStruct::printStruct() {}
int NodeStruct::get(int a) {return(0);}
std::vector<int> NodeStruct::get2(int a)
  { std::vector<int> b; return(b); }
std::vector<std::vector<int> > NodeStruct::get3(int a)
  { std::vector<std::vector<int> > b; return(b); }
bool NodeStruct::checkEqual(NodeStruct* n) {return(0);}
void NodeStruct::setTimeRange(int lower, int upper) {}
void NodeStruct::setTimeProbs(Eigen::VectorXd newProbs) {}
Eigen::VectorXd NodeStruct::getTimeProbs() {}

/**
 * @brief Construct a new DLNMStruct::DLNMStruct object
 * 
 * @param xmin_in set to 0 for new tree
 * @param xmax_in set to 1 for new DLM tree or larger for DLNM
 * @param tmin_in set to 1 for new tree
 * @param tmax_in set to number of exposure measurements
 * @param Xp_in probability of selecting expsoure splitting locations
 * @param Tp_in 
 */
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


Eigen::VectorXd DLNMStruct::getTimeProbs()  { return(this->Tp); }
void DLNMStruct::setTimeProbs(Eigen::VectorXd newProbs) { Tp = newProbs; }




bool DLNMStruct::proposeSplit()
{
  if (((xmin >= xmax - 1) && (tmin >= tmax)) ||
      ((totXp <= 0.0) && (totTp <= 0.0))) {
    return(0);
  }
  if ((xmin >= xmax - 1) || (totXp <= 0.0)) { // sample T
    tsplit = sampleInt(Tp.segment(tmin - 1, tmax - tmin)) + tmin;
    xsplit = 0;

  } else if ((tmin >= tmax) || (totTp <= 0.0)) { // sample X
    xsplit = sampleInt(Xp.segment(xmin, xmax - xmin - 1)) + xmin + 1;
    tsplit = 0;

  } else { // sample X or T
    if (R::runif(0, 1) < (totXp / (totXp + totTp))) {
      xsplit = sampleInt(Xp.segment(xmin, xmax - xmin - 1)) + xmin + 1;
      tsplit = 0;

    } else {
      tsplit = sampleInt(Tp.segment(tmin - 1, tmax - tmin)) + tmin;
      xsplit = 0;

    }
  }

  return(1);
}

void DLNMStruct::dropSplit()
{
  xsplit = 0;
  tsplit = 0;
}

bool DLNMStruct::valid()
{
  if ((xmin >= xmax) || (tmin > tmax)) {
    return(0);
  }
  return(1);
}

bool DLNMStruct::checkEqual(NodeStruct* n) {return(0);}

void DLNMStruct::updateStruct(NodeStruct* parStruct, bool left)
{
  Tp = parStruct->getTimeProbs();
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
  totXp = -1.0;
  totTp = -1.0;
  if (xmax > xmin + 1)
    totXp = Xp.segment(xmin, xmax - xmin - 1).sum();
  if (tmax > tmin)
    totTp = Tp.segment(tmin - 1, tmax - tmin).sum();
}

void DLNMStruct::setTimeRange(int lower, int upper)
{
  tmin = lower;
  tmax = upper;
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

  return(ns);
}


NodeStruct* DLNMStruct::clone() { return new DLNMStruct(*this); }



double DLNMStruct::logPRule()
{
  if (xsplit != 0) {
    return(log(Xp(xsplit - 1) / (totXp + totTp)));
  } else if (tsplit != 0) {
    return(log(Tp(tsplit - 1) / (totXp + totTp)));
  } else {
    return(0);
  }
}

void DLNMStruct::printStruct()
{
  Rcout << "Struct:" << "x in [" << xmin << "," << xmax << "] split at " <<
    xsplit << ", t in [" <<
      tmin << "," << tmax << "] split at " << tsplit << ", logPRule = " <<
        logPRule() << ", totXp = " << totXp << ", totTp = " << totTp << "\n";
}

int DLNMStruct::get(int a)
{
  switch (a) {
    case 1: return(xmin);
    case 2: return(xmax);
    case 3: return(tmin);
    case 4: return(tmax);
    case 5: return(xsplit);
    case 6: return(tsplit);
  }
  stop("incorrect call to DLNMStruct::get");
}






ModStruct::ModStruct(modDat* md)
{
  modFncs = md;
  availMod = md->availMod;
  splitVal = -1;
  splitVar = -1;
}

ModStruct::ModStruct(modDat* md, std::vector<std::vector<int> > am)
{
  modFncs = md;
  availMod = am;
  splitVal = -1;
  splitVar = -1;
}

// Delete (do not delete modDat object, set pointer to zero)
ModStruct::~ModStruct()
{
  modFncs = 0;
}

// Clone function
ModStruct::ModStruct(const ModStruct& ns)
{
  modFncs = ns.modFncs;
  availMod = ns.availMod;
  splitVal = ns.splitVal;
  splitVar = ns.splitVar;
  splitVec = ns.splitVec;
}

// Proposed new split based on availMod
bool ModStruct::proposeSplit()
{
  std::size_t i;

  // Determine which modifiers are available for split
  std::vector<int> whichAvail;
  std::vector<double> modProbAvail;
  double totProbAvail = 0;
  for (i = 0; i < availMod.size(); ++i) {
    if ((availMod[i].size() > 0) && (modFncs->modProb(i) > 0)) {
      whichAvail.push_back(i);
      modProbAvail.push_back(modFncs->modProb(i));
      totProbAvail += modFncs->modProb(i);
    }
  }
  if ((whichAvail.size() == 0) || (totProbAvail == 0))
    return(0);

  // Cycle through available modifiers until split is found
  for (i = 0; i < whichAvail.size(); ++i) {
    splitVar = whichAvail[sampleInt(modProbAvail, totProbAvail)];
    std::vector<int> am = availMod[splitVar];

    // Continuous split
    if (modFncs->varIsNum[splitVar]) {
      splitVal = am[floor(R::runif(0, am.size() - 1))];
      // Rcout << "\nSel:" << splitVar << " " << splitVal;
      return(1);

    // Categorical split
    } else {
      if (am.size() == 1)
        return(0);
      std::vector<int> unavailMod;
      std::size_t j = 0;
      for (i = 0; i < (std::size_t) modFncs->nModSplit[splitVar]; ++i) {
        if (i != am[j]) {
          unavailMod.push_back(i);
        } else if (j < (am.size() - 1)) {
          ++j;
        }
      }

      // Select random categories from available
      splitVec.clear();
      int nCat = floor(R::runif(1.0, am.size()));
      std::random_shuffle(am.begin(), am.end());
      for (i = 0; i < (std::size_t) nCat; ++i) {
        splitVec.push_back(am[i]);
      }

      // Select random categories from unavailable
      if (unavailMod.size() > 0) {
        int nOther = floor(R::runif(0.0, unavailMod.size() + 1.0));
        if (nOther > 0) {
          std::random_shuffle(unavailMod.begin(), unavailMod.end());
          for (i = 0; i < (std::size_t) nOther; ++i) {
            splitVec.push_back(unavailMod[i]);
          }
        }
      }
      if ((splitVec.size() == 0) || 
          (splitVec.size() == modFncs->nModSplit[splitVar]))
        return(0);
        
      std::sort(splitVec.begin(), splitVec.end());
      return(1);
    } // end categorical split
  }
  return(0);
}

// Reset split
void ModStruct::dropSplit()
{
  splitVal = -1;
  splitVar = -1;
  splitVec.clear();
}

// Check that split is valid
bool ModStruct::valid()
{
  if (splitVar == -1)
    return(1);
  if (availMod[splitVar].size() == 0)
    return(0);
  if (splitVal == -1) {
    // std::vector<std::vector<int> > 
    std::sort(splitVec.begin(), splitVec.end());
    std::vector<int> v_intersection;
    std::set_intersection(availMod[splitVar].begin(), availMod[splitVar].end(),
                          splitVec.begin(), splitVec.end(),
                          std::back_inserter(v_intersection));
    // std::pair<std::vector<int>, std::vector<int> > iD =
    //   intersectAndDiff(availMod[splitVar], splitVec);
    if (v_intersection.size() == 0)
      return(0);
  } else {
    for (int i : availMod[splitVar]) {
      if (i == splitVal)
        return(1);
    }
    return(0);
  }
  return(1);
}

bool ModStruct::checkEqual(NodeStruct* ns)
{
  if ((splitVar == ns->get(1)) &&
      (splitVal == ns->get(2)) &&
      (splitVec == ns->get2(1))) {
    return(1);
  } else {
    return(0);
  }
}

void ModStruct::updateStruct(NodeStruct* parStruct, bool left)
{
  availMod = modFncs->getAvailMods(parStruct->get(1),
                                   parStruct->get(2),
                                   parStruct->get2(1),
                                   parStruct->get3(1), left);
}




NodeStruct* ModStruct::subStruct(bool left)
{
  // Rcout << "makeSubSt";
  std::vector<std::vector<int> > am =
    modFncs->getAvailMods(splitVar, splitVal, splitVec, availMod, left);
  NodeStruct *ns;
  ns = new ModStruct(modFncs, am);
  return(ns);
}

NodeStruct* ModStruct::clone() { return new ModStruct(*this); }


double ModStruct::logPRule()
{
  if (splitVar == -1)
    return(0);

  if (modFncs->varIsNum[splitVar]) { // Continuous split
    return(log(modFncs->modProb(splitVar)) -
           log(modFncs->totalProb(availMod)) -
           log(availMod[splitVar].size()));

  } else { // Categorical split
    return(log(modFncs->modProb(splitVar)) -
           log(modFncs->totalProb(availMod)) -
           log(pow(2.0, double(availMod[splitVar].size()) - 1.0) - 1.0));
  }

}

int ModStruct::get(int a)
{
  switch(a) {
    case 1: return(splitVar);
    case 2: return(splitVal);
  }
  stop("incorrect call to ModStruct::get");
}

std::vector<int> ModStruct::get2(int a)
{
  switch(a) {
    case 1: return(splitVec);
  }
  stop("incorrect call to ModStruct::get2");
}

std::vector<std::vector<int> > ModStruct::get3(int a)
{
  switch(a) {
    case 1: return(availMod);
  }
  stop("incorrect call to ModStruct::get3");
}


void ModStruct::printStruct()
{
  Rcout << "\nStruct: splitVar = " << splitVar <<
    " splitVal = " << splitVal << " splitVec = ";
  for (int i : splitVec)
    Rcout << i << " ";
}
