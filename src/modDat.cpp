#include <Rcpp.h>
#include <RcppEigen.h>
#include "modDat.h"
#include "Fncs.h"
#include "Node.h"
#include "NodeStruct.h"
using namespace Rcpp;

modDat::modDat(std::vector<int> isnum,
               Rcpp::List spIdx,
               std::vector<int> fidx)
{
  varIsNum = isnum;
  nMods = varIsNum.size();
  fullIdx = fidx;
  n = fullIdx.size();
  for (int i = 0; i < nMods; ++i) {
    std::vector<std::vector<int> > temp;
    std::vector<int> temp2;
    Rcpp::List splits = as<Rcpp::List>(spIdx[i]);
    for (int j = 0; j < splits.size(); ++j) {
      temp.push_back(as<std::vector<int> >(splits[j]));
      temp2.push_back(j);
    }
    splitIdx.push_back(temp);
    availMod.push_back(temp2);
    nModSplit.push_back(temp2.size());
  }
  modProb.resize(nMods); modProb.setOnes();
  modProb /= nMods;
}

modDat::~modDat() {}

double modDat::totalProb(std::vector<std::vector<int> > am)
{
  double tp = 0;
  for (int i = 0; i < nMods; ++i) {
    if (am[i].size() > 0)
      tp += modProb(i);
  }
  return(tp);
}

std::vector<std::vector<int> >
  modDat::getAvailMods(int splitVar, int splitVal, std::vector<int> splitVec,
                       std::vector<std::vector<int> > am, bool left)
{
  std::vector<std::vector<int> > newAvailMod = am;
  std::size_t i;

  if (splitVar == -1)
    return(newAvailMod);

  if (splitVal != -1) {
    std::vector<int> idx;
    if (left) {
      for (i = 0; i < newAvailMod[splitVar].size(); ++i) {
        if (splitVal > newAvailMod[splitVar][i]) {
          idx.push_back(newAvailMod[splitVar][i]);
        }
      }
    } else {
      for (i = 0; i < newAvailMod[splitVar].size(); ++i) {
        if (splitVal < newAvailMod[splitVar][i]) {
          idx.push_back(newAvailMod[splitVar][i]);
        }
      }
    }
    newAvailMod[splitVar] = idx;

  } else {
    std::pair<std::vector<int>, std::vector<int> > iD =
      intersectAndDiff(newAvailMod[splitVar], splitVec);
    if (left) {
      newAvailMod[splitVar] = iD.first;
    } else {
      newAvailMod[splitVar] = iD.second;
    }
  }
  return(newAvailMod);
}


void modDat::updateNodeVals(Node *n)
{
  // stop if no update needed
  if (n->update == 0)
    return;
    
  // load NodeVals struct
  if (n->nodevals == 0)
    n->nodevals = new NodeVals(this->n);

  // update parent nodes first, if needed
  Node *sib, *parent;
  sib = 0; parent = 0;

  if (n->depth > 0) {
    parent = n->parent;
    sib = n->sib();

    if (parent->update) {
      updateNodeVals(parent);
      parent->nodevals->updateXmat = 1;
    }

    if (sib->nodevals == 0)
      sib->nodevals = new NodeVals(this->n);

  } else {
    n->nodevals->idx = this->fullIdx;
    n->update = 0;
    n->nodevals->updateXmat = 1;
    return;
  }

  if (parent->nodevals->idx.size() == 0){
    n->nodevals->idx = parent->nodevals->idx;
    sib->nodevals->idx = parent->nodevals->idx;
    n->update = 0;
    n->nodevals->updateXmat = 1;
    sib->update = 0;
    sib->nodevals->updateXmat = 1;
    return;
  }

  int splitVar = (parent->nodestruct)->get(1);
  
  // Find intersection and difference of parent and rule indices
  std::pair<std::vector<int>, std::vector<int> > iD;
  if (varIsNum[splitVar]) {
    int splitVal = (parent->nodestruct)->get(2);
    iD = intersectAndDiff((parent->nodevals)->idx,
                          splitIdx[splitVar][splitVal]);
  } else {
    std::vector<int> splitVec = (parent->nodestruct)->get2(1);
    int vecLen = 0;
    for (auto i : splitVec)
      vecLen += splitIdx[splitVar][i].size();
    std::vector<int> splitVecIdx;
    splitVecIdx.reserve(vecLen);
    for (auto i : splitVec)
      splitVecIdx.insert(splitVecIdx.end(),
                         splitIdx[splitVar][i].begin(),
                         splitIdx[splitVar][i].end());
    std::sort(splitVecIdx.begin(), splitVecIdx.end());
    iD = intersectAndDiff((parent->nodevals)->idx, splitVecIdx);
  }

  // Assign indices to node and sibling
  if (n == parent->c1) {
    (n->nodevals)->idx = iD.first;
    (sib->nodevals)->idx = iD.second;
  } else if (n == parent->c2) {
    (n->nodevals)->idx = iD.second;
    (sib->nodevals)->idx = iD.first;
  } else if (parent->proposed != 0) {
    if (n == (parent->proposed)->c1) {
      (n->nodevals)->idx = iD.first;
      (sib->nodevals)->idx = iD.second;
    } else {
      (n->nodevals)->idx = iD.second;
      (sib->nodevals)->idx = iD.first;
    }
  }

  n->update = 0;
  n->nodevals->updateXmat = 1;
  sib->update = 0;
  sib->nodevals->updateXmat = 1;
}
