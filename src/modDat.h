#define MODDAT_H
#include <Rcpp.h>
#include <RcppEigen.h>

class Node;

class modDat {
public:
  modDat(std::vector<int> varIsNum, 
         Rcpp::List splitIdx, 
         std::vector<int> fullIdx);
  ~modDat();

  int nMods, n;  // number of modifiers, data size
  std::vector<int> varIsNum; // index of continuous modifiers
  std::vector<int> nModSplit; // index of number of available splits
  Eigen::VectorXd modProb;    // probability of selecting each modifier
  std::vector<std::vector<int> > availMod; // list of modifier splits
  std::vector<std::vector<std::vector<int> > > splitIdx;
      // vectors of split indices (1.modifier, 2.split, 3.index of observations)
  std::vector<int> fullIdx;
  double totalProb(std::vector<std::vector<int> >);
      // total probability of available modifiers
  std::vector<std::vector<int> > getAvailMods(int, int, std::vector<int>,
                                             std::vector<std::vector<int> >, bool);
      // get new list of available modifiers for subStructs


  void updateNodeVals(Node*); // update node indices based on rules
};
