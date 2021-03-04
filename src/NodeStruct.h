#define NODESTRUCT_H
#include <RcppEigen.h>
#include <Rcpp.h>

class modDat;

class NodeStruct {
public:
  NodeStruct();
  virtual ~NodeStruct();
  NodeStruct(const NodeStruct&);
  virtual NodeStruct* clone();
  virtual NodeStruct* subStruct(bool);
  virtual bool valid();
  virtual double logPRule();
  virtual bool proposeSplit();
  virtual void dropSplit();
  virtual void printStruct();
  virtual int get(int);
  virtual std::vector<int> get2(int);
  virtual std::vector<std::vector<int> > get3(int);
  virtual void updateStruct(NodeStruct*, bool);
  virtual bool checkEqual(NodeStruct*);
};

class DLNMStruct: public NodeStruct {
public:
  DLNMStruct(int, int, int, int, Eigen::VectorXd, Eigen::VectorXd);
  ~DLNMStruct();
  DLNMStruct(const DLNMStruct&);

  int xmin, xmax, tmin, tmax; // DLNM limits
  int xsplit, tsplit; // set to 0 until split is created
  Eigen::VectorXd Xp, Tp; // probabilities of selecting split location
  double totXp, totTp; // sum of Xp, Tp
  
  // DLNMstruct info
  bool valid();
  double logPRule();
  NodeStruct* clone();
  bool checkEqual(NodeStruct*);
  void printStruct();
  int get(int);
  
  // proposal functions
  bool proposeSplit();
  NodeStruct* subStruct(bool);
  void dropSplit();
  void updateStruct(NodeStruct*, bool);
};

class ModStruct: public NodeStruct {
public:
  ModStruct(modDat*);
  ModStruct(modDat* modFncs, std::vector<std::vector<int> > availMod);
  ~ModStruct();
  ModStruct(const ModStruct&);

  int splitVar;   // splitting modifier
  int splitVal;   // splitting value (continuous modifier)
  std::vector<int> splitVec;
                  // splitting vector (categorical modifier)
  std::vector<std::vector<int> > availMod;
                  // index of remaining splits
  modDat *modFncs;// pointer to modifier functions (get avail modifiers)
    // varIsNum;  // index of continuous modifiers
    // nModSplit; // index of number of available splits
    // modProb;   // probability of selecting each modifier
    // totalProb(availMod);
                  // total probability of available modifiers
    // getAvailMod(splitVar, splitVal, splitVec, availMod, left);
                  // get new list of available modifiers for subStructs

  bool proposeSplit();
  void dropSplit();
  bool valid();
  bool checkEqual(NodeStruct*);
  void updateStruct(NodeStruct*, bool);

  NodeStruct* subStruct(bool);
  NodeStruct* clone();

  double logPRule();

  int get(int);
  std::vector<int> get2(int);
  std::vector<std::vector<int> > get3(int);
  void printStruct();

};
