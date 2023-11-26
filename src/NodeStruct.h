#define NODESTRUCT_H
#include <RcppEigen.h>
#include <Rcpp.h>

class modDat;

class NodeStruct {
public: //virtual = override parent functions
  NodeStruct();
  virtual ~NodeStruct(); // Destructor
  NodeStruct(const NodeStruct&);
  virtual NodeStruct* clone();
  virtual NodeStruct* subStruct(bool);
  virtual bool valid();
  virtual double logPRule();
  virtual bool proposeSplit();
  virtual void dropSplit();
  virtual void printStruct();
  virtual int get(int);
  virtual Eigen::VectorXd getTimeProbs();
  virtual void setTimeProbs(Eigen::VectorXd);
  virtual std::vector<int> get2(int);
  virtual std::vector<std::vector<int> > get3(int);
  virtual void updateStruct(NodeStruct*, bool);
  virtual bool checkEqual(NodeStruct*);
  virtual void setTimeRange(int, int);
};

class DLNMStruct: public NodeStruct {
public:
  DLNMStruct(int xmin_in, int xmax_in, int tmin_in, int tmax_in,
                       Eigen::VectorXd Xp_in, Eigen::VectorXd Tp_in);
  ~DLNMStruct();
  DLNMStruct(const DLNMStruct& ns);

  int xmin, xmax, tmin, tmax; // DLNM limits
  int xsplit, tsplit; // set to 0 until split is created
  Eigen::VectorXd Xp, Tp; // probabilities of selecting split location
  double totXp, totTp; // sum of Xp, Tp
  
  // DLNMstruct info
  bool valid();
  double logPRule();
  NodeStruct* clone();
  bool checkEqual(NodeStruct* n);
  void printStruct();
  int get(int a);
  Eigen::VectorXd getTimeProbs();
  
  // proposal functions
  bool proposeSplit();
  NodeStruct* subStruct(bool left);
  void dropSplit();
  void updateStruct(NodeStruct* parStruct, bool left);
  void setTimeRange(int lower, int upper);
  void setTimeProbs(Eigen::VectorXd newProbs);
};

class ModStruct: public NodeStruct {
public:
  ModStruct(modDat* md);
  ModStruct(modDat* modFncs, std::vector<std::vector<int> > availMod);
  ~ModStruct();
  ModStruct(const ModStruct& ns);

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
  bool checkEqual(NodeStruct* ns);
  void updateStruct(NodeStruct* parStruct, bool left);

  NodeStruct* subStruct(bool left);
  NodeStruct* clone();

  double logPRule();

  int get(int a);
  std::vector<int> get2(int a);
  std::vector<std::vector<int> > get3(int a);
  void printStruct();

};
