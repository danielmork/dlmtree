#include <RcppEigen.h>
#include <Rcpp.h>


class NodeStruct {
  // int xmin, xmax, tmin, tmax;
public:
  // int xsplit, tsplit;
  NodeStruct();
  virtual ~NodeStruct();
  NodeStruct(const NodeStruct&);
  virtual NodeStruct* clone();
  virtual NodeStruct* subStruct(bool);
  // virtual NodeStruct* copy();
  virtual bool valid();
  virtual double logPRule();
  virtual bool proposeSplit();
  virtual void dropSplit();
  virtual void printStruct();
  virtual int get(int);
  // virtual void setSplit(NodeStruct*);
  virtual void updateStruct(NodeStruct*, bool);
};

class DLNMStruct: public NodeStruct {
public:
  DLNMStruct(int, int, int, int, Eigen::VectorXd, Eigen::VectorXd);
  ~DLNMStruct();
  DLNMStruct(const DLNMStruct&);

  int xmin, xmax, tmin, tmax;
  double totXp, totTp;
  Eigen::VectorXd Xp, Tp;
  // Eigen::VectorXi aXs, aTs;
  // Eigen::VectorXd aXp, aTp;
  int xsplit, tsplit; // set to 0 until split is created
  bool valid();
  double logPRule();
  bool proposeSplit();
  NodeStruct* subStruct(bool);
  NodeStruct* clone();
  void dropSplit();
  void printStruct();
  // void setSplit(NodeStruct*);
  // NodeStruct* copy();
  int get(int);
  void updateStruct(NodeStruct*, bool);

};
