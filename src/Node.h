#include <RcppEigen.h>
#include "NodeStruct.h"

struct NodeVals {
  NodeVals(int n, int pZ = 0);
  ~NodeVals();
  NodeVals(const NodeVals&);
public:
  Eigen::VectorXd X;
  Eigen::VectorXd ZtX;
  Eigen::VectorXd VgZtX;
  Eigen::MatrixXd tempV;
  // Eigen::VectorXd idx;
};

class Node {
public:
  Node(int, bool);
  ~Node();
  Node(const Node&);

  // Node variables
  int depth;
  bool update;
  Node *c1, *c2, *parent, *proposed;
  NodeStruct *nodestruct;
  NodeVals *nodevals;

  // Node functions
  // Node* clone();
  bool grow();
  void prune();
  bool change();
  // Node* copy(Node*, bool, bool);
  void accept();
  void reject();
  int nTerminal();
  int nGen2();
  std::vector<Node*> listTerminal(bool follow_proposed = 0);
  std::vector<Node*> listInternal();
  std::vector<Node*> listGen2();
  bool isProposed();
  bool isGen2();
  Node* sib();
  void replaceNodeVals(Node*);
  void setUpdate(bool);
  bool updateStruct();
};

std::vector<Node*> CombineNodeLists(std::vector<Node*> x, std::vector<Node*> y);

