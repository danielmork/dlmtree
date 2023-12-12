#define NODE_H
#include <RcppEigen.h>
class NodeStruct;
class Node;

struct NodeVals {
  NodeVals(int n, int pZ); 
  NodeVals(int n);
  ~NodeVals();
  NodeVals(const NodeVals&);
  
public:
  Eigen::VectorXd X;
  Eigen::MatrixXd XtX;
  Eigen::VectorXd ZtX;
  Eigen::MatrixXd ZtXmat;
  Eigen::VectorXd VgZtX;
  Eigen::MatrixXd VgZtXmat;
  Eigen::MatrixXd tempV;
  bool updateXmat;
  std::vector<int> idx;
  Node* nestedTree;
};

class Node {
public:
  Node(int depth_in = 0, bool update_in = 1);
  ~Node();
  Node(const Node&);

  // node variables
  int depth;
  bool update;
  Node *c1, *c2, *parent, *proposed;
  NodeStruct *nodestruct;
  NodeVals *nodevals;

  // tree proposal functions
  bool grow();
  void prune();
  bool change();
  bool swap(Node* c);
  void accept();
  void reject();

  // tree information functions
  int nTerminal();
  int nGen2();
  std::vector<Node*> listTerminal(bool follow_proposed = 0);
  std::vector<Node*> listInternal();
  std::vector<Node*> listGen2();
  // double logPTree(double alpha, double beta);
  
  // current node information functions
  Node* sib();
  bool isProposed();
  bool isGen2();

  // nodevals and nodestruct functions
  void replaceTree(Node* newTree);
  void replaceNodeVals(Node* n);
  void setUpdate(bool update);
  void setUpdateXmat(bool update);
  bool updateStruct();
};

std::vector<Node*> CombineNodeLists(std::vector<Node*> x, std::vector<Node*> y);

