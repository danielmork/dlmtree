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
  
  Eigen::MatrixXd XplProposed;
  Eigen::MatrixXd ZtXmatProposed;
  Eigen::MatrixXd VgZtXmatProposed;
  Eigen::MatrixXd Xpl;
};

class Node {
public:
  Node(int depth_in = 0, bool update_in = 1);
  ~Node();
  Node(const Node&);

  // node variables
  int depth; // depth of node beginning at 0
  bool update; // boolean to indicate update needed following proposal
  Node *c1, *c2, *parent, *proposed; // pointers to child1 node, child2 node, parent node, proposed replacement
  NodeStruct *nodestruct; // pointer to node structure
  NodeVals *nodevals; // pointer to class NodeVals containig information for current node

  // tree proposal functions
  bool grow(); // Add a rule and two proposed child nodes at current node
  void prune(); // function to remove child nodes and make current node terminal
  bool change(); // functino to propose new splitting rule at current node
  bool swap(Node* c); // function to swap rule with input child node
  void accept(); // function to accept proposal (replacement)
  void reject(); // function to reject proposal (replacement)

  // tree information functions
  int nTerminal(); // Count terminal nodes of current tree
  int nGen2(); // count generation 2 internal nodes (nodes with two child terminal nodes)
  std::vector<Node*> listTerminal(bool follow_proposed = 0); // list all terminal nodes of current (0) or proposed (1)=
  std::vector<Node*> listInternal(); // list internal nodes
  std::vector<Node*> listGen2(); // list gen 2 nodes
  // double logPTree(double alpha, double beta);
  
  // current node information functions
  Node* sib(); // finds a sibling node
  bool isProposed(); // Check if the node is a proposed node
  bool isGen2(); // Check if the node is a Gen2 node

  // nodevals and nodestruct functions
  void replaceTree(Node* newTree);
  void replaceNodeVals(Node* n); // Replace current node values and all child node values with those of node in input pointer
  void setUpdate(bool update); // Set update value for current and child nodes
  void setUpdateXmat(bool update); // set NodeVals updateXmat value for current and child nodes
  bool updateStruct(); // Update NodeStruct after change in tree structure
  // Updates possible splits, splitting probability, and check that splits are valid.
};

std::vector<Node*> CombineNodeLists(std::vector<Node*> x, std::vector<Node*> y);

