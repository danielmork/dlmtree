#include <Rcpp.h>
#include "Node.h"
#include "NodeStruct.h"
using namespace std;
using namespace Rcpp;


NodeVals::NodeVals(int n, int pZ)
{
  X.resize(n);
  ZtX.resize(pZ);
  VgZtX.resize(pZ);
  tempV.resize(0, 0);
  updateXmat = 1;
  nestedTree = 0;
}
NodeVals::NodeVals(int n)
{
  updateXmat = 1;
  nestedTree = 0;
}
NodeVals::~NodeVals()
{
  X.resize(0);
  ZtX.resize(0);
  VgZtX.resize(0);
  XtX.resize(0, 0);
  ZtXmat.resize(0, 0);
  VgZtXmat.resize(0, 0);
  tempV.resize(0, 0);
  if (nestedTree != 0) {
    delete nestedTree;
    nestedTree = 0;
  }
}
NodeVals::NodeVals(const NodeVals& x)
{
  X = x.X;
  Xpl = x.Xpl;
  XplProposed = x.XplProposed;
  XtX = x.XtX;
  ZtX = x.ZtX;
  ZtXmat = x.ZtXmat;
  ZtXmatProposed = x.ZtXmatProposed;
  VgZtX = x.VgZtX;
  VgZtXmat = x.VgZtXmat;
  VgZtXmatProposed = x.VgZtXmatProposed;
  tempV = x.tempV;
  idx = x.idx;
  updateXmat = x.updateXmat;
  nestedTree = 0;
  if (x.nestedTree != 0) {
    nestedTree = new Node(*(x.nestedTree));
  }
}

Node::Node(int depth_in,
           bool update_in)
{
  depth = depth_in;
  update = update_in;
  c1 = 0; c2 = 0; parent = 0;
  proposed = 0;
  nodevals = 0;
  nodestruct = 0;
}

Node::~Node()
{
  if (nodevals != 0) {
    delete nodevals;
    nodevals = 0;
  }
  if (nodestruct != 0) {
    delete nodestruct;
    nodestruct = 0;
  }

  if (c1 != 0) {
    delete c1;
    c1 = 0;
  }

  if (c2 != 0) {
    delete c2;
    c2 = 0;
  }

  if (proposed != 0) {
    delete proposed;
    proposed = 0;
  }
  
  parent = 0;
}

Node::Node(const Node& n)
{
  depth = n.depth;
  update = n.update;
  c1 = 0; c2 = 0; parent = 0;
  proposed = 0;
  nodevals = 0;
  nodestruct = 0;
  nodestruct = (n.nodestruct)->clone();
  if (n.nodevals != 0) {
    nodevals = new NodeVals(*(n.nodevals));
  }
  if (n.c1 != 0) {
    c1 = new Node(*(n.c1));
    c2 = new Node(*(n.c2));
    c1->parent = this;
    c2->parent = this;
  }
}


bool Node::grow()
{
  NodeStruct* new_ns = nodestruct->clone();

  if (new_ns->proposeSplit()) {
    proposed = new Node(depth, update);
    proposed->nodestruct = new_ns;

    proposed->c1 = new Node(depth + 1, 1);
    proposed->c2 = new Node(depth + 1, 1);
    (proposed->c1)->nodestruct = new_ns->subStruct(1);
    (proposed->c2)->nodestruct = new_ns->subStruct(0);

    proposed->nodevals = nodevals;
    (proposed->c1)->parent = proposed;
    (proposed->c2)->parent = proposed;

    return(1);
  }

  delete new_ns;
  new_ns = 0;
  return(0);
}

void Node::prune()
{
  if (c1) {
    proposed = new Node(depth, update);
    proposed->nodestruct = nodestruct->clone();
    proposed->nodevals = nodevals;
    (proposed->nodestruct)->dropSplit();
  }
}

bool Node::change()
{
  if (c1 == 0) {
    return(0);
  }

  NodeStruct* new_ns = nodestruct->clone();

  if (new_ns->proposeSplit()) {
    proposed = new Node(depth, update);
    proposed->nodestruct = new_ns;
    proposed->c1 = new Node(*c1);
    proposed->c2 = new Node(*c2);
    if (!(proposed)->updateStruct()) {
      delete proposed;
      proposed = 0;
      return(0);
    }

    proposed->nodevals = nodevals;
    (proposed->c1)->parent = proposed;
    (proposed->c2)->parent = proposed;
    (proposed->c1)->setUpdate(1);
    (proposed->c2)->setUpdate(1);


    return(1);
  }

  delete new_ns;
  proposed = 0;
  return(0);
}

// Swap node rule with child node rule
bool Node::swap(Node* child)
{
  if (child->c1 == 0) {
    return(0);
  }

  // Create proposed nodes
  proposed = new Node(depth, update);
  proposed->nodestruct = child->nodestruct->clone();
  proposed->c1 = new Node(*c1);
  proposed->c2 = new Node(*c2);
  if (c1 == child) {
    proposed->c1->nodestruct = nodestruct->clone();
    if (c2->c1 != 0) {
      if (c2->nodestruct->checkEqual(c1->nodestruct))
        proposed->c2->nodestruct = nodestruct->clone();
    }
  } else {
    proposed->c2->nodestruct = nodestruct->clone();
    if (c1->c1 != 0) {
      if (c1->nodestruct->checkEqual(c2->nodestruct))
        proposed->c1->nodestruct = nodestruct->clone();
    }
  }

  // Update structures of proposed nodes, check for issues
  if (proposed->updateStruct()) {
    proposed->nodevals = nodevals;
    (proposed->c1)->parent = proposed;
    (proposed->c2)->parent = proposed;
    (proposed->c1)->setUpdate(1);
    (proposed->c2)->setUpdate(1);
    return(1);
  }
  delete proposed;
  proposed = 0;
  return(0);
}

void Node::accept()
{
  // If update occured at top of tree, change out structure and subnodes
  if ((depth == 0) && (proposed != 0)) {

    delete nodestruct;
    nodestruct = proposed->nodestruct;
    proposed->nodestruct = 0;
    update = proposed->update;

    if (c1 != 0) {
      delete c1;
      delete c2;
      c1 = 0;
      c2 = 0;
    }

    if ((proposed->c1) != 0) {
      c1 = proposed->c1;
      c2 = proposed->c2;
      c1->parent = this;
      c2->parent = this;
    }

    proposed->c1 = 0;
    proposed->c2 = 0;
    if (nodevals != proposed->nodevals)
      delete proposed->nodevals;
    proposed->nodevals = 0;

    delete proposed;
    proposed = 0;

  // Depth is not zero, check if child nodes have new proposals
  } else {
    if (c1 != 0) {
      if (c1->proposed != 0) {

        Node* temp = c1;
        c1 = (c1->proposed);
        temp->proposed = 0;
        if (temp->nodevals != c1->nodevals)
          delete temp->nodevals;
        temp->nodevals = 0;
        delete temp;
        temp = 0;
        c1->parent = this;

      } else {

        c1->accept();

      }
    }

    if (c2 != 0) {
      if (c2->proposed != 0) {

        Node* temp = c2;
        c2 = (c2->proposed);
        temp->proposed = 0;
        if (temp->nodevals != c2->nodevals)
          delete temp->nodevals;
        temp->nodevals = 0;
        delete temp;
        temp = 0;
        c2->parent = this;

      } else {

        c2->accept();

      }
    }
  }
}

void Node::reject()
{
  if (proposed != 0) {
    if ((proposed->nodevals != 0)) {
      if (proposed->nodevals != nodevals)
        delete proposed->nodevals;
    }
    proposed->nodevals = 0;
      
    if (proposed->c1 != 0) {
      delete proposed->c1;
      delete proposed->c2;
    }
    proposed->c1 = 0;
    proposed->c2 = 0;
    
    if (proposed->nodestruct != 0) {
      delete proposed->nodestruct;
    }
    proposed->nodestruct = 0;
    
    delete proposed;
    proposed = 0;

    return;
  } else {

    if (c1 != 0) {
      c1->reject();
      c2->reject();
    }
    return;
  }
}




int Node::nTerminal()
{
  if (c1 == 0) {
    return(1);

  } else {
    return(c1->nTerminal() + c2->nTerminal());
  }
}

int Node::nGen2()
{
  if (c1 == 0) {
    return(0);

  } else if ((c1->c1 == 0) && (c2->c2 == 0)) {
    return(1);

  } else {
    return(c1->nGen2() + c2->nGen2());
  }
}



std::vector<Node*> Node::listTerminal(bool follow_proposed)
{
  Node* follow;
  if ((follow_proposed) && (proposed != 0)) {
    follow = proposed;
  } else {
    follow = this;
  }
  if (follow->c1 == 0) {
    std::vector<Node*> v = {follow};
    return(v);

  } else {
    return(CombineNodeLists((follow->c1)->listTerminal(follow_proposed),
                            (follow->c2)->listTerminal(follow_proposed)));
  }
}

std::vector<Node*> Node::listInternal()
{
  std::vector<Node*> v;
  if (c1 != 0) {
    v = CombineNodeLists(c1->listInternal(), c2->listInternal());
    v.push_back(this);
  }
  return(v);
}

std::vector<Node*> Node::listGen2()
{
  if (c1 == 0) {
    std::vector<Node*> v;
    return(v);

  } else if ((c1->c1 == 0) && (c2->c2 == 0)) {
    std::vector<Node*> v = {this};
    return(v);

  } else {
    return(CombineNodeLists(c1->listGen2(), c2->listGen2()));
  }
}



Node* Node::sib()
{
  if (depth > 0) {

    if (this == (parent->c1)) {
      if ((parent->c2)->proposed == 0) {
        return(parent->c2);
      } else {
        return((parent->c2)->proposed);
      }
    }

    if (this == parent->c2) {
      if ((parent->c1)->proposed == 0) {
        return(parent->c1);
      } else {
        return((parent->c1)->proposed);
      }
    }

    if (parent->proposed != 0) {
      if (this == (parent->proposed)->c1) {
        return((parent->proposed)->c2);
      } else if (this == (parent->proposed)->c2) {
        return((parent->proposed)->c1);
      }
    }


    if (((parent->c1)->proposed != 0) |
        ((parent->c2)->proposed != 0)) {
      if (this == (parent->c1)->proposed) {
        return(parent->c2);
      }
      if (this == (parent->c2)->proposed) {
        return(parent->c1);
      }
    }
  }
  return(0);
}







bool Node::isProposed()
{
  if (proposed != 0) {
    return(1);
  } else if (c1 != 0) {
    if ((c1->isProposed()) || (c2->isProposed())) {
      return(1);
    }
  }
  return(0);
}

bool Node::isGen2()
{
  if (c1 != 0) {
    if ((c1->c1 == 0) && (c2->c1 == 0)) {
      return(1);
    } else {
      return(0);
    }
  }
  return(0);
}


void Node::replaceTree(Node* newTree)
{
  // Delete previous attributes
  if (nodevals != 0) {
      delete nodevals;
      nodevals = 0;
    }
  if (nodestruct != 0) {
    delete nodestruct;
    nodestruct = 0;
  }

  if (c1 != 0) {
    delete c1;
    c1 = 0;
  }

  if (c2 != 0) {
    delete c2;
    c2 = 0;
  }

  if (proposed != 0) {
    delete proposed;
    proposed = 0;
  }

  parent = 0;

  // Update to new information
  depth = newTree->depth;
  update = newTree->update;
  nodestruct = (newTree->nodestruct)->clone();
  if (newTree->nodevals != 0) {
    nodevals = new NodeVals(*(newTree->nodevals));
  }
  if (newTree->c1 != 0) {
    c1 = new Node(*(newTree->c1));
    c2 = new Node(*(newTree->c2));
    c1->parent = this;
    c2->parent = this;
  }
}



void Node::replaceNodeVals(Node* newTree)
{
  if (nodevals != 0) {
    delete nodevals;
    nodevals = 0;
  }
  nodevals = newTree->nodevals;
  newTree->nodevals = 0;
  if ((c1 != 0) && (newTree->c1 != 0)) {
    c1->replaceNodeVals(newTree->c1);
    c2->replaceNodeVals(newTree->c2);
  }
}

void Node::setUpdate(bool x) {
  update = x;
  if (c1 != 0) {
    c1->setUpdate(x);
    c2->setUpdate(x);
  }
}

void Node::setUpdateXmat(bool x) {
  if (nodevals != 0)
    nodevals->updateXmat = x;
  if (c1 != 0) {
    c1->setUpdateXmat(x);
    c2->setUpdateXmat(x);
  }
}

bool Node::updateStruct() {
  if (c1 == 0) {
    return(1);

  } else {
    (c1->nodestruct)->updateStruct(this->nodestruct, 1);

    if ((c1->nodestruct)->valid()) {
      (c2->nodestruct)->updateStruct(this->nodestruct, 0);

      if ((c2->nodestruct)->valid()) {
        c1->update = 1;
        c2->update = 1;

        if (c1->updateStruct()) {
          if (c2->updateStruct()) {
            return(1);

          }
        }
      }
    }
  }
  return(0);
}


std::vector<Node*> CombineNodeLists(std::vector<Node*> x, std::vector<Node*> y)
{
  if (x.size() >= y.size()) {
    if (!y.empty()) {
      for (size_t i = 0; i < y.size(); ++i) {
        x.push_back(y[i]);
      }
    }
    return(x);

  } else if (y.size() > x.size()) {
    if (!x.empty()) {
      for (size_t i = 0; i < x.size(); ++i) {
        y.push_back(x[i]);
      }
    }
    return(y);
  }

  return(x);
}



