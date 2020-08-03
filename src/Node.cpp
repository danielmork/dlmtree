#include <Rcpp.h>
#include "Node.h"
using namespace std;
using namespace Rcpp;


NodeVals::NodeVals(int n, int pZ)
{
  Eigen::VectorXd X(n);
  Eigen::VectorXd ZtX(pZ);
  Eigen::VectorXd VgZtX(pZ);
  Eigen::MatrixXd tempV(0, 0);
}
NodeVals::~NodeVals()
{
  X.resize(0);
  ZtX.resize(0);
  VgZtX.resize(0);
  tempV.resize(0, 0);
}
NodeVals::NodeVals(const NodeVals& x)
{
  X = x.X;
  ZtX = x.ZtX;
  VgZtX = x.VgZtX;
  tempV = x.tempV;
}

Node::Node(int depth_in = 0,
           bool update_in = 1)
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

  if (c2 != 0) {
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
}

Node::Node(const Node& n)
{
  // Rcout << "\nCopy Node\n";
  depth = n.depth;
  update = n.update;
  c1 = 0; c2 = 0; parent = 0;
  proposed = 0;
  nodevals = 0;
  nodestruct = 0;
  nodestruct = (n.nodestruct)->clone();
  if (n.nodevals != 0) {
    // Rcout << "\nCopy NodeVals\n";
    nodevals = new NodeVals(*(n.nodevals));
  }
  if (n.c1 != 0) {
    // Rcout << "\nCopy c1 c2\n";
    c1 = new Node(*(n.c1));
    c2 = new Node(*(n.c2));
    // Rcout << "\nSet Parent\n";
    c1->parent = this;
    c2->parent = this;
  }
  // if (depth == 0)
    // Rcout << "\nCopy complete\n";
}
// Node* Node::clone() { return new Node(*this); }

// Node* Node::copy(Node* copyParent, bool left, bool set_update)
// {
//   NodeStruct *ns;
//   if ((depth == 0) || (copyParent == 0)) {
//     // Rcout << "nopar\n";
//     ns = nodestruct->clone();
//   } else {
//     // Rcout << "par\n";
//     if (left) {
//       ns = (copyParent->nodestruct)->subStruct(nodestruct, 1);
//     } else {
//       ns = (copyParent->nodestruct)->subStruct(nodestruct, 0);
//     }
//   }
//
//   // Rcout << "\norig\n";
//   // nodestruct->printStruct();
//   // Rcout << "\ncopy\n";
//   // ns->printStruct();
//   if (!ns->valid()) {
//     delete ns;
//     ns = 0;
//     return(0);
//   }
//
//   Node *n = new Node(depth, set_update);
//   n->nodestruct = ns;
//   // if (set_update == 0)
//   //   n->nodevals = nodevals;
//   if (copyParent != 0)
//     n->parent = copyParent;
//
//   if (c1 != 0) {
//     // Rcout << "copych\n" << c1;
//     n->c1 = c1->copy(n, 1, set_update);
//     if (n->c1 != 0) {
//       n->c2 = c2->copy(n, 0, set_update);
//       if (n->c2 == 0) {
//         delete n->c1;
//         n->c1 = 0;
//         delete n;
//         n = 0;
//         return(0);
//       }
//     } else {
//       delete n;
//       n = 0;
//       return(0);
//     }
//     // c1->parent = n;
//     // c2->parent = n;
//   }
//   return(n);
// }

bool Node::grow()
{
  NodeStruct* new_ns = nodestruct->clone();
  // Rcout << depth << '\n';
  if (new_ns->proposeSplit()) {
    // Rcout << new_ns->get(1) << new_ns->get(2) <<
    //   new_ns->get(3) << new_ns->get(4) <<
    //     new_ns->get(5) << new_ns->get(6);
    proposed = new Node(depth, update);
    proposed->nodestruct = new_ns;

    proposed->c1 = new Node(depth + 1, 1);
    proposed->c2 = new Node(depth + 1, 1);
    (proposed->c1)->nodestruct = new_ns->subStruct(1);
    (proposed->c2)->nodestruct = new_ns->subStruct(0);

    // if (depth == 0) {
    proposed->nodevals = nodevals;
    (proposed->c1)->parent = proposed;
    (proposed->c2)->parent = proposed;
    // } else {
    //   proposed->parent = parent;
    //   (proposed->c1)->parent = proposed;
    //   (proposed->c2)->parent = proposed;
    // }

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
    // nodestruct->printStruct();
    // c1->nodestruct->printStruct();
    // c2->nodestruct->printStruct();
    // new_ns->printStruct();
    proposed = new Node(depth, update);
    proposed->nodestruct = new_ns;
    proposed->c1 = new Node(*c1);
    proposed->c2 = new Node(*c2);
    if (!(proposed)->updateStruct()) {
      delete proposed;
      proposed = 0;
      return(0);
    }
    // proposed->c1->nodestruct->printStruct();
    // proposed->c2->nodestruct->printStruct();

    proposed->nodevals = nodevals;
    (proposed->c1)->parent = proposed;
    (proposed->c2)->parent = proposed;
    (proposed->c1)->setUpdate(1);
    (proposed->c2)->setUpdate(1);


    return(1);
  }

  // Rcout << "propose";
  delete new_ns;
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

    if (proposed->nodevals != nodevals)
      delete proposed->nodevals;
    proposed->nodevals = 0;
    delete proposed;
    proposed = 0;

  } else {

    if (c1 != 0) {
      c1->reject();
      c2->reject();
    }

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
  // Rcout << "d" << depth << "-" << update << " ";
  if (c1 != 0) {
    c1->setUpdate(x);
    c2->setUpdate(x);
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



