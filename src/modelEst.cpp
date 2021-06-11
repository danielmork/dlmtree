/**
 * @file modelEst.cpp
 * @author Daniel Mork (danielmork.github.io)
 * @brief Functions for model estimation to accompany Bayesian treed distributed lag methods
 * @version 1.0
 * @date 2021-02-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include "RcppEigen.h"
#include "modelCtr.h"
#include "exposureDat.h"
#include "modDat.h"
#include "Fncs.h"
#include "Node.h"
#include "NodeStruct.h"
using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Lower;
using Eigen::Upper;

/**
 * @brief function to update fixed effect coef., sigma^2, polya-gamma
 * 
 * @param ctr model control data
 */
void tdlmModelEst(modelCtr *ctr)
{
  const VectorXd ZR = ctr->Zw.transpose() * ctr->R;
  ctr->gamma =        ctr->Vg * ZR;
   
  // * Update sigma^2 and xi_sigma2
  if (!(ctr->binomial)) {
    rHalfCauchyFC(&(ctr->sigma2), (double)ctr->n + (double)ctr->totTerm, 
                  ctr->R.dot(ctr->R) - ZR.dot(ctr->gamma) + 
                    ctr->sumTermT2 / ctr->nu, &(ctr->xiInvSigma2));
    if ((ctr->sigma2 != ctr->sigma2)) // ! stop if infinte or nan variance
      stop("\nNaN values (sigma) occured during model run, rerun model.\n");
  }
  
  // * Draw fixed effect coefficients
  ctr->gamma.noalias() += ctr->VgChol * 
    as<VectorXd>(rnorm(ctr->pZ, 0, sqrt(ctr->sigma2)));
    
  // * Update polya gamma vars
  if (ctr->binomial) {
    VectorXd psi =   ctr->fhat;
    psi.noalias() += ctr->Z * ctr->gamma;
    ctr->Omega =     rcpp_pgdraw(ctr->binomialSize, psi);
    ctr->Zw =        ctr->Omega.asDiagonal() * ctr->Z;
    MatrixXd VgInv(ctr->pZ, ctr->pZ);
    VgInv.triangularView<Lower>() =   ctr->Z.transpose() * ctr->Zw;
    VgInv.diagonal().array() +=       1 / 100000.0;
    VgInv.triangularView<Upper>() =   VgInv.transpose().eval();
    ctr->Vg.triangularView<Lower>() = VgInv.inverse();
    ctr->Vg.triangularView<Upper>() = ctr->Vg.transpose().eval();
    ctr->VgChol = ctr->Vg.llt().matrixL();
    ctr->Y =      ctr->kappa.array() / ctr->Omega.array();
    ctr->R =      ctr->Y - ctr->fhat; // Recalc R using new Y
  }
} // end tdlmModelEst function

/**
 * @brief Construct a new progress Meter::progress Meter object
 * 
 * @param c model control data to track meter progress
 */
progressMeter::progressMeter(modelCtr* c)
{
  ctr = c;
  startTime = time(NULL);
  if (ctr->verbose)
    Rcout << "Burn-in % complete \n" <<
      "[0--------25--------50--------75--------100]\n '";
  burnProgInc = (ctr->burn / 42.0);
  // burnProgMark = burnProgInc;
  iterProgInc = (ctr->iter / 42.0);
  // iterProgMark = double(ctr->burn) + iterProgInc;
  burnProgMark = 1.0;
  iterProgMark = 1.0;
}
progressMeter::~progressMeter()
{
  ctr = 0;
}
/**
 * @brief print progress mark "'" at given intervals thorughout model run 
 */
void progressMeter::printMark()
{
  if (ctr->verbose) {
    int completedMarks = 0.0;
    if (ctr->b > ctr->burn) {
      completedMarks = floor(42 * (ctr->b - ctr->burn) / ctr->iter);
      if (completedMarks > iterProgMark) {
        do {
          Rcout << "'";
          ++iterProgMark;
        } while (iterProgMark < completedMarks);
      }
    } else {
      completedMarks = floor(42 * ctr->b / ctr->burn);
      if (completedMarks > burnProgMark) {
        do {
          Rcout << "'";
          ++burnProgMark;
        } while (burnProgMark < completedMarks);
      }
      if (ctr->b == ctr->burn) {
        timediff = difftime(time(NULL), startTime);
        timediff = timediff * ctr->iter / ctr->burn;
        if (timediff > 3600) {
          Rprintf("\nMCMC iterations (est time: %.2g hours)\n",
                  round(100 * timediff / 3600) / 100);
        } else if (timediff > 60) {
          Rprintf("\nMCMC iterations (est time: %.2g minutes)\n",
                  round(100 * timediff / 60) / 100);
        } else {
          Rprintf("\nMCMC iterations (est time: %.2g seconds)\n",
                  round(100 * timediff) / 100);
        }
        Rcout << "[0--------25--------50--------75--------100]\n '";
      }
    } // end iter <= burn
  } // end if verbose
} // end progressMeter::printMark()


/**
 * @brief propose new treed DLM
 * 
 * @param tree pointer to current tree
 * @param Exp pointer to exposureDat containing exposure data
 * @param ctr pointer to model control
 * @param step grow (0), prune (1), change (2)
 * @return double MH ratio of current proposal
 */
double tdlmProposeTree(Node* tree, exposureDat* Exp, modelCtr* ctr, int step)
{
  int no = 0;
  double stepMhr = 0;
  std::vector<Node*> dlnmTerm, tempNodes;

  // List current tree terminal nodes
  dlnmTerm = tree->listTerminal();

  // Grow
  if (step == 0) {
    // select node to grow
    no = (std::size_t) floor(R::runif(0, dlnmTerm.size())); 

    if (dlnmTerm[no]->grow()) { // propose new split
      double nGen2 = double(tree->nGen2());
      if (dlnmTerm[no]->depth == 0) { // depth == 0
        ++nGen2;
      } else {
        if (!(dlnmTerm[no]->parent->isGen2())) {
          ++nGen2;
        }
      }
      
      stepMhr = log((double)tree->nTerminal()) - log(nGen2) +
        2 * logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth + 1, 1) +
        logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth, 0) -
        logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth, 1);

      Exp->updateNodeVals((dlnmTerm[no]->proposed)->c1); // update node values
      
      // newDlnmTerm = tree->listTerminal(1); // list proposed terminal nodes
    }


    // Prune
  } else if (step == 1) {
    tempNodes = tree->listGen2();
    no = floor(R::runif(0, tempNodes.size())); // select gen2 node to prune

    stepMhr = log((double)tree->nGen2()) - log((double)tree->nTerminal() - 1.0) -
      2 * logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], tempNodes[no]->depth + 1, 1) -
      logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], tempNodes[no]->depth, 0) +
      logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], tempNodes[no]->depth, 1);

    tempNodes[no]->prune(); // prune nodes
    // newDlnmTerm = tree->listTerminal(1); // list proposed terminal nodes


    // Change
  } else {
    tempNodes = tree->listInternal();
    no = floor(R::runif(0, tempNodes.size())); // select internal nodes to change 
    if (tempNodes[no]->change()) { // propose new split
      for (Node* tn : tempNodes[no]->proposed->listTerminal())
        Exp->updateNodeVals(tn);
        
      if ((tempNodes[no]->c1)->c1 != 0) { // calculate mhr if splits on c1 nodes
        for (Node* tn : tempNodes[no]->c1->listInternal())
          stepMhr -= (tn->nodestruct)->logPRule();
        for (Node* tn : tempNodes[no]->proposed->c1->listInternal())
          stepMhr += (tn->nodestruct)->logPRule();
      }

      if ((tempNodes[no]->c2)->c1 != 0) { // calculate mhr if splits on c2 nodes
        for (Node* tn : tempNodes[no]->c2->listInternal())
          stepMhr -= (tn->nodestruct)->logPRule();
        for (Node* tn : tempNodes[no]->proposed->c2->listInternal())
          stepMhr += (tn->nodestruct)->logPRule();
      }

    }
  }

  return(stepMhr);
} // end tdlmProposeTree function

/**
 * @brief propose new modifier tree
 * 
 * @param tree pointer to current tree
 * @param Mod pointer to modDat, containing modifier functions
 * @param ctr pointer to model control data
 * @param step grow (0), prune (1), change (2), or swap (3)
 * @return double MH ratio for current proposal
 */
double modProposeTree(Node* tree, modDat* Mod, dlmtreeCtr* ctr, int step)
{
  int no = 0;
  double stepMhr = 0;
  std::vector<Node*> modTerm, tempNodes;


  // List current tree terminal nodes
  modTerm = tree->listTerminal();
  
  // Grow
  if (step == 0) {
    // select node to grow
    no = (std::size_t) floor(R::runif(0, modTerm.size())); 

    if (modTerm[no]->grow()) { // propose new split
      double nGen2 = double(tree->nGen2());
      if (modTerm[no]->depth == 0) { // depth == 0
        ++nGen2;
      } else {
        if (!(modTerm[no]->parent->isGen2())) {
          ++nGen2;
        }
      }
      
      stepMhr = log((double)tree->nTerminal()) - log(nGen2) +
        2 * logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                      modTerm[no]->depth + 1, 1) +
        logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                  modTerm[no]->depth, 0) -
        logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                  modTerm[no]->depth, 1);

      Mod->updateNodeVals((modTerm[no]->proposed)->c1); // update node values
      if ((modTerm[no]->proposed)->c1->nodevals->idx.size() == 0 ||
          (modTerm[no]->proposed)->c2->nodevals->idx.size() == 0) {
        tree->reject();
        return(0);
      } // end reject if empty
    } // end grow proposal

  // Prune
  } else if (step == 1) {
    tempNodes = tree->listGen2();
    no = floor(R::runif(0, tempNodes.size())); // select gen2 node to prune

    stepMhr = log((double)tree->nGen2()) - log((double)tree->nTerminal() - 1.0) -
      2 * logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                    tempNodes[no]->depth + 1, 1) -
      logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                tempNodes[no]->depth, 0) +
      logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                tempNodes[no]->depth, 1);

    tempNodes[no]->prune(); // prune nodes
    if (tempNodes[no]->nodevals->nestedTree != 0)
      delete tempNodes[no]->nodevals->nestedTree;
    tempNodes[no]->nodevals->nestedTree = 0;

  // Change
  } else if (step == 2) {
    tempNodes = tree->listInternal();
    no = floor(R::runif(0, tempNodes.size())); // select internal nodes to change

    if (tempNodes[no]->change()) { // propose new split
      for (Node* tn : tempNodes[no]->proposed->listTerminal()) {
        Mod->updateNodeVals(tn);
        if (tn->nodevals->idx.size() == 0) {
          tree->reject();
          return(0);
        } // end reject if empty
      }

      if ((tempNodes[no]->c1)->c1 != 0) { // calculate mhr if splits on c1 nodes
        for (Node* tn : tempNodes[no]->c1->listInternal())
          stepMhr -= (tn->nodestruct)->logPRule();
        for (Node* tn : tempNodes[no]->proposed->c1->listInternal())
          stepMhr += (tn->nodestruct)->logPRule();
      }

      if ((tempNodes[no]->c2)->c1 != 0) { // calculate mhr if splits on c2 nodes
        for (Node* tn : tempNodes[no]->c2->listInternal())
          stepMhr -= (tn->nodestruct)->logPRule();
        for (Node* tn : tempNodes[no]->proposed->c2->listInternal())
          stepMhr += (tn->nodestruct)->logPRule();
      }
    }

  // Swap
  } else {
    tempNodes = tree->listInternal();
    no = floor(R::runif(0, tempNodes.size() - 1)); // dont select top node
    if (tempNodes[no]->parent->swap(tempNodes[no])) {
      for (Node* tn : tempNodes[no]->parent->proposed->listTerminal()) {
        Mod->updateNodeVals(tn);
        if (tn->nodevals->idx.size() == 0) {
          tree->reject();
          return(0);
        } // end reject if empty
      }

      for (Node* tn : tempNodes[no]->parent->listInternal())
        stepMhr -= (tn->nodestruct)->logPRule();
      for (Node* tn : tempNodes[no]->parent->proposed->listInternal())
        stepMhr += (tn->nodestruct)->logPRule();

    }
  }

  return(stepMhr);
} // end modProposeTree function

// void dlmtreeRecDLM(dlmtreeCtr* ctr, dlmtreeLog* dgn)
// {
//   double sumDLM;
//   for (int i = 0; i < ctr->n; ++i) {
//     dgn->exDLM.col(i) += ctr->exDLM.col(i);
//     dgn->ex2DLM.col(i) += ctr->exDLM.col(i).array().square().matrix();
//     sumDLM = ctr->exDLM.col(i).sum();
//     dgn->cumDLM(i) += sumDLM;
//     dgn->cum2DLM(i) += pow(sumDLM, 2);
//   }
// } // end dlmtreeRecDLM function

/**
 * @brief return a string describing current series of modifier rules leading to terminal node
 * 
 * @param n pointer to node
 * @param Mod pointer to modDat
 * @return std::string 
 */
std::string modRuleStr(Node* n, modDat* Mod)
{
  std::string rule = "";
  if (n->depth == 0)
    return(rule);
  
  Node* parent = n->parent;
  int splitVar = parent->nodestruct->get(1);
  int splitVal = parent->nodestruct->get(2);
  std::vector<int> splitVec = parent->nodestruct->get2(1);
  
  rule += std::to_string(splitVar);
  if (Mod->varIsNum[splitVar]) {
    if (parent->c1 == n)
      rule += "<";
    else
      rule += ">=";
    rule += std::to_string(splitVal);
  } else {
    if (parent->c1 == n)
      rule += "[]";
    else
      rule += "][";
    for (int i : splitVec)
      rule += std::to_string(i) + ",";
    rule.pop_back();
  }
  if (parent->depth != 0)
    rule += "&" + modRuleStr(parent, Mod);
  
  return(rule);
} // end modRuleStr function
                 
/**
 * @brief count of modifier variables used in current tree
 * 
 * @param tree pointer to tree
 * @param Mod pointer to modDat
 * @return VectorXd 
 */
VectorXd countMods(Node* tree, modDat* Mod)
{
  VectorXd modCount(Mod->nMods); modCount.setZero();
  VectorXd unavailProb(Mod->nMods); unavailProb.setZero();
  std::vector<int> unavail;  
  for (Node* tn : tree->listInternal()) {
    modCount(tn->nodestruct->get(1)) += 1.0;
    unavail.clear();
    unavailProb.setZero();
    for (int i = 0; i < Mod->nMods; ++i) {
      if (tn->nodestruct->get3(1)[i].size() == 0) {
        unavail.push_back(i);
        unavailProb(i) = Mod->modProb[i];
      }
    }
    if (unavail.size() > 0) {
      std::random_shuffle(unavail.begin(), unavail.end());
      double totProb = unavailProb.sum();
      int pseudoDraw = R::rgeom(std::max(0.00000001, 1 - totProb));
      int binomDraw = 0;
      if (pseudoDraw > 0) {
        for (int i : unavail) {
          binomDraw = R::rbinom(pseudoDraw, unavailProb(i) / totProb);
          if (binomDraw > 0)
            modCount(i) += binomDraw * 1.0;
          totProb -= unavailProb(i);
          pseudoDraw -= binomDraw;
          if (pseudoDraw < 1)
            break;
        } // end multinom
      } // end pseudoDraw
    } // end unavail
  } // end modCount
  return(modCount);
} // end countMods function

/**
 * @brief draw tree structure from tree prior distribution
 * 
 * @param tree pointer to tree being grown
 * @param n pointer to specific node in tree
 * @param alpha alpha parameter (0 to 1)
 * @param beta beta parameter (>0)
 */
void drawTree(Node* tree, Node* n, double alpha, double beta)
{
  double logProb = log(alpha) - beta * log(1.0 + n->depth);
  if (log(R::runif(0, 1)) < logProb) {
    if (n->grow()) {
      if (n->depth > 0)
        n = n->proposed;
      tree->accept();
      drawTree(tree, n->c1, alpha, beta);
      drawTree(tree, n->c2, alpha, beta);
    }
  } // end grow tree
  return;
} // end drawTree function

/**
 * @brief update design matrices for subgroup Gaussian process DLM
 * 
 * @param n pointer to node
 * @param ctr pointer to model control
 */
void updateGPMats(Node* n, dlmtreeCtr* ctr)
{
  if (n->nodevals->updateXmat == 0)
    return;
  if (n->depth == 0) {
    n->nodevals->XtX = ctr->XtXall;
    n->nodevals->ZtXmat = ctr->ZtXall;
    n->nodevals->VgZtXmat = ctr->VgZtXall;
    n->nodevals->updateXmat = 0;
    return;
  }
  Node* par = n->parent;
  if (par->nodevals->updateXmat)
    updateGPMats(par, ctr);
  Node* sib = n->sib();
  std::vector<int> idx;
  if (n->nodevals->idx.size() <= sib->nodevals->idx.size()) {
    idx = n->nodevals->idx;
  } else {
    idx = sib->nodevals->idx;
  }
  
  MatrixXd Xtemp(idx.size(), ctr->pX); Xtemp.setZero();
  MatrixXd Ztemp(idx.size(), ctr->pZ); Ztemp.setZero();
  
  for (std::size_t i = 0; i < idx.size(); ++i) {
    Xtemp.row(i) = ctr->X.row(idx[i]);
    Ztemp.row(i) = ctr->Z.row(idx[i]);
  }
  
  if (n->nodevals->idx.size() <= sib->nodevals->idx.size()) {
    n->nodevals->XtX = Xtemp.transpose() * Xtemp;
    n->nodevals->ZtXmat = Ztemp.transpose() * Xtemp;
    n->nodevals->VgZtXmat = ctr->Vg * n->nodevals->ZtXmat;
    sib->nodevals->XtX = par->nodevals->XtX - n->nodevals->XtX;
    sib->nodevals->ZtXmat = par->nodevals->ZtXmat - n->nodevals->ZtXmat;
    sib->nodevals->VgZtXmat = par->nodevals->VgZtXmat - n->nodevals->VgZtXmat;
  } else {
    sib->nodevals->XtX = Xtemp.transpose() * Xtemp;
    sib->nodevals->ZtXmat = Ztemp.transpose() * Xtemp;
    sib->nodevals->VgZtXmat = ctr->Vg * sib->nodevals->ZtXmat;
    n->nodevals->XtX = par->nodevals->XtX - sib->nodevals->XtX;
    n->nodevals->ZtXmat = par->nodevals->ZtXmat - sib->nodevals->ZtXmat;
    n->nodevals->VgZtXmat = par->nodevals->VgZtXmat - sib->nodevals->VgZtXmat;
  }
  n->nodevals->updateXmat = 0;
  sib->nodevals->updateXmat = 0;
  
}