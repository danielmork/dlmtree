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
#include <random>
#include <iostream>
#include <algorithm>
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
  if(!(ctr->zinb)){ 
    const VectorXd ZR = ctr->Zw.transpose() * ctr->R; 
    ctr->gamma = ctr->Vg * ZR; 
    
    // * Update sigma^2 and xi_sigma2
    if (!(ctr->binomial)) {
      rHalfCauchyFC(&(ctr->sigma2), (double)ctr->n + (double)ctr->totTerm, 
                    ctr->R.dot(ctr->R) - ZR.dot(ctr->gamma) + ctr->sumTermT2 / ctr->nu, &(ctr->xiInvSigma2));

      if ((ctr->sigma2 != ctr->sigma2)) {// ! stop if infinite or nan variance
        // Rcout << ctr->sigma2 << " " << ctr->totTerm << " " << 
        //   ctr->R.dot(ctr->R) << " " << ZR.dot(ctr->gamma) << " " << ctr->R << " " <<
        //   " " << ctr->Z << " " << ZR << " " << (ctr->gamma) << " " << ctr->sumTermT2 / ctr->nu << " " << ctr->xiInvSigma2;
        stop("\nNaN values (sigma) occured during model run, rerun model.\n");
      }
    }

    // * Draw fixed effect coefficients' variance
    ctr->gamma.noalias() += ctr->VgChol * as<VectorXd>(rnorm(ctr->pZ, 0, sqrt(ctr->sigma2))); 

    // * Update polya gamma vars
    if (ctr->binomial) {
      VectorXd psi =   ctr->fhat + ctr->Z * ctr->gamma;
      
      // Latent variable, Omega
      ctr->Omega = rcpp_pgdraw(ctr->binomialSize, psi); 
      ctr->Zw = ctr->Omega.asDiagonal() * ctr->Z;

      // Constructing V_gamma Inverse 
      Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ); 
      VgInv.triangularView<Eigen::Lower>() = ctr->Z.transpose() * ctr->Zw;
      VgInv.diagonal().array() += 1 / 100000.0; 
      VgInv.triangularView<Eigen::Upper>() = VgInv.transpose().eval(); 

      // Constructing V_gamma = Inverse of V_gamma Inverse
      ctr->Vg.triangularView<Eigen::Lower>() = VgInv.inverse();
      ctr->Vg.triangularView<Eigen::Upper>() = ctr->Vg.transpose().eval();
      
      // Update the V_gamma cholesky using LLT Decomposition, Lower triangular part of matrix L
      ctr->VgChol = ctr->Vg.llt().matrixL();
      // recalculate 'pseudo-Y' = kappa / omega, kappa = (y - n_b)/2
      ctr->Y =      ctr->kappa.array() / ctr->Omega.array();
      ctr->R = ctr->Y - ctr->fhat; // Recalc R using new Y
    }


  } else { // ZINB
    // *** Step 1: Update gamma_1 (ZI coefficient) ***
    // 1-1: Calculate eta1 and logit1
    Eigen::VectorXd eta1 = (ctr->Z1 * ctr->b1);

    // 1-2: Sample PG(1, eta1)
    ctr->omega1 = rcpp_pgdraw(ctr->ones, eta1);
    ctr->Zw1 = ctr->omega1.asDiagonal() * ctr->Z1;   

    // 1-3: Update Vg and z1
    // Vg
    Eigen::MatrixXd VgInv1(ctr->pZ1, ctr->pZ1); 
    VgInv1.triangularView<Eigen::Lower>() = ctr->Z1.transpose() * ctr->Zw1;  
    VgInv1.diagonal().array() += 1 / 100.0;                                
    VgInv1.triangularView<Eigen::Upper>() = VgInv1.transpose().eval(); 
    ctr->Vg1.triangularView<Eigen::Lower>() = VgInv1.inverse();
    ctr->Vg1.triangularView<Eigen::Upper>() = ctr->Vg1.transpose().eval();   
    ctr->VgChol1 = ctr->Vg1.llt().matrixL();

    // z1
    ctr->z1 = ((ctr->w).array() - 0.5).array() / ctr->omega1.array();

    // 1-4: Calculate the mean of gamma_1 & Variance of gamma_1 with cholesky
    const Eigen::VectorXd ZR1 = ctr->Zw1.transpose() * (ctr->z1); 
    ctr->b1 = ctr->Vg1 * ZR1;
    ctr->b1.noalias() += ctr->VgChol1 * as<Eigen::VectorXd>(rnorm(ctr->pZ1, 0, sqrt(ctr->sigma2)));

    // *** Step 2: Update ZI indicator auxiliary variable, w ***
    // 2-1: Calculate eta1 & logit1 (ZI), eta2 & logit2 (NB)
    eta1 = (ctr->Z1 * ctr->b1);
    Eigen::VectorXd eta2 = (ctr->Z * ctr->b2) + (ctr->fhat);
    Eigen::VectorXd logit1 = 1 / (1 + exp(-(eta1).array()));
    Eigen::VectorXd logit2 = 1 / (1 + exp(-(eta2).array()));
    
    // 2-2: Sampling w only for observations with y = 0
    for (int i = 0; i < (ctr->yZeroN); i++){ 
      int idx = (ctr->yZeroIdx)[i];
      double prob = log(logit1[idx]) - log(pow(1 - logit2[idx], ctr->r) * (1 - logit1[idx]) + logit1[idx]);
      (ctr->w)[idx] = R::rbinom(1, exp(prob));  // Update the index with the probability
    }

    // Update the number of at-risk individuals
    ctr->nStar = ctr->n - (ctr->w).sum();       

    // 2-3: Update the indices of at-risk
    ctr->NBidx.clear();  
    for (int j = 0; j < (ctr->n); j++){ 
      if(!(ctr->w)[j]){ 
        ctr->NBidx.push_back(j);
      }
    }

    // *** Step 3: Update the dispersion parameter, r ***
    // 3-1: Propose r with a random walk
    int rP;
    if(R::runif(0, 1) < 0.5){ 
      rP = ctr->r - 1;
    } else { 
      rP = ctr->r + 1;
    }

    // 3-2: Metropolis-Hastings
    ctr->MHratio = 0; 
    if(rP > 0){
      for(int k = 0; k < (ctr->nStar); k++){ 
        int idx_aR = (ctr->NBidx)[k];
        ctr->MHratio += R::dnbinom((ctr->Y0)[idx_aR], rP, logit2[idx_aR], true) -
                        R::dnbinom((ctr->Y0)[idx_aR], ctr->r, logit2[idx_aR], true);
      }

      // Accept / Reject 
      if(log(R::runif(0, 1)) < ctr->MHratio){
        ctr->r = rP;
        ctr->rVec = (ctr->ones).array() * (ctr->r);
      }
    }

    // *** Step 4: Update gamma_2 (Negative Binomial coefficients) ***
    // 4-1: Update omega2 ~ PG(y + r, eta2 + f)
    (ctr->omega2).setOnes(); 
    for(int l = 0; l < ctr->nStar; l++){
      int idx_NB = (ctr->NBidx)[l];
      (ctr->omega2)[idx_NB] = samplepg_na((ctr->Y0)[idx_NB] + ctr->r, eta2[idx_NB]);
    }

    // 4-2: Update Zstar, Zw, Vg, z2, R
    ctr->Zstar = (ctr->Z).array().colwise() * (1 - ctr->w.array());

    // Zw
    ctr->Zw = (ctr->omega2).asDiagonal() * ctr->Zstar; 

    // Vg
    Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ); 
    VgInv.setZero(); 
    VgInv.triangularView<Eigen::Lower>() = ctr->Zstar.transpose() * (ctr->Zw); 
    VgInv.diagonal().array() += 1 / 100.0; 
    VgInv.triangularView<Eigen::Upper>() = VgInv.transpose().eval(); 
    ctr->Vg.triangularView<Eigen::Lower>() = VgInv.inverse();
    ctr->Vg.triangularView<Eigen::Upper>() = ctr->Vg.transpose().eval();   
    ctr->VgChol = ctr->Vg.llt().matrixL();

    // z2 (= Ystar)
    ctr->z2 = (ctr->Y0 - ctr->rVec).array() / (2*(ctr->omega2).array()).array();
    ctr->Ystar = (ctr->z2).array() * (1 - ctr->w.array());

    // R (partial residual)
    Eigen::VectorXd fhatStar = ctr->fhat.array() * (1 - ctr->w.array()); 
    ctr->R = ctr->Ystar - fhatStar;

    // 4-3: Sample gamma_2
    const Eigen::VectorXd ZR = ctr->Zw.transpose() * (ctr->R);
    ctr->b2 = ctr->Vg * ZR;
    ctr->b2.noalias() += ctr->VgChol * as<Eigen::VectorXd>(rnorm(ctr->pZ, 0, sqrt(ctr->sigma2)));
  } // End ZINB
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
double tdlmProposeTree(Node* tree, exposureDat* Exp, modelCtr* ctr, int step,
                       double depth)
{
  int no = 0;
  double stepMhr = 0;
  std::vector<Node*> dlnmTerm, tempNodes;

  // List current tree terminal nodes
  dlnmTerm = tree->listTerminal();

  // Grow
  if (step == 0) {
    // select node to grow
    no = (std::size_t) floor(R::runif(0, dlnmTerm.size())); // Uniform selection among terminal nodes

    if (dlnmTerm[no]->grow()) { // propose new split
      double nGen2 = double(tree->nGen2());
      if (dlnmTerm[no]->depth == 0) { // If selected terminal node is the root node,
        ++nGen2;
      } else {
        if (!(dlnmTerm[no]->parent->isGen2())) {
          ++nGen2;
        }
      }
      
      stepMhr = log((double)tree->nTerminal()) - log(nGen2) +
        2 * logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth + depth + 1, 1) +
        logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth + depth, 0) -
        logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth + depth, 1);
      if (Exp != 0)
        Exp->updateNodeVals((dlnmTerm[no]->proposed)->c1); // update node values
      
      // newDlnmTerm = tree->listTerminal(1); // list proposed terminal nodes
    }


  // Prune
  } else if (step == 1) {
    tempNodes = tree->listGen2();
    no = floor(R::runif(0, tempNodes.size())); // select gen2 node to prune

    stepMhr = log((double)tree->nGen2()) - log((double)tree->nTerminal() - 1.0) -
      2 * logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], tempNodes[no]->depth + depth + 1, 1) -
      logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], tempNodes[no]->depth + depth, 0) +
      logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], tempNodes[no]->depth + depth, 1);

    tempNodes[no]->prune(); // prune nodes
    // newDlnmTerm = tree->listTerminal(1); // list proposed terminal nodes


  // Change
  } else {
    tempNodes = tree->listInternal();
    no = floor(R::runif(0, tempNodes.size())); // select internal nodes to change 
    if (tempNodes[no]->change()) { // propose new split
      for (Node* tn : tempNodes[no]->proposed->listTerminal()) {
        if (Exp != 0)
          Exp->updateNodeVals(tn);
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
  if (n->depth == 0)  // One root node -> no rules to the terminal nodes -> return
    return(rule);
  
  // If there is a node to ther terminal nodes, store a parent node
  Node* parent = n->parent;
  int splitVar = parent->nodestruct->get(1);  // What variable of the parent?
  int splitVal = parent->nodestruct->get(2);  // What value of the parent?
  std::vector<int> splitVec = parent->nodestruct->get2(1); // Return split vector
  
  rule += std::to_string(splitVar);       // Convert to a string and concatenate
  if (Mod->varIsNum[splitVar]) {          // [If continuous], 
    if (parent->c1 == n)                  // If the first child node is the same as the node, (c1 = child node 1)
      rule += "<";                        // Add a rule
    else                                  
      rule += ">=";                       // Add a rule
    rule += std::to_string(splitVal);     // Add the splitting value
  } else {
    if (parent->c1 == n)                  // [If categorical],
      rule += "[]";                       // Add a subsetting rule
    else
      rule += "][";                       // Add a subsetting rule
    for (int i : splitVec)
      rule += std::to_string(i) + ",";    // 
    rule.pop_back();                      // Remove the last character of the string
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
    // if (unavail.size() > 0) {
    //   std::shuffle(unavail.begin(), unavail.end(), std::default_random_engine());
    //   double totProb = unavailProb.sum();
    //   int pseudoDraw = R::rgeom(std::max(0.00000001, 1 - totProb));
    //   int binomDraw = 0;
    //   if (pseudoDraw > 0) {
    //     for (int i : unavail) {
    //       binomDraw = R::rbinom(pseudoDraw, unavailProb(i) / totProb);
    //       if (binomDraw > 0)
    //         modCount(i) += binomDraw * 1.0;
    //       totProb -= unavailProb(i);
    //       pseudoDraw -= binomDraw;
    //       if (pseudoDraw < 1)
    //         break;
    //     } // end multinom
    //   } // end pseudoDraw
    // } // end unavail
    // if (unavail.size() > 0) {
    //   std::shuffle(unavail.begin(), unavail.end(), std::default_random_engine());
    //   double totProb = unavailProb.sum();
    //   int pseudoDraw = R::rgeom(std::max(0.00000001, 1 - totProb));
    //   int binomDraw = 0;
    //   if (pseudoDraw > 0) {
    //     for (int i : unavail) {
    //       binomDraw = R::rbinom(pseudoDraw, unavailProb(i) / totProb);
    //       if (binomDraw > 0)
    //         modCount(i) += binomDraw * 1.0;
    //       totProb -= unavailProb(i);
    //       pseudoDraw -= binomDraw;
    //       if (pseudoDraw < 1)
    //         break;
    //     } // end multinom
    //   } // end pseudoDraw
    // } // end unavail
  } // end modCount
  return(modCount);
} // end countMods function







void updateZirtGamma(std::vector<Node*> trees, modelCtr* ctr) {
  ctr->zirtSplitCounts.setZero();

  // count terminal nodes in time trees
  int nTerm = 0;
  for (Node* tree : trees)
    nTerm += tree->nTerminal();

  // Create split outcome (zY) and design matrix (zX)
  VectorXd zY(nTerm); zY.array() = -0.5; // polya-gamma y 
  MatrixXd zX(nTerm, ctr->pX); zX.setZero(); // times
  int i = 0;

  // loop over tree terminal nodes to determine times and nonzero effect
  for (Node* tree : trees) {
    for (Node* eta : tree->listTerminal(0)) { 
      // set indices of x matrix equal to 1 corresponding to terminal node times
      int tmin = eta->nodestruct->get(3);
      int tmax = eta->nodestruct->get(4);
      zX.row(i).segment(tmin - 1, tmax - tmin + 1).array() = 1.0 / (tmax - tmin + 1);

      // determine if node has nonzero effect, set y = 1 (kappa = 1/2)
      if (eta->nodevals->nestedTree->c1 != 0) {
        zY(i) += 1.0;
        ctr->zirtSplitCounts.segment(tmin - 1, tmax - tmin + 1).array() += 1.0;
      }

      // increase iterator for X and y
      ++i;
    }
  }
      
  // draw polya-gamma for CW var selection
  VectorXd zOnes(nTerm); zOnes.setOnes();
  MatrixXd zV;
  VectorXd psi, zPG;
  psi = zX * ctr->zirtGamma;
  zPG = rcpp_pgdraw(zOnes, psi);
  zV = zX.transpose() * zPG.asDiagonal() * zX + ctr->zirtSigma;
  ctr->zirtGamma = zV.inverse() * (zX.transpose() * zY + ctr->zirtSigma * ctr->zirtGamma0);
  ctr->zirtGamma += zV.inverse().llt().matrixL() * as<VectorXd>(rnorm(ctr->pX, 0, 1));
} // end updateZirtGamma







int updateZirtSigma(std::vector<Node*> trees, modelCtr* ctr, 
 int curCov, std::vector<MatrixXd> zirtSigmaInv, 
 std::vector<double> zirtSigmaDet) {

  double covMHR = 0.0;
  int newCov = curCov;

  // step to next covariance matrix unless at limit (0 and 18)
  if (curCov == 0) {
    newCov = 1;
    covMHR += log(0.5);
  } else if (curCov == 18) {
    newCov = 17;
    covMHR += log(0.5);
  } else {
    if (R::runif(0.0, 1.0) < 0.5)
      ++newCov;
    else
      --newCov;
  } // end step to next covariance matrix

  // calculate MHR for change in covariance matrix
  covMHR += -zirtSigmaDet[newCov] - 
    (ctr->zirtGamma - ctr->zirtGamma0).dot(zirtSigmaInv[newCov] * (ctr->zirtGamma - ctr->zirtGamma0)) + 
    zirtSigmaDet[curCov] + (ctr->zirtGamma - ctr->zirtGamma0).dot(zirtSigmaInv[curCov] * (ctr->zirtGamma - ctr->zirtGamma0));
  if (log(R::runif(0.0, 1.0)) < covMHR) {
    curCov = newCov;
    ctr->zirtSigma = zirtSigmaInv[curCov];
  }
  return(curCov);
} // end updateZirtSigma



void updateTimeSplitProbs(std::vector<Node*> trees, modelCtr* ctr) {
  ctr->timeSplitCounts.setZero();

  // count split locations
  // Rcout << 1;
  for (Node* tree : trees) {
    for (Node* tn : tree->listInternal()) {
      if (tn->nodestruct->get(6) == 0)
        continue;
      ctr->timeSplitCounts(tn->nodestruct->get(6) - 1) += 1.0;
    }
  }

  // update time split probabilities
  // Rcout << 2;
  ctr->timeSplitProbs = rDirichlet(ctr->timeSplitCounts + ctr->timeKappa * ctr->timeSplitProb0);
  // Rcout << ctr->timeSplitProbs;
  for (Node* tree : trees) {
    tree->nodestruct->setTimeProbs(ctr->timeSplitProbs);
    tree->updateStruct();
  }

  // update dirichlet scaling factor (if not fixed)
  // Rcout << 3;
  if (ctr->updateTimeKappa) {
    double beta = R::rbeta(1.0, 1.0);
    double timeKappaNew = beta * (ctr->pX - 1.0)/ (1 - beta);
    double mhrDir = 
      logDirichletDensity(ctr->timeSplitProbs, ctr->timeSplitCounts + timeKappaNew * ctr->timeSplitProb0) - 
      logDirichletDensity(ctr->timeSplitProbs, ctr->timeSplitCounts + ctr->timeKappa * ctr->timeSplitProb0);
    if (log(R::runif(0, 1)) < mhrDir)
      ctr->timeKappa = timeKappaNew;
  }
} // end updateTimeSplitProbs


/**
 * @brief draw tree structure from tree prior distribution
 * 
 * @param tree pointer to tree being grown
 * @param n pointer to specific node in tree
 * @param alpha alpha parameter (0 to 1)
 * @param beta beta parameter (>0)
 * @param depth added to depth in p_split
 */
void drawTree(Node* tree, Node* n, double alpha, double beta, 
              double depth)
{
  double logProb = log(alpha) - beta * log(1.0 + depth + n->depth);
  if (log(R::runif(0, 1)) < logProb) {
    if (n->grow()) {
      if (n->depth > 0)
        n = n->proposed;
      tree->accept();
      drawTree(tree, n->c1, alpha, beta, depth);
      drawTree(tree, n->c2, alpha, beta, depth);
    }
  } // end grow tree
  return;
} // end drawTree function

/**
 * @brief 
 * 
 * @param eta 
 * @param ctr 
 * @param nsX 
 */
void drawZirt(Node* eta, tdlmCtr* ctr, NodeStruct* nsX)
{
  int tmin = eta->nodestruct->get(3);
  int tmax = eta->nodestruct->get(4);
  eta->nodevals->nestedTree = new Node(0, 1);
  eta->nodevals->nestedTree->nodestruct = nsX->clone();
  eta->nodevals->nestedTree->nodestruct->setTimeRange(tmin, tmax);
  
  double logProb = logZIPSplit(ctr->zirtGamma, tmin, tmax, ctr->nTrees, 0);  
  if (log(R::runif(0, 1)) < logProb) {
    if (eta->nodevals->nestedTree->grow()) {
      eta->nodevals->nestedTree->accept();
      drawTree(eta->nodevals->nestedTree, eta->nodevals->nestedTree->c1, ctr->treePrior2[0], ctr->treePrior2[1], 0.0);
      drawTree(eta->nodevals->nestedTree, eta->nodevals->nestedTree->c2,   ctr->treePrior2[0], ctr->treePrior2[1], 0.0);
    }
  } // end grow tree
  return;
} // end drawTree function


double zeroInflatedTreeMHR(VectorXd timeProbs, std::vector<Node*> trees, int t, double newProb)
{
  double mhr = 0.0;
  int tmin, tmax;
  VectorXd timeProbsNew = timeProbs;
  timeProbsNew(t) = newProb;
  
  for (Node* tree : trees) { // loop over all trees
    for (Node* eta : tree->listTerminal(0)) { // loop over tree terminal nodes
      tmin = eta->nodestruct->get(3);
      tmax = eta->nodestruct->get(4);
      
      if ((t >= tmin - 1) && (t < tmax)) { // check time within range        
        if (eta->nodevals->nestedTree->c1 == 0) { // single node tree
          mhr += logZIPSplit(timeProbsNew, tmin, tmax, trees.size(), 1) -
            logZIPSplit(timeProbs, tmin, tmax, trees.size(), 1);
            
        } else {
          mhr += logZIPSplit(timeProbsNew, tmin, tmax, trees.size(), 0) -
            logZIPSplit(timeProbs, tmin, tmax, trees.size(), 0);
            
        } // end if single node tree
      } // end time within range
    } // end loop over terminal nodes
  } // end loop over trees
  return (mhr);
} // end zeroInflatedTreeMHR function

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