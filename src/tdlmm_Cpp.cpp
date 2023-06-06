/**
 * @file tdlmmGaussian.cpp
 * @author Daniel Mork (danielmork.github.io)
 * @brief Method for TDLMM
 * @version 1.0
 * @date 2021-03-02
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <RcppEigen.h> // Linear algebra package
#include "modelCtr.h"
#include "exposureDat.h"
#include "Node.h"
#include "NodeStruct.h"
#include "Fncs.h" // Helpful functions
#include <random>
#include <iostream>
using namespace Rcpp;
using namespace Eigen;

/**
 * @brief 
 * 
 * @param nodes1 // Terminal nodes of tree1
 * @param nodes2 // Terminal nodes of tree2
 * @param ctr 
 * @param ZtR 
 * @param treeVar 
 * @param m1Var
 * @param m2Var
 * @param mixVar 
 * @param tree 
 * @param newTree 
 * @return treeMHR 
 */
treeMHR mixMHR(std::vector<Node*> nodes1, std::vector<Node*> nodes2,
                  tdlmCtr *ctr, Eigen::VectorXd ZtR,
                  double treeVar, double m1Var, double m2Var, double mixVar,
                  Node* tree, bool newTree)
{
  //Rcout << "Initiating mixMHR \n";
  treeMHR out; // treeMHR object

  // This part refers to Eq(6) & (7)
  int pX1 = nodes1.size();  // Number of terminal nodes for tree1
  int pX2 = nodes2.size();  // Number of terminal nodes for tree2
  int pXd = pX1 + pX2;      // pXd = sum of sizes
    
  int interaction = 0;      // interaction
  
  if (mixVar != 0) {  // If mix variable is not zero (there is an interaction),
    pXd += pX1 * pX2; // Add (nodes x nodes) to the parameter # (the rectangle figure)
    interaction = 1;  // Interaction flag
  }
  
  // Rcout << ".";
  out.Xd.resize(ctr->n, pXd);         out.Xd.setZero();   // This is X*delta vector for mean in Eq(8)
  Eigen::MatrixXd ZtX(ctr->pZ, pXd);  ZtX.setZero();      // Create a matrix of (predictor# x (total terminal node#))
  Eigen::VectorXd diagVar(pXd);       diagVar.setZero();  // Create a vector of length (parameter size)

  // ****** Constructing U_a (Exposure specific variance diagonal matrix) and computing Eq. (9) ****** 
  // Each column of out.Xd represents a terminal node 
  int i, j, k;
  // ----- Update the calculations for the tree pair 1: mu ----- 
  // Rcout << "Updating the calculations for Nodes 1 \n";
  for (i = 0; i < pX1; ++i) { // Partitions of design matrix
    out.Xd.col(i) = (nodes1[i]->nodevals)->X; // Each column of Xd is a vector X associated with a node from nodes1
    diagVar(i) = 1.0 / (m1Var * treeVar);

    // For binomial, we have to implement omega, For ZINB, Zw is already At-risk individuals only
    if (ctr->binomial || ctr->zinb) { 
      ZtX.col(i) = ctr->Zw.transpose() * (nodes1[i]->nodevals)->X; // (pZ x 1) = (n x pZ).transpose x (n x 1)
    } else { // Gaussian
      ZtX.col(i) = (nodes1[i]->nodevals)->ZtX;
    }
  }

  // Rcout << "Updating the calculations for the tree pair 2 \n";
  // ----- Update the calculations for Nodes2 ----- 
  for (j = 0; j < pX2; ++j) {
    k = pX1 + j;                                // k is pushed back because pX2 needs to be after pX1.
    out.Xd.col(k) = (nodes2[j]->nodevals)->X;
    diagVar(k) = 1.0 / (m2Var * treeVar);
    if (ctr->binomial || ctr->zinb){
      ZtX.col(k) = ctr->Zw.transpose() * (nodes2[j]->nodevals)->X;
    } else {
      ZtX.col(k) = (nodes2[j]->nodevals)->ZtX;
    }
  }

  // ----- Update the calculations for interaction ----- 
  if(interaction) {
    for (i = 0; i < pX1; ++i) {
      for (j = 0; j < pX2; ++j) {
        k = pX1 + pX2 + i * pX2 + j;
        out.Xd.col(k) = (((nodes1[i]->nodevals)->X).array() * ((nodes2[j]->nodevals)->X).array()).matrix();
        diagVar(k) = 1.0 / (mixVar * treeVar);
        ZtX.col(k) = ctr->Zw.transpose() * out.Xd.col(k);
      }
    }
  }

  // Rcout << ".";
  // ----- calculate MHR -----
  const Eigen::MatrixXd VgZtX = ctr->Vg * ZtX; // V_gamma * Z_t * X = (p x p) * (p x n) * (n x 1), w multiplied for ZINB
  Eigen::MatrixXd tempV(pXd, pXd); // 
  Eigen::VectorXd XtVzInvR(ctr->n);
  if (ctr->binomial) { // V_theta = V_delta_a from eq(11) & Xd = X_a in paper
    const Eigen::MatrixXd Xdw = (ctr->Omega).asDiagonal() * out.Xd;     // (nxn) x (nx2) = (nx2) 
    tempV = Xdw.transpose() * out.Xd;                                   // (2xn) x (nx2) = (2x2)
    tempV.noalias() -= ZtX.transpose() * VgZtX;                         // (2x2) - (2x5) x (5x2)
    XtVzInvR = Xdw.transpose() * ctr->R;                                // (2xn) x (nx1) = (2x1)

  } else if (ctr->zinb){ // ZINB (subsetting at-risk observations)
    Eigen::MatrixXd outXdstar = selectIndM(out.Xd, ctr->atRiskIdx);        // At-risk observations only
    const Eigen::MatrixXd Xdw = (selectInd(ctr->omega2, ctr->atRiskIdx)).asDiagonal() * outXdstar; // (n*xn*) x (n*x2) = (n*x2) 
    tempV = Xdw.transpose() * outXdstar;                                   // (2xn*) x (n*x2) = (2x2)
    tempV.noalias() -= ZtX.transpose() * VgZtX;                            // (2x2) - (2x5) x (5x2)
    XtVzInvR = Xdw.transpose() * selectInd(ctr->R, ctr->atRiskIdx);       // (2xn*) x (n*x1) = (2x1)

    // Eigen::MatrixXd outXdstar = out.Xd.array().colwise() * (ctr->w).array();        // At-risk observations only
    // const Eigen::MatrixXd Xdw = ctr->omega2.asDiagonal() * outXdstar; // (n*xn*) x (n*x2) = (n*x2) 
    // tempV = Xdw.transpose() * outXdstar;                                   // (2xn*) x (n*x2) = (2x2)
    // tempV.noalias() -= ZtX.transpose() * VgZtX;                            // (2x2) - (2x5) x (5x2)
    // XtVzInvR = Xdw.transpose() * ctr->R;       // (2xn*) x (n*x1) = (2x1)
  } else {
    if (newTree) { // V_theta = V_delta_a from eq(11) & Xd = X_a in paper
      tempV.triangularView<Eigen::Lower>() = out.Xd.transpose() * out.Xd;
      tempV.noalias() -= ZtX.transpose() * VgZtX;
      out.tempV = tempV;
    } else {
      tempV = tree->nodevals->tempV;
    }
    XtVzInvR = out.Xd.transpose() * ctr->R;
  }

  // Finalize Eq.(9)
  XtVzInvR.noalias() -= VgZtX.transpose() * ZtR; // (2x1) - (2x5) * (5x1) = (2x1)
  tempV.diagonal().noalias() += diagVar;


  // Rcout << ".";
  Eigen::MatrixXd VTheta(pXd, pXd); // V_theta = V_delta_a from eq(11) & Xd = X_a in paper
  VTheta.triangularView<Eigen::Lower>() =
    tempV.selfadjointView<Eigen::Lower>().llt().solve(
        Eigen::MatrixXd::Identity(pXd, pXd));
  const Eigen::MatrixXd VThetaChol =
    VTheta.selfadjointView<Eigen::Lower>().llt().matrixL();

  // This part is sampling from Supplemental Eq(11) & (12) (Theta = Delta_a)
  // Mean
  const Eigen::VectorXd ThetaHat =
    VTheta.selfadjointView<Eigen::Lower>() * XtVzInvR; // Eq. (11) mean
  Eigen::VectorXd ThetaDraw = ThetaHat;

  // Add variance
  ThetaDraw.noalias() += VThetaChol * as<Eigen::VectorXd>(rnorm(pXd, 0, sqrt(ctr->sigma2)));

  out.drawAll = ThetaDraw;                     // (px1)
  out.draw1 = ThetaDraw.head(pX1);
  out.term1T2 = (out.draw1).dot(out.draw1);
  out.nTerm1 = double(pX1);
  out.draw2 = ThetaDraw.segment(pX1, pX2);
  out.term2T2 = (out.draw2).dot(out.draw2);
  out.nTerm2 = double(pX2);

  if (interaction) {
    out.drawMix = ThetaDraw.tail(pXd - pX1 - pX2);  // Extract last pX1 x pX2 element for mixture
    out.mixT2 = (out.drawMix).dot(out.drawMix);     // dot product
  }

  out.beta = ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
  out.pXd = pXd;
  return(out);
}

/**
 * @brief 
 * 
 * @param t     // Index for 't'th tree
 * @param tree1 // Tree1 from Trees1[t]
 * @param tree2 // Tree2 from Trees2[t]
 * @param ctr   // model control object
 * @param dgn   // Model logs
 * @param Exp   // Exposure data
 */
void tdlmmTreeMCMC(int t, Node *tree1, Node *tree2, tdlmCtr *ctr, tdlmLog *dgn,
                   std::vector<exposureDat*> Exp)
{
  // Rcout << "tdlmmTreeMCMC: Initialize \n";
  int m1, m2, newExp, success, step1, step2;
  double stepMhr, ratio;
  double m1Var, m2Var, mixVar, newExpVar, newMixVar, treeVar;
  double RtR = -1.0;
  double RtZVgZtR = 0;
  std::vector<Node*> term1, term2, newTerm;
  Node* newTree = 0;
  treeMHR mhr0, mhr;

  // Rcout << "Mixture calculation \n";
  term1 = tree1->listTerminal();        // List the terminal nodes of tree 1
  term2 = tree2->listTerminal();        // List the terminal nodes of tree 2
  treeVar = (ctr->nu) * (ctr->tau[t]);  // Global shrinkage x Local Shrinkage
  m1 = ctr->tree1Exp[t];                // Exposure m1 of tree 1 of t th tree pair // Note that tree1Exp is a vector of exposures with a designated number
  m2 = ctr->tree2Exp[t];                // Exposure m2 of tree 2 of t th tree pair
  m1Var = ctr->muExp(m1);               // Exposure-specific variance
  m2Var = ctr->muExp(m2);               // Exposure-specific variance
  mixVar = 0;                           // No interaction variance for now

  // If there is an interaction, extract mixture variance value from the muMix matrix
  if ((ctr->interaction) && ((ctr->interaction == 2) || (m1 != m2))) { 
    if (m1 <= m2)
      mixVar = ctr->muMix(m2, m1); // Interaction rectangles
    else
      mixVar = ctr->muMix(m1, m2); // Interaction rectangles
  }

  // Compute t(Z) * (Y - fhat)
  Eigen::VectorXd ZtR = (ctr->Zw).transpose() * (ctr->R); // (pxn) x (nx1) = (px1)

  // * Update tree 1 ------------------------------------------------------------
  newExp    = m1; 
  newExpVar = m1Var;
  newMixVar = mixVar;
  stepMhr   = 0;
  success   = 0;
  
  // Rcout << ".";
  // * List current tree terminal nodes
  // 1) grow/prune, 2) change, 3) switch exposure, defaults to (0.25, 0.25, 0.25)
  // Choose a step by sampling one integer
  step1 = sampleInt(ctr->stepProb, 1); 

  // If there is only one terminal node and grow/prune/change is selected, the tree can only grow hence step1 = 0 
  if ((term1.size() == 1) && (step1 < 3))
    step1 = 0;

  // * Propose update:
  // Grow(0) / Prune(1) / Change(2)
  if (step1 < 3) { 
    stepMhr = tdlmProposeTree(tree1, Exp[m1], ctr, step1); // tdlmProposeTree returns MH ratio for each possible step
    success = tree1->isProposed();
    newTerm = tree1->listTerminal(1);

  // * Switch exposures (3)
  } else {
    newExp = sampleInt(ctr->expProb); // Propose a new exposure
    if (newExp != m1) { // If the proposed exposure is different from the current one,
      // Update the information with new exposure
      success   = 1;
      newExpVar = ctr->muExp(newExp);
      newTree   = new Node(*tree1);
      newTree->setUpdate(1);
      newTerm   = newTree->listTerminal();

      for (Node* nt : newTerm)
        Exp[newExp]->updateNodeVals(nt);

      // Update the interaction using the new exposure as well
      if ((ctr->interaction) && ((ctr->interaction == 2) || (newExp != m2))) {
        if (newExp <= m2)
          newMixVar = ctr->muMix(m2, newExp);
        else
          newMixVar = ctr->muMix(newExp, m2);
      } else {
        newMixVar = 0;
      }
    }
  }  // Propose update end

  // * Tree 1 MHR
  // StepMHR: Transition ratio
  // mhr0 & mhr: Likelihood of R
  // mhr0: Current MHR
  // Rcout << "a";
  if ((tree1->nodevals->tempV).rows() == 0)
    mhr0 = mixMHR(term1, term2, ctr, ZtR, treeVar,  // 0 = denominator(old state)
                  m1Var, m2Var, mixVar, tree1, 1);
  else
    mhr0 = mixMHR(term1, term2, ctr, ZtR, treeVar, 
                  m1Var, m2Var, mixVar, tree1, 0);

  //Rcout << "Tree 1 MHR: if success\n";
  if (success) {
    // New tree terminal -> newTerm & newExpVar
    mhr = mixMHR(newTerm, term2, ctr, ZtR, treeVar,  // numerator(new or updated state)
                 newExpVar, m2Var, newMixVar, tree1, 1);
    // combine mhr parts into log-MH ratio
    //Rcout << "Tree 1 MHR: log-MH ratio\n";
    if (ctr->binomial || ctr->zinb) {
      ratio = stepMhr +                                 
              mhr.logVThetaChol - mhr0.logVThetaChol +  // Combine mhr and mhr0 to obtain Eq (10)
              0.5 * (mhr.beta - mhr0.beta) -
              (0.5 * ((log(treeVar * newExpVar) * mhr.nTerm1) -
              (log(treeVar * m1Var) * mhr0.nTerm1)));
    } else { // Gaussian
      // Rcout << "Tree 1 MHR: Gaussian\n";
      // Equation (9)
      if (RtR < 0) {
        RtR = (ctr->R).dot(ctr->R);                                             // Rt * I * R
        RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);    // Rt * (-ZtVgZ) * R
      }
      ratio = stepMhr +                                           // stepMhr = p(T|T*) / p(T*|T)
                mhr.logVThetaChol - mhr0.logVThetaChol -          // |V_delta_a|, 1/2 canceled out in Eq(10)
                (0.5 * (ctr->n + 1.0) *                           // RtVinvXVXtVintR + 1/xi in Eq.(9) & Eq.(10)
                (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
                log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
                (0.5 * ((log(treeVar * newExpVar) * mhr.nTerm1) - // -0.5 * (B1 * (log(nu) + log(mu of new Exposure)))
                (log(treeVar * m1Var) * mhr0.nTerm1)));           // -0.5 * (B1 * (log(nu) + log(mu1)))
    }
    // Mixture
    if (newMixVar != 0)
      ratio -= 0.5 * log(treeVar * newMixVar) * mhr.nTerm1 * mhr0.nTerm2; // new mu1 (new B1 x old B2) - Mixture update with switch exposure
    if (mixVar != 0)
      ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2;   // old mu1 (old B1 x old B2) - Mixture

    if (log(R::runif(0, 1)) < ratio) { // Evaluate the ratio
      mhr0 = mhr; // Update the MHR
      success = 2;

      // If accepted with switch-exposure transition, update a tree more
      if (step1 == 3) { 
        m1      = newExp;
        m1Var   = newExpVar;
        mixVar  = newMixVar;
        tree1->replaceNodeVals(newTree);
      } else {
        tree1->accept();
      }
      if (!(ctr->binomial) && !(ctr->zinb)) { // For Gaussian approach,
        (tree1->nodevals->tempV).resize(mhr0.pXd, mhr0.pXd);
        tree1->nodevals->tempV = mhr0.tempV;
      }
      term1 = tree1->listTerminal();

    } else {
      tree1->reject();
    }

  } else if (step1 < 3) {
      tree1->reject();
  }

  if (newTree != 0)
    delete newTree;
  newTree = 0;
  //Rcout << "Recording Tree 1\n";
  // * Record tree 1
  if (ctr->diagnostics) {
    Eigen::VectorXd acc(7);
    acc << 1, step1, success, m1, term1.size(), stepMhr, ratio;
    (dgn->TreeAccept).push_back(acc);
  }
  //Rcout << "Updating Tree 2\n";
  // * Update tree 2 ------------------------------------------------------------
  newExp    = m2;
  newExpVar = m2Var;
  newMixVar = mixVar;
  stepMhr   = 0;
  success   = 0;

  // * List current tree terminal nodes
  step2 = sampleInt(ctr->stepProb, 1);
  if ((term2.size() == 1) && (step2 < 3))
    step2 = 0;

  // * Propose update
  if (step2 < 3) {
    // Rcout << "g" << step2;
    stepMhr = tdlmProposeTree(tree2, Exp[m2], ctr, step2);
    success = tree2->isProposed();
    newTerm = tree2->listTerminal(1);

  // * Switch exposures
  } else {
    newExp = sampleInt(ctr->expProb);
    if (newExp != m2) {
      success   = 1;
      newExpVar = ctr->muExp(newExp);
      newTree   = new Node(*tree2);
      newTree->setUpdate(1);
      newTerm   = newTree->listTerminal();
      for (Node* nt : newTerm)
        Exp[newExp]->updateNodeVals(nt);


      if ((ctr->interaction) && ((ctr->interaction == 2) || (newExp != m1))) {
        if (newExp <= m1) {
          newMixVar = ctr->muMix(m1, newExp);
        } else {
          newMixVar = ctr->muMix(newExp, m1);
        }
      } else {
        newMixVar = 0;
      }
    }
  }

  if (success) {
    // calculate new mhr part
    // Rcout << "c";
    mhr = mixMHR(term1, newTerm, ctr, ZtR, treeVar, 
                 m1Var, newExpVar, newMixVar, tree1, 1);
    // combine mhr parts into log-MH ratio
    if (ctr->binomial || ctr->zinb) {
      ratio = stepMhr + mhr.logVThetaChol - mhr0.logVThetaChol +
        0.5 * (mhr.beta - mhr0.beta) -
        (0.5 * ((log(treeVar * newExpVar) * mhr.nTerm2) -
         (log(treeVar * m2Var) * mhr0.nTerm2)));
    } else {
      if (RtR < 0) {
        RtR = (ctr->R).dot(ctr->R);
        RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);
      }
      ratio = stepMhr + mhr.logVThetaChol - mhr0.logVThetaChol -
        (0.5 * (ctr->n + 1.0) *
         (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
          log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
        (0.5 * ((log(treeVar * newExpVar) * mhr.nTerm2) -
         (log(treeVar * m2Var) * mhr0.nTerm2)));
    }

    if (newMixVar != 0)
      ratio -= 0.5 * log(treeVar * newMixVar) * mhr0.nTerm1 * mhr.nTerm2;
    if (mixVar != 0)
      ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2;

    if (log(R::runif(0, 1)) < ratio) {
      mhr0    = mhr;
      success = 2;

      if (step2 == 3) {
        m2      = newExp;
        m2Var   = newExpVar;
        mixVar  = newMixVar;
        tree2->replaceNodeVals(newTree);

      } else {
        tree2->accept();
      }
      if (!(ctr->binomial) && !(ctr->zinb)) {
        (tree1->nodevals->tempV).resize(mhr0.pXd, mhr0.pXd);
        tree1->nodevals->tempV = mhr0.tempV;
      }
      term2 = tree2->listTerminal();

    } else {
      tree2->reject();
    }

  } else if (step2 < 3) {
    tree2->reject();
  }

  if (newTree != 0)
    delete newTree;
  newTree = 0;

  // Rcout << "Recording tree 2\n";
  // * Record tree 2
  if (ctr->diagnostics) {
    Eigen::VectorXd acc(7);
    acc << 2, step2, success, m2, term2.size(), stepMhr, ratio;
    (dgn->TreeAccept).push_back(acc);
  }

  // * Update variance and residuals ----------------------------------------
  double tauT2 = mhr0.term1T2 / m1Var + mhr0.term2T2 / m2Var;
  int totTerm = mhr0.nTerm1 + mhr0.nTerm2;
  if (mixVar != 0) {
    tauT2 += mhr0.mixT2 / mixVar;
    totTerm += mhr0.nTerm1 * mhr0.nTerm2;
  }
  // Rcout << "d";
  if (ctr->shrinkage > 1)
    rHalfCauchyFC(&(ctr->tau(t)), totTerm, tauT2 / (ctr->sigma2 * ctr->nu));
  
  // if ((ctr->tau)(t) != (ctr->tau)(t)) 
  //   stop("\nNaN values (tau) occured during model run, rerun model.\n");
    
  (ctr->nTerm)(t) = mhr0.nTerm1;
  (ctr->nTerm2)(t) = mhr0.nTerm2;
  (ctr->tree1Exp)(t) = m1;
  (ctr->tree2Exp)(t) = m2;
  (ctr->expCount)(m1)++;
  (ctr->expCount)(m2)++;
  (ctr->expInf)(m1) += ((ctr->tau)(t));
  (ctr->expInf)(m2) += ((ctr->tau)(t));
  (ctr->totTermExp)(m1) += mhr0.nTerm1;
  (ctr->totTermExp)(m2) += mhr0.nTerm2;
  (ctr->sumTermT2Exp)(m1) += mhr0.term1T2 / (ctr->tau)(t);
  (ctr->sumTermT2Exp)(m2) += mhr0.term2T2 / (ctr->tau)(t);
  if (mixVar != 0) {
    if (m1 <= m2) {
      (ctr->mixCount)(m2, m1)++;
      (ctr->totTermMix)(m2, m1) += mhr0.nTerm1 * mhr0.nTerm2;
      (ctr->sumTermT2Mix)(m2, m1) += mhr0.mixT2 / ctr->tau[t];
      (ctr->mixInf)(m2, m1) += ((ctr->tau)(t));
    } else {
      (ctr->mixCount)(m1, m2)++;
      (ctr->totTermMix)(m1, m2) += mhr0.nTerm1 * mhr0.nTerm2;
      (ctr->sumTermT2Mix)(m1, m2) += mhr0.mixT2 / ctr->tau[t];
      (ctr->mixInf)(m1, m2) += ((ctr->tau)(t));
    }
  }

  // Sample R_a (Eq. (8) in Supplementary)
  // t is an index for tree Ensemble
  (ctr->Rmat).col(t) = mhr0.Xd * mhr0.drawAll; 

  // * Record
  if (ctr->record > 0) {
    Eigen::VectorXd rec(8);
    Eigen::VectorXd mix(10);
    rec << ctr->record, t, 0, 0, 0, 0, 0, 0;
    mix << ctr->record, t, 0, 0, 0, 0, 0, 0, 0, 0;
    int k = 0;
    for(int i = 0; i < mhr0.nTerm1; ++i) {
      rec[2] = 0; // First of the tree pair
      rec[3] = m1;
      rec[4] = (term1[i]->nodestruct)->get(3);
      rec[5] = (term1[i]->nodestruct)->get(4);
      rec[6] = mhr0.draw1(i);
      rec[7] = (ctr->tau)(t) * m1Var;
      (dgn->DLMexp).push_back(rec);
      for (int j = 0; j < mhr0.nTerm2; ++j) {
        if (i == 0) {
          rec[2] = 1; // Second of the tree pair
          rec[3] = m2;
          rec[4] = (term2[j]->nodestruct)->get(3);
          rec[5] = (term2[j]->nodestruct)->get(4);
          rec[6] = mhr0.draw2(j);
          rec[7] = (ctr->tau)(t) * m2Var;
          (dgn->DLMexp).push_back(rec);
        }
        if (mixVar != 0) {
          if (m1 <= m2) {
            mix[2] = m1;
            mix[3] = (term1[i]->nodestruct)->get(3);
            mix[4] = (term1[i]->nodestruct)->get(4);
            mix[5] = m2;
            mix[6] = (term2[j]->nodestruct)->get(3);
            mix[7] = (term2[j]->nodestruct)->get(4);
          } else {
            mix[5] = m1;
            mix[6] = (term1[i]->nodestruct)->get(3);
            mix[7] = (term1[i]->nodestruct)->get(4);
            mix[2] = m2;
            mix[3] = (term2[j]->nodestruct)->get(3);
            mix[4] = (term2[j]->nodestruct)->get(4);
          }
          mix[8] = mhr0.drawMix(k);
          (dgn->MIXexp).push_back(mix);
          ++k;
        }
      }
    }
  }
} // end function tdlmmTreeMCMC

/**
 * @brief 
 * 
 * @param model 
 * @return Rcpp::List 
 */
// [[Rcpp::export]]
Rcpp::List tdlmm_Cpp(const Rcpp::List model)
{
  // ****** Set up model control parameters converting from R to C++ ******
  // Rcout << "Initiating tdlmm_Cpp function \n";

  // Model control(ctr) object with pointer
  tdlmCtr *ctr = new tdlmCtr; 
  
  // MCMC parameters
  ctr->iter = as<int>(model["nIter"]);        // Total MCMC iteration
  ctr->burn = as<int>(model["nBurn"]);        // MCMC burn-in iteration
  ctr->thin = as<int>(model["nThin"]);        // Thinning iteration
  ctr->nRec = floor(ctr->iter / ctr->thin);   // rounding to calculate thinned iteration

  // Tree parameters
  ctr->nTrees = as<int>(model["nTrees"]);                        // Number of trees for the ensemble
  ctr->stepProb = as<std::vector<double> >(model["stepProb"]);   // A vector of step probability for the tree structure update
  ctr->swapStep = as<bool>(model["swapStep"]); 
  ctr->treePrior = as<std::vector<double> >(model["treePrior"]); // Tree prior of vector length 2 c(alpha, beta)
  ctr->verbose = as<bool>(model["verbose"]);                     // Boolean for returning output
  ctr->diagnostics = as<bool>(model["diagnostics"]);             // Store diagnostic or not

  // Model selection
  ctr->binomial = as<bool>(model["binomial"]);  // Binomial
  ctr->zinb = as<bool>(model["zinb"]);          // ZINB

  // Mixture & Shrinkage
  ctr->modKappa = as<double>(model["mixPrior"]);  // Positive scalar hyperparameter for sparsity of exposures (Linero): Set as 1
  bool updateKappa = false;                       // Check to update modKappa: if negative, set to 1 like the default.
  if (ctr->modKappa < 0) {
    updateKappa = true;
    ctr->modKappa = 1;
  }
  ctr->shrinkage = as<int>(model["shrinkage"]);   // Store shrinkage
                                                  // 3 = all, 2 = trees, 1 = exposures/interactions, 0 = none
  
  // Fixed effect data
  ctr->Y0 = as<Eigen::VectorXd>(model["Y"]);      // Response variable
  ctr->Ystar = as<Eigen::VectorXd>(model["Y"]);   // For Gaussian default,
  ctr->n = (ctr->Y0).size();                      // Sample size, n

  ctr->Z = as<Eigen::MatrixXd>(model["Z"]);       // Fixed effect design matrix, Z (n x p)
  ctr->Zw = ctr->Z;                               // Copy Z to Zw, will be updated throughout iteration (Z Omega)
  ctr->pZ = (ctr->Z).cols();                      // Number of fixed effect covariates (ncol(Z))

  // Binary component
  ctr->Z1 = as<Eigen::MatrixXd>(model["Z_zi"]);   // (Only for ZINB)Fixed effect design matrix of at-risk model, Z (n x p)
  ctr->Zw1 = ctr->Z1;                             // Copy Z to Zw, will be updated throughout iteration (Z Omega)
  ctr->pZ1 = (ctr->Z1).cols();                    // Number of fixed effect covariates (ncol(Z))

  // (At-risk model) Full conditional initialization for logistic model: V_gamma
  Eigen::MatrixXd VgInv1(ctr->pZ1, ctr->pZ1);     // var-cov matrix
  VgInv1 = (ctr->Z1).transpose() * (ctr->Z1);     // Zt*Z
  VgInv1.diagonal().array() += 1.0 / 100.0;       // Zt*Z + I/c where c = 100
  ctr->Vg1 = VgInv1.inverse();                    // V_gamma
  VgInv1.resize(0,0);                             // Clear VgInv now that we obtained Vg
  ctr->VgChol1 = (ctr->Vg1).llt().matrixL();      // Compute the cholesky of Vg: V_gamma = L*Lt = VgChol * t(VgChol)
                                                  // (Ref: https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html)

  // (Gaussian / Binary / NB) Full conditional initialization for logistic model: V_gamma
  Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ);        // var-cov matrix
  VgInv = (ctr->Z).transpose() * (ctr->Z);        // Zt*Z
  VgInv.diagonal().array() += 1.0 / 100.0;        // Zt*Z + I/c where c = 100
  ctr->Vg = VgInv.inverse();                      // V_gamma
  VgInv.resize(0,0);                              // Clear VgInv now that we obtained Vg
  ctr->VgChol = (ctr->Vg).llt().matrixL();        // Compute the cholesky of Vg: V_gamma = L*Lt = VgChol * t(VgChol)
                                                  // (Ref: https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html)


  // ****** Set up parameters for logistic model ******
  // In the paper, lambda = kappa/omega where kappa = (y_i - n_i/2)
  ctr->binomialSize.resize(ctr->n);               ctr->binomialSize.setZero(); // n of Binom(n, p): Zero vector (n x 1)
  ctr->kappa.resize(ctr->n);                      ctr->kappa.setOnes();        // Kappa: One vector (n x 1)
  ctr->Omega.resize(ctr->n);                      ctr->Omega.setOnes();        // Omega: One vector (n x 1)
  if (ctr->binomial) { // If binomial,
    ctr->binomialSize = as<Eigen::VectorXd>(model["binomialSize"]);   // Save the binomialSize parameter
    ctr->kappa = ctr->Y0 - 0.5 * (ctr->binomialSize);                 // Kappa = Ystar = y - n_i/2
    ctr->Ystar = ctr->kappa;                                          // Save the current Kappa as ystar
  }                                                                  

  // ****** Set up parameters for ZINB model ******
  // Initialize values
  ctr->r = 5;   // Dispersion parameter: Starting at the center of Unif(0, 10)

  // Useful vectors for PG variable sampling 
  ctr->ones.resize(ctr->n);   ctr->ones.setOnes();        // Vector of ones
  ctr->rVec = (ctr->r) * (ctr->ones).array();             // Dispersion parameter as a vector: rep(r, n)

  // If ZINB, calculate the fixed values
  if(ctr->zinb){                                      
    ctr->z2 = 0.5*((ctr->Y0).array() - ctr->r).array();   // Compute z2
    ctr->Ystar = ctr->z2;                                 // Ystar in tdlmm_Cpp.cpp, z2 in modelEst.cpp
  }

  // Initialize parameters
  ctr->w.resize(ctr->n);                                          // At-risk binary latent variable
  ctr->b1 = as<Eigen::VectorXd>(rnorm(ctr->pZ1, 0, sqrt(100)));   // Prior sampling of coefficients for binary component 
  ctr->b2 = as<Eigen::VectorXd>(rnorm(ctr->pZ, 0, sqrt(100)));    // Prior sampling of coefficients for NB component

  ctr->omega1.resize(ctr->n);        ctr->omega1.setOnes();       // Initiate omega1 (binary) as an identity matrix
  ctr->omega2.resize(ctr->n);        ctr->omega2.setOnes();       // Initiate omega2 (NB) as an identity matrix
  ctr->z1.resize(ctr->n);            ctr->z1.setZero();           // z1 for binary component

  // Store the indices of y == 0 (yZeroIdx) & y != 0 (atRiskIdx)
  for(int j = 0; j < ctr->n; j++){
    if((ctr->Y0)[j] == 0){
      (ctr->yZeroIdx).push_back(j);     // y == 0: Could be structural or at-risk zero
      (ctr->w)[j] = 0.5;
    } else {
      (ctr->atRiskIdx).push_back(j);    // y != 0: Determined to be at-risk as y != 0
      (ctr->w)[j] = 1;
    }
  }

  // Fixed effect matrix with [w == 0] zeroed out
  ctr->Zstar = (ctr->Z).array().colwise() * (ctr->w).array();    

  // Calculate the size: yZeroN + nStar should equal to n
  ctr->yZeroN = (ctr->yZeroIdx).size();   // Fixed as we are just counting y == 0
  ctr->nStar = (ctr->atRiskIdx).size();   // Random as some of y == 0 can be at-risk

  // ****** Create exposure data management ******
  std::vector<exposureDat*> Exp;                      // exposureDat object
  Rcpp::List exp_dat = as<Rcpp::List>(model["X"]);    // List of exposures
  ctr->nExp = exp_dat.size();                         // Number of exposures
  for (int i = 0; i < ctr->nExp; ++i) { // For each exposure,
    if (ctr->binomial || ctr->zinb) // For logistic / ZINB,
      Exp.push_back(
        new exposureDat(
          as<Eigen::MatrixXd>(
            as<Rcpp::List>(exp_dat[i])["Tcalc"]))); // Store Tcalc
    else // Gaussian
      Exp.push_back(
        new exposureDat(
          as<Eigen::MatrixXd>(
            as<Rcpp::List>(exp_dat[i])["Tcalc"]), ctr->Z, ctr->Vg)); // Store Tcalc, Z, Vg
  }

  // ****** Mixture/interaction management ******
  ctr->pX = Exp[0]->pX;                                     // Total time lag: T
  ctr->nSplits = 0;
  ctr->interaction = as<int>(model["interaction"]);         // Interaction: TDLMMadd(1) / TDLMMns(2) / TDLMMall(3)
  ctr->nMix = 0;                                            // Default: No interaction terms
  if (ctr->interaction) { // TDLMMns
    ctr->nMix += int (ctr->nExp * (ctr->nExp - 1.0) / 2.0); // (M choose 2) for TDLMMns
    if (ctr->interaction == 2) { // TDLMMall
      ctr->nMix += ctr->nExp;                               // (M choose 2) + M for TDLMMall
    }
  }

  // ****** Create trees ******
  int t;                                                             // Tree index (a = 1, ..., A in the papers)
  (ctr->tree1Exp).resize(ctr->nTrees);    (ctr->tree1Exp).setZero(); // Trees1 vector (Number of trees for Ensemble x 1)
  (ctr->tree2Exp).resize(ctr->nTrees);    (ctr->tree2Exp).setZero(); // Trees2 vector (Number of trees for Ensemble x 1)
  ctr->expProb = as<Eigen::VectorXd>(model["expProb"]);              // Uniform probability for exposure selection
  (ctr->expCount).resize((ctr->expProb).size());                     // Store number of exposures
  (ctr->expInf).resize((ctr->expProb).size());                       

  // Create root nodes to start trees
  std::vector<Node*> trees1; // (Root Node, ..., Root Node)
  std::vector<Node*> trees2; // (Root Node, ..., Root Node)
  NodeStruct *ns;            // Construct NodeStruct object
  ns = new DLNMStruct(0,                          // xmin = 0
                      ctr->nSplits + 1,           // xmax = 1
                      1,                          // tmin = 1
                      int (ctr->pX),              // tmax = pX (total time lag T)
                      as<Eigen::VectorXd>(model["splitProb"]), // Split probability
                      as<Eigen::VectorXd>(model["timeProb"])); // Uniform probability for time split

  for (t = 0; t < ctr->nTrees; ++t) { // For each tree pair a in the ensemble,
    if(ctr->swapStep){
      ctr->tree1Exp(t) = sampleInt(ctr->expProb); // Randomly assign one exposure to Tree1
      ctr->tree2Exp(t) = sampleInt(ctr->expProb); // Randomly assign one expousre to Tree2
    } else {
      if(ctr->interaction == 1){
        // Only works for two exposures now.
        ctr->tree1Exp(t) = 0;
        ctr->tree2Exp(t) = 1;
      } else if (ctr->interaction == 2){
        if(t % 3 == 0){
          ctr->tree1Exp(t) = 0; 
          ctr->tree2Exp(t) = 0; 
        } else if (t % 3 == 1) {
          ctr->tree1Exp(t) = 1; 
          ctr->tree2Exp(t) = 1; 
        } else {
          ctr->tree1Exp(t) = 0; 
          ctr->tree2Exp(t) = 1; 
        }
      }
    }
    
    trees1.push_back(new Node(0, 1));           // Fill in the tree vector 1 with a node of 0 depth and not terminal node
    trees2.push_back(new Node(0, 1));           // Fill in the tree vector 2 with a node of 0 depth and not terminal node
    trees1[t]->nodestruct = ns->clone();        // Tree 1: Copy NodeStruct object to Node parameter
    trees2[t]->nodestruct = ns->clone();        // Tree 2: Copy NodeStruct object to Node parameter
    Exp[ctr->tree1Exp(t)]->updateNodeVals(trees1[t]); // Update trees in tree vector 1
    Exp[ctr->tree2Exp(t)]->updateNodeVals(trees2[t]); // Update trees in tree vector 2
  }
  delete ns;                                    // Delete NodeStruct ns as we cloned it
  
  // ****** Setup model logs ******
  tdlmLog *dgn = new tdlmLog;                                                   // Log object
  (dgn->gamma).resize(ctr->pZ, ctr->nRec);          (dgn->gamma).setZero();     // Fixed effect coefficient
  (dgn->sigma2).resize(ctr->nRec);                  (dgn->sigma2).setZero();    // Gaussian variance
  (dgn->kappa).resize(ctr->nRec);                   (dgn->kappa).setZero();     // Logistic z1
  (dgn->nu).resize(ctr->nRec);                      (dgn->nu).setZero();        // Shrinkage parameter 1
  (dgn->tau).resize(ctr->nTrees, ctr->nRec);        (dgn->tau).setZero();       // Shrinkage parameter 2
  (dgn->muExp).resize(ctr->nExp, ctr->nRec);        (dgn->muExp).setZero();     // Exposure-specific variance
  if (ctr->interaction > 0) { // For TDLMMns / TDLMMall
    (dgn->muMix).resize(ctr->nMix, ctr->nRec);      (dgn->muMix).setZero();
    (dgn->mixInf).resize(ctr->nMix, ctr->nRec);     (dgn->mixInf).setZero();
    (dgn->mixCount).resize(ctr->nMix, ctr->nRec);   (dgn->muMix).setZero();
  } else { // For TDLMMadd -> No interaction terms
    (dgn->muMix).resize(1, 1);                      (dgn->muMix).setZero();
    (dgn->mixInf).resize(1, 1);                     (dgn->mixInf).setZero();
    (dgn->mixCount).resize(1, 1);                   (dgn->mixCount).setZero();
  }
  (dgn->expProb).resize(ctr->nExp, ctr->nRec);      (dgn->expProb).setZero();     // Exposure selection probability
  (dgn->expCount).resize(ctr->nExp, ctr->nRec);     (dgn->expCount).setZero();    // Number of exposures
  (dgn->expInf).resize(ctr->nExp, ctr->nRec);       (dgn->expInf).setZero();
  (dgn->fhat).resize(ctr->n);                       (dgn->fhat).setZero();        // fhat: Sum of DLM effect except tree pair a
  (dgn->termNodes).resize(ctr->nTrees, ctr->nRec);  (dgn->termNodes).setZero();   // Terminal nodes of tree vector 1
  (dgn->termNodes2).resize(ctr->nTrees, ctr->nRec); (dgn->termNodes2).setZero();  // Terminal nodes of tree vector 2
  (dgn->tree1Exp).resize(ctr->nTrees, ctr->nRec);   (dgn->tree1Exp).setZero();    // Exposure vector of tree vector 1
  (dgn->tree2Exp).resize(ctr->nTrees, ctr->nRec);   (dgn->tree2Exp).setZero();    // Exposure vector of tree vector 2

  // ZINB log
  (dgn->b1).resize(ctr->pZ1, ctr->nRec);             (dgn->b1).setZero();         // Binary component coeffient (p x MCMC)
  (dgn->b2).resize(ctr->pZ, ctr->nRec);              (dgn->b2).setZero();         // NB component coeffient (p x MCMC)
  (dgn->r).resize(ctr->nRec);                        (dgn->r).setZero();          // Dispersion parameter
  (dgn->wMat).resize(ctr->n, ctr->nRec);             (dgn->wMat).setZero();       // At-risk Auxiliary: Each column is w at each iteration

  // ****** Initial draws ******
  (ctr->fhat).resize(ctr->n);                         (ctr->fhat).setZero();          // fhat: A vector of zeros
  ctr->R = ctr->Ystar;                                                                // Partial residual
  (ctr->gamma).resize(ctr->pZ);                                                       // Logistic model coefficient
  (ctr->totTermExp).resize(ctr->nExp);                (ctr->totTermExp).setZero();    
  (ctr->sumTermT2Exp).resize(ctr->nExp);              (ctr->sumTermT2Exp).setZero();  
  (ctr->muExp).resize(ctr->nExp);                     (ctr->muExp).setOnes();                             // Exposure-specific variance

  // If there is an interaction, define more parameters
  if (ctr->interaction) { 
    (ctr->totTermMix).resize(ctr->nExp, ctr->nExp);   (ctr->totTermMix).setZero();    // Interaction parameters
    (ctr->sumTermT2Mix).resize(ctr->nExp, ctr->nExp); (ctr->sumTermT2Mix).setZero();  //
    (ctr->muMix).resize(ctr->nExp, ctr->nExp);        (ctr->muMix).setOnes();         // interaction-specific variance
    (ctr->mixInf).resize(ctr->nExp, ctr->nExp);       (ctr->mixInf).setZero();        // 
    (ctr->mixCount).resize(ctr->nExp, ctr->nExp);     (ctr->mixCount).setZero();      
  }

  // Hyperparameter initial values
  ctr->totTerm = 0;       // This is represented by "b" in the paper
  ctr->sumTermT2 = 0;
  ctr->nu = 1.0; 
  ctr->sigma2 = 1.0;

  // model estimation with the initial values
  tdlmModelEst(ctr); 

  // node-specific full-conditional with half-Cauchy: C+(0, 1)
  rHalfCauchyFC(&(ctr->nu), ctr->nTrees, 0.0);
  (ctr->tau).resize(ctr->nTrees);                 (ctr->tau).setOnes();
  for (t = 0; t < ctr->nTrees; ++t) {
    if (ctr->shrinkage > 1) // shrinkage on all (tree & exposure) or just tree
      rHalfCauchyFC(&(ctr->tau(t)), 0.0, 0.0);
  }
  ctr->nTerm.resize(ctr->nTrees);                 (ctr->nTerm).setOnes();
  ctr->nTerm2.resize(ctr->nTrees);                (ctr->nTerm2).setOnes();
  (ctr->Rmat).resize(ctr->n, ctr->nTrees);        (ctr->Rmat).setZero();

  // ****** Create Progress Meter ******
  progressMeter* prog = new progressMeter(ctr);

  // ****** Beginning of MCMC ******
  double sigmanu;
  std::size_t s;
  // For all the MCMC iteration + burn-in
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) { 

    // Record the iteration number if the iteration is after the burn-in and not one of thinning iteration
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) { 
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin); 
    } else {
      ctr->record = 0;
    }

    // [Update trees]
    // Rcout << "Update MCMC trees \n";
    // Add the first column of Rmat to the partial residual
    ctr->R += (ctr->Rmat).col(0); 

    // Reset the parameters
    (ctr->fhat).setZero();                
    ctr->totTerm = 0;                    
    ctr->sumTermT2 = 0;                   
    (ctr->totTermExp).setZero();
    (ctr->sumTermT2Exp).setZero();
    (ctr->expCount).setZero();
    (ctr->mixCount).setZero();
    (ctr->expInf).setZero();
    (ctr->mixInf).setZero();
    if (ctr->interaction > 0) {
      (ctr->totTermMix).setZero();                (ctr->sumTermT2Mix).setZero();
    }

    // For each tree pair, perform one iteration of MCMC
    for (t = 0; t < ctr->nTrees; ++t) {                       // t = index for tree pair
      tdlmmTreeMCMC(t, trees1[t], trees2[t], ctr, dgn, Exp);  // Tree transitions
      ctr->fhat += (ctr->Rmat).col(t);
      if (t < ctr->nTrees - 1) 
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t); 
    }
    // Rcout << " end trees \n";

    // * Pre-calculations for control and variance
    ctr->R = ctr->Ystar - ctr->fhat;             
    ctr->sumTermT2 = (ctr->sumTermT2Exp).sum();
    ctr->totTerm = (ctr->totTermExp).sum();
    if(ctr->interaction) {
      ctr->sumTermT2 += (ctr->sumTermT2Mix).sum();
      ctr->totTerm += (ctr->totTermMix).sum();
    }

    // * Update model (Binomial model should go through MCMC as well)
    tdlmModelEst(ctr); // Go to modelEst.cpp::tdlmModelEst
    
    // * Update variance parameters from supplemental material -> Horseshoe sampling with IG & Half-cauchy relationship
    rHalfCauchyFC(&(ctr->nu), ctr->totTerm, ctr->sumTermT2 / ctr->sigma2);  
    // if ((ctr->nu) != (ctr->nu)) 
    //   stop("\nNaN values (nu) occured during model run, rerun model.\n");
    sigmanu = ctr->sigma2 * ctr->nu;
    if ((ctr->shrinkage == 3) || (ctr->shrinkage == 1)) { // Exposure shrinkage
      for (int i = 0; i < ctr->nExp; ++i) {
        rHalfCauchyFC(&(ctr->muExp(i)), ctr->totTermExp(i), ctr->sumTermT2Exp(i) / sigmanu);
        // if ((ctr->muExp)(i) != (ctr->muExp)(i)) 
        //   stop("\nNaN values (muExp) occured during model run, rerun model.\n");
        if (ctr->interaction) { // mu_{Sa1Sa2}
          for (int j = i; j < ctr->nExp; ++j) {
            if ((j > i) || (ctr->interaction == 2))
              rHalfCauchyFC(&(ctr->muMix(j, i)), ctr->totTermMix(j, i),
                            ctr->sumTermT2Mix(j, i) / sigmanu);
            // if ((ctr->muMix)(j,i) != (ctr->muMix)(j,i)) 
            //   stop("\nNaN values (muMix) occured during model run, rerun model.\n");
          } // end for loop updating interaction variances
        } // end if interactions
      } // end for loop updating exposure variances
    } // end if shrinkage == 3 or 1R

    // * Update exposure selection probability from supplemental material
    if ((ctr->b > 1000) || (ctr->b > (0.5 * ctr->burn))) {
      if (updateKappa) {
        // double modKappaNew = exp(log(ctr->modKappa) + R::rnorm(0, 0.5));
        double modKappaNew = R::rgamma(1.0, ctr->nTrees/4.0);
        double mhrDir =
          logDirichletDensity(ctr->expProb,
                              ((ctr->expCount).array() + 
                               modKappaNew).matrix()) -
          logDirichletDensity(ctr->expProb,
                              ((ctr->expCount).array() + 
                               ctr->modKappa).matrix());

        if (log(R::runif(0, 1)) < mhrDir)
          ctr->modKappa = modKappaNew;
      }
      
      ctr->expProb = rDirichlet(((ctr->expCount).array() + ctr->modKappa).matrix());
    }

    // * Record
    if (ctr->record > 0) {
      dgn->fhat += ctr->fhat;
      (dgn->gamma).col(ctr->record - 1) = ctr->gamma;
      (dgn->sigma2)(ctr->record - 1) = ctr->sigma2;
      (dgn->nu)(ctr->record - 1) = ctr->nu;
      (dgn->tau).col(ctr->record - 1) = ctr->tau;
      (dgn->termNodes).col(ctr->record - 1) = ctr->nTerm;
      (dgn->termNodes2).col(ctr->record - 1) = ctr->nTerm2;
      (dgn->tree1Exp).col(ctr->record - 1) = ctr->tree1Exp;
      (dgn->tree2Exp).col(ctr->record - 1) = ctr->tree2Exp;
      (dgn->expProb).col(ctr->record - 1) = ctr->expProb;
      (dgn->expCount).col(ctr->record - 1) = ctr->expCount;
      (dgn->expInf).col(ctr->record - 1) = ctr->expInf;
      (dgn->muExp).col(ctr->record - 1) = ctr->muExp;
      (dgn->kappa)(ctr->record - 1) = ctr->modKappa;

      // ZINB
      (dgn->b1).col(ctr->record - 1) = ctr->b1;
      (dgn->b2).col(ctr->record - 1) = ctr->b2;
      (dgn->r)(ctr->record - 1) = ctr->r;
      (dgn->wMat).col(ctr->record - 1) = ctr->w;

      if (ctr->interaction) {
        int k = 0;
        for (int i = 0; i < ctr->nExp; ++i) {
          for (int j = i; j < ctr->nExp; ++j) {
            if ((j > i) || (ctr->interaction == 2)) {
              dgn->muMix(k, ctr->record - 1) = ctr->muMix(j, i);
              dgn->mixInf(k, ctr->record - 1) = ctr->mixInf(j, i);
              dgn->mixCount(k, ctr->record - 1) = ctr->mixCount(j, i);
              ++k;
            }
          }
        }
      }
    }


    if (ctr->diagnostics){
      Rcout << "MCMC Iteration : ";
      Rcout << ctr->b;
      Rcout << "\n";

      Rcout << "Sigma2 : ";
      Rcout << ctr->sigma2;
      Rcout << "\n";
      Rcout << "nu : ";
      Rcout << ctr->nu;
      Rcout << "\n";

      for(int p = 0; p < ctr->nTrees; p++){
        Rcout << ctr->tau[p];
        Rcout << "|";
      }
      Rcout << "\n";

      Rcout << "muExp(0) : ";
      Rcout << ctr->muExp(0);
      Rcout << "\n";
      Rcout << "muExp(1) : ";
      Rcout << ctr->muExp(1);
      Rcout << "\n";

      if(ctr->interaction > 0){
        Rcout << "muMix(0-1) : ";
        Rcout << ctr->muMix(1, 0);
        Rcout << "\n";
        Rcout << "muMix(1-2) : ";
        Rcout << ctr->muMix(2, 1);
        Rcout << "\n";
        Rcout << "muMix(0-2) : ";
        Rcout << ctr->muMix(2, 0);
        Rcout << "\n";
      }

      Rcout << "--------------------------------\n";
    }

    // * Progress
    prog->printMark();
  } // End MCMC -----------------------------------------------------------------


  // * Setup data for return
  Eigen::MatrixXd DLM((dgn->DLMexp).size(), 8);
  for (s = 0; s < (dgn->DLMexp).size(); ++s)
    DLM.row(s) = dgn->DLMexp[s];
  Eigen::VectorXd sigma2 = dgn->sigma2;
  Eigen::VectorXd nu = dgn->nu;
  Eigen::VectorXd kappa = dgn->kappa;
  Eigen::VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  Eigen::MatrixXd gamma = (dgn->gamma).transpose();
  Eigen::MatrixXd tau = (dgn->tau).transpose();
  Eigen::MatrixXd termNodes = (dgn->termNodes).transpose();
  Eigen::MatrixXd termNodes2 = (dgn->termNodes2).transpose();
  Eigen::MatrixXd tree1Exp = (dgn->tree1Exp).transpose();
  Eigen::MatrixXd tree2Exp = (dgn->tree2Exp).transpose();
  Eigen::MatrixXd expProb = (dgn->expProb).transpose();
  Eigen::MatrixXd expCount = (dgn->expCount).transpose();
  Eigen::MatrixXd expInf = (dgn->expInf).transpose();
  Eigen::MatrixXd mixInf = (dgn->mixInf).transpose();
  Eigen::MatrixXd muExp = (dgn->muExp).transpose();
  Eigen::MatrixXd mixCount = (dgn->mixCount).transpose();
  Eigen::MatrixXd muMix(1, 1); muMix.setZero();
  Eigen::MatrixXd MIX(0, 10); MIX.setZero();

  // ZINB return 
  Eigen::MatrixXd b1 = (dgn->b1).transpose(); // binary component coefficients
  Eigen::MatrixXd b2 = (dgn->b2).transpose(); // count component coefficients (This gets return for NB)
  Eigen::VectorXd r = dgn->r;                 // dispersion parameter
  Eigen::MatrixXd wMat = dgn->wMat;           // What percentage of iteration was an individual at risk?

  if (ctr->interaction) {
    muMix.resize((dgn->muMix).cols(), (dgn->muMix).rows());
    muMix = (dgn->muMix).transpose();
    MIX.resize((dgn->MIXexp).size(), 10);
    for (s = 0; s < (dgn->MIXexp).size(); ++s)
      MIX.row(s) = dgn->MIXexp[s];
  }
  Eigen::MatrixXd Accept((dgn->TreeAccept).size(), 7);
  for (s = 0; s < (dgn->TreeAccept).size(); ++s)
    Accept.row(s) = dgn->TreeAccept[s];
  delete prog;
  delete ctr;
  delete dgn;
  for (s = 0; s < Exp.size(); ++s)
    delete Exp[s];
  for (s = 0; s < trees1.size(); ++s) {
    delete trees1[s];
    delete trees2[s];
  }

  return(Rcpp::List::create(Named("DLM") = wrap(DLM),
                            Named("MIX") = wrap(MIX),
                            Named("gamma") = wrap(gamma),
                            // Named("fhat") = wrap(fhat),
                            Named("sigma2") = wrap(sigma2),
                            //Named("nu") = wrap(nu),
                            //Named("tau") = wrap(tau),
                            Named("termNodes") = wrap(termNodes),
                            Named("termNodes2") = wrap(termNodes2), 
                            // Named("tree1Exp") = wrap(tree1Exp),
                            // Named("tree2Exp") = wrap(tree2Exp),
                            Named("expProb") = wrap(expProb),
                            Named("expInf") = wrap(expInf),
                            Named("expCount") = wrap(expCount),
                            Named("mixInf") = wrap(mixInf),
                            Named("mixCount") = wrap(mixCount), 
                            Named("muExp") = wrap(muExp),
                            Named("muMix") = wrap(muMix),
                            //Named("kappa") = wrap(kappa),
                            Named("treeAccept") = wrap(Accept),
                            Named("b1") = wrap(b1),
                            Named("b2") = wrap(b2),
                            Named("r") = wrap(r),
                            Named("wMat") = wrap(wMat)));
} // end tdlmm_Cpp






