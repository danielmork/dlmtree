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
#include <RcppEigen.h>
#include "modelCtr.h"
#include "exposureDat.h"
#include "Node.h"
#include "NodeStruct.h"
#include "Fncs.h"
using namespace Rcpp;

/**
 * @brief 
 * 
 * @param nodes1 
 * @param nodes2 
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
  treeMHR out;
  int pX1 = nodes1.size();
  int pX2 = nodes2.size();
  int pXd = pX1 + pX2;
  int interaction = 0;
  if (mixVar != 0) {
    pXd += pX1 * pX2;
    interaction = 1;
  }
  // Rcout << ".";
  out.Xd.resize(ctr->n, pXd);         out.Xd.setZero();
  Eigen::MatrixXd ZtX(ctr->pZ, pXd);  ZtX.setZero();
  Eigen::VectorXd diagVar(pXd);       diagVar.setZero();

  int i, j, k;
  for (i = 0; i < pX1; ++i) {
    out.Xd.col(i) = (nodes1[i]->nodevals)->X;
    diagVar(i) = 1.0 / (m1Var * treeVar);
    if (ctr->binomial)
      ZtX.col(i) = ctr->Zw.transpose() * (nodes1[i]->nodevals)->X;
    else 
      ZtX.col(i) = (nodes1[i]->nodevals)->ZtX;
  }

  for (j = 0; j < pX2; ++j) {
    k = pX1 + j;
    out.Xd.col(k) = (nodes2[j]->nodevals)->X;
    diagVar(k) = 1.0 / (m2Var * treeVar);
    if (ctr->binomial)
      ZtX.col(k) = ctr->Zw.transpose() * (nodes2[j]->nodevals)->X;
    else 
      ZtX.col(k) = (nodes2[j]->nodevals)->ZtX;
  }

  if (interaction) {
    for (i = 0; i < pX1; ++i) {
      for (j = 0; j < pX2; ++j) {
        k = pX1 + pX2 + i * pX2 + j;
        out.Xd.col(k) =
          (((nodes1[i]->nodevals)->X).array() *
          ((nodes2[j]->nodevals)->X).array()).matrix();
        diagVar(k) = 1.0 / (mixVar * treeVar);
        ZtX.col(k) = ctr->Zw.transpose() * out.Xd.col(k);
      }
    }
  }

  // Rcout << ".";
  // calculate MHR
  const Eigen::MatrixXd VgZtX = ctr->Vg * ZtX;
  Eigen::MatrixXd tempV(pXd, pXd);
  Eigen::VectorXd XtVzInvR(pXd);
  if (ctr->binomial) {
    tempV = out.Xd.transpose() * 
      (out.Xd.array().colwise() * ctr->Omega.array()).matrix();
    tempV.noalias() -= ZtX.transpose() * VgZtX;
    XtVzInvR = (out.Xd.array().colwise() * 
                ctr->Omega.array()).matrix().transpose() * ctr->R;
    
  } else {
    if (newTree) {
      tempV.triangularView<Eigen::Lower>() = out.Xd.transpose() * out.Xd;
      tempV.noalias() -= ZtX.transpose() * VgZtX;
      out.tempV = tempV;
    } else {
      tempV = tree->nodevals->tempV;
    }
    XtVzInvR = out.Xd.transpose() * ctr->R;
  }
  XtVzInvR.noalias() -= VgZtX.transpose() * ZtR;
  tempV.diagonal().noalias() += diagVar;

  // Rcout << ".";
  Eigen::MatrixXd VTheta(pXd, pXd);
  VTheta.triangularView<Eigen::Lower>() =
    tempV.selfadjointView<Eigen::Lower>().llt().solve(
        Eigen::MatrixXd::Identity(pXd, pXd));
  const Eigen::MatrixXd VThetaChol =
    VTheta.selfadjointView<Eigen::Lower>().llt().matrixL();
  const Eigen::VectorXd ThetaHat =
    VTheta.selfadjointView<Eigen::Lower>() * XtVzInvR;


  // Rcout << ".";
  Eigen::VectorXd ThetaDraw = ThetaHat;
  ThetaDraw.noalias() +=
    VThetaChol * as<Eigen::VectorXd>(rnorm(pXd, 0, sqrt(ctr->sigma2)));

  out.drawAll = ThetaDraw;
  out.draw1 = ThetaDraw.head(pX1);
  out.term1T2 = (out.draw1).dot(out.draw1);
  out.nTerm1 = double(pX1);
  out.draw2 = ThetaDraw.segment(pX1, pX2);
  out.term2T2 = (out.draw2).dot(out.draw2);
  out.nTerm2 = double(pX2);
  if (interaction) {
    out.drawMix = ThetaDraw.tail(pXd - pX1 - pX2);
    out.mixT2 = (out.drawMix).dot(out.drawMix);
  }
  out.beta = ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();

  out.pXd = pXd;
  return(out);
}

/**
 * @brief 
 * 
 * @param t 
 * @param tree1 
 * @param tree2 
 * @param ctr 
 * @param dgn 
 * @param Exp 
 */
void tdlmmTreeMCMC(int t, Node *tree1, Node *tree2, tdlmCtr *ctr, tdlmLog *dgn,
                   std::vector<exposureDat*> Exp)
{

  int m1, m2, newExp, success, step1, step2;
  double stepMhr, ratio;
  double m1Var, m2Var, mixVar, newExpVar, newMixVar, treeVar;
  double RtR = -1.0;
  double RtZVgZtR = 0;
  std::vector<Node*> term1, term2, newTerm;
  Node* newTree = 0;
  treeMHR mhr0, mhr;

  // Rcout << ".";
  term1 = tree1->listTerminal();
  term2 = tree2->listTerminal();
  treeVar = (ctr->nu) * (ctr->tau[t]);
  m1 = ctr->tree1Exp[t];
  m2 = ctr->tree2Exp[t];
  m1Var = ctr->muExp(m1);
  m2Var = ctr->muExp(m2);
  mixVar = 0;
  if ((ctr->interaction) && ((ctr->interaction == 2) || (m1 != m2))) {
    if (m1 <= m2)
      mixVar = ctr->muMix(m2, m1);
    else
      mixVar = ctr->muMix(m1, m2);
  }
  // Rcout << ".";
  Eigen::VectorXd ZtR = (ctr->Zw).transpose() * (ctr->R);

  // * Update tree 1
  newExp    = m1;
  newExpVar = m1Var;
  newMixVar = mixVar;
  stepMhr   = 0;
  success   = 0;

  // Rcout << ".";
  // * List current tree terminal nodes
  step1 = sampleInt(ctr->stepProb, 1);
  if ((term1.size() == 1) && (step1 < 3))
    step1 = 0;

  // * Propose update
  if (step1 < 3) {
    // Rcout << "g" << step1;
    stepMhr = tdlmProposeTree(tree1, Exp[m1], ctr, step1);
    // Rcout << ".";
    success = tree1->isProposed();
    newTerm = tree1->listTerminal(1);

  // * Switch exposures
  } else {
    // Rcout << "s";
    newExp = sampleInt(ctr->expProb);
    if (newExp != m1) {
      success   = 1;
      newExpVar = ctr->muExp(newExp);
      newTree   = new Node(*tree1);
      newTree->setUpdate(1);
      newTerm   = newTree->listTerminal();
      for (Node* nt : newTerm)
        Exp[newExp]->updateNodeVals(nt);


      if ((ctr->interaction) && ((ctr->interaction == 2) || (newExp != m2))) {
        if (newExp <= m2)
          newMixVar = ctr->muMix(m2, newExp);
        else
          newMixVar = ctr->muMix(newExp, m2);
      } else {
        newMixVar = 0;
      }
    }
  }

  // * Tree 1 MHR
  // Rcout << "a";
  if ((tree1->nodevals->tempV).rows() == 0)
    mhr0 = mixMHR(term1, term2, ctr, ZtR, treeVar, 
                  m1Var, m2Var, mixVar, tree1, 1);
  else
    mhr0 = mixMHR(term1, term2, ctr, ZtR, treeVar, 
                  m1Var, m2Var, mixVar, tree1, 0);

  if (success) {
    // calculate new tree part of MHR and draw node effects
    // Rcout << "b";
    mhr = mixMHR(newTerm, term2, ctr, ZtR, treeVar, 
                 newExpVar, m2Var, newMixVar, tree1, 1);
    // combine mhr parts into log-MH ratio
    if (ctr->binomial) {
      ratio = stepMhr + mhr.logVThetaChol - mhr0.logVThetaChol +
        0.5 * (mhr.beta - mhr0.beta) -
        (0.5 * ((log(treeVar * newExpVar) * mhr.nTerm1) -
         (log(treeVar * m1Var) * mhr0.nTerm1)));
    } else {
      if (RtR < 0) {
        RtR = (ctr->R).dot(ctr->R);
        RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);
      }
      ratio = stepMhr + mhr.logVThetaChol - mhr0.logVThetaChol -
        (0.5 * (ctr->n + 1.0) *
         (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
          log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
        (0.5 * ((log(treeVar * newExpVar) * mhr.nTerm1) -
         (log(treeVar * m1Var) * mhr0.nTerm1)));
    }

    if (newMixVar != 0)
      ratio -= 0.5 * log(treeVar * newMixVar) * mhr.nTerm1 * mhr0.nTerm2;
    if (mixVar != 0)
      ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2;

    if (log(R::runif(0, 1)) < ratio) {
      mhr0 = mhr;
      success = 2;

      if (step1 == 3) {
        m1      = newExp;
        m1Var   = newExpVar;
        mixVar  = newMixVar;
        tree1->replaceNodeVals(newTree);
      } else {
        tree1->accept();
      }
      if (!(ctr->binomial)) {
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

  // * Record tree 1
  if (ctr->diagnostics) {
    Eigen::VectorXd acc(7);
    acc << 1, step1, success, m1, term1.size(), stepMhr, ratio;
    (dgn->TreeAccept).push_back(acc);
  }

  // * Update tree 2
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
    // Rcout << ".";
    success = tree2->isProposed();
    newTerm = tree2->listTerminal(1);

  // * Switch exposures
  } else {
    // Rcout << "s";
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
    if (ctr->binomial) {
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
      if (!(ctr->binomial)) {
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


  // * Record tree 2
  if (ctr->diagnostics) {
    Eigen::VectorXd acc(7);
    acc << 1, step2, success, m2, term2.size(), stepMhr, ratio;
    (dgn->TreeAccept).push_back(acc);
  }


  // * Update variance and residuals
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

  (ctr->Rmat).col(t) = mhr0.Xd * mhr0.drawAll;

  // * Record
  if (ctr->record > 0) {
    Eigen::VectorXd rec(7);
    Eigen::VectorXd mix(10);
    rec << ctr->record, t, 0, 0, 0, 0, 0;
    mix << ctr->record, t, 0, 0, 0, 0, 0, 0, 0, 0;//(ctr->tau)(t) *
    int k = 0;
    for(int i = 0; i < mhr0.nTerm1; ++i) {
      rec[2] = m1;
      rec[3] = (term1[i]->nodestruct)->get(3);
      rec[4] = (term1[i]->nodestruct)->get(4);
      rec[5] = mhr0.draw1(i);
      rec[6] = (ctr->tau)(t) * m1Var;
      (dgn->DLMexp).push_back(rec);
      for (int j = 0; j < mhr0.nTerm2; ++j) {
        if (i == 0) {
          rec[2] = m2;
          rec[3] = (term2[j]->nodestruct)->get(3);
          rec[4] = (term2[j]->nodestruct)->get(4);
          rec[5] = mhr0.draw2(j);
          rec[6] = (ctr->tau)(t) * m2Var;
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
  // * Set up model control
  tdlmCtr *ctr = new tdlmCtr;
  ctr->iter = as<int>(model["nIter"]);
  ctr->burn = as<int>(model["nBurn"]);
  ctr->thin = as<int>(model["nThin"]);
  ctr->nRec = floor(ctr->iter / ctr->thin);
  ctr->nTrees = as<int>(model["nTrees"]);
  ctr->verbose = as<bool>(model["verbose"]);
  ctr->diagnostics = as<bool>(model["diagnostics"]);
  ctr->binomial = as<bool>(model["binomial"]);
  ctr->stepProb = as<std::vector<double> >(model["stepProb"]);
  ctr->treePrior = as<std::vector<double> >(model["treePrior"]);
  ctr->modKappa = as<double>(model["mixPrior"]);
  bool updateKappa = false;
  if (ctr->modKappa < 0) {
    updateKappa = true;
    ctr->modKappa = 1;
  }
  ctr->shrinkage = as<int>(model["shrinkage"]); 
  // 3 = all, 2 = trees, 1 = exposures/interactions, 0 = none
  
  // * Set up model data
  ctr->Y = as<Eigen::VectorXd>(model["Y"]);
  ctr->n = (ctr->Y).size();
  ctr->Z = as<Eigen::MatrixXd>(model["Z"]);
  ctr->Zw = ctr->Z;
  ctr->pZ = (ctr->Z).cols();
  Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ);
  VgInv = (ctr->Z).transpose() * (ctr->Z);
  VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg = VgInv.inverse();
  VgInv.resize(0,0);
  ctr->VgChol = (ctr->Vg).llt().matrixL();
  
  // * Set up data for logistic model
  ctr->binomialSize.resize(ctr->n);               ctr->binomialSize.setZero();
  ctr->kappa.resize(ctr->n);                      ctr->kappa.setOnes();
  ctr->Omega.resize(ctr->n);                      ctr->Omega.setOnes();
  if (ctr->binomial) {
    ctr->binomialSize = as<Eigen::VectorXd>(model["binomialSize"]);
    ctr->kappa = ctr->Y - 0.5 * (ctr->binomialSize);
    ctr->Y = ctr->kappa;
  }

  // * Create exposure data management
  std::vector<exposureDat*> Exp;
  Rcpp::List exp_dat = as<Rcpp::List>(model["X"]);
  ctr->nExp = exp_dat.size();
  for (int i = 0; i < ctr->nExp; ++i) {
    if (ctr->binomial)
      Exp.push_back(
        new exposureDat(
          as<Eigen::MatrixXd>(as<Rcpp::List>(exp_dat[i])["Tcalc"])));
    else
      Exp.push_back(
        new exposureDat(
          as<Eigen::MatrixXd>(
            as<Rcpp::List>(exp_dat[i])["Tcalc"]), ctr->Z, ctr->Vg));
  }
  ctr->pX = Exp[0]->pX;
  ctr->nSplits = 0;
  ctr->interaction = as<int>(model["interaction"]);
  // 0 = no interation, 1 = no-self, 2 = all
  ctr->nMix = 0;
  if (ctr->interaction) {
    ctr->nMix += int (ctr->nExp * (ctr->nExp - 1.0) / 2.0);
    if (ctr->interaction == 2) {
      ctr->nMix += ctr->nExp;
    }
  }
  
  // * Create trees
  int t;
  (ctr->tree1Exp).resize(ctr->nTrees);              (ctr->tree1Exp).setZero();
  (ctr->tree2Exp).resize(ctr->nTrees);              (ctr->tree2Exp).setZero();
  ctr->expProb = as<Eigen::VectorXd>(model["expProb"]);
  (ctr->expCount).resize((ctr->expProb).size());
  (ctr->expInf).resize((ctr->expProb).size());
  std::vector<Node*> trees1;
  std::vector<Node*> trees2;
  NodeStruct *ns;
  ns = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                      as<Eigen::VectorXd>(model["splitProb"]),
                      as<Eigen::VectorXd>(model["timeProb"]));
  for (t = 0; t < ctr->nTrees; ++t) {
    ctr->tree1Exp(t) = sampleInt(ctr->expProb);
    ctr->tree2Exp(t) = sampleInt(ctr->expProb);
    trees1.push_back(new Node(0, 1));
    trees2.push_back(new Node(0, 1));
    trees1[t]->nodestruct = ns->clone();
    trees2[t]->nodestruct = ns->clone();
    Exp[ctr->tree1Exp(t)]->updateNodeVals(trees1[t]);
    Exp[ctr->tree2Exp(t)]->updateNodeVals(trees2[t]);
  }
  delete ns;
  
  // * Setup model logs
  tdlmLog *dgn = new tdlmLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec);          (dgn->gamma).setZero();
  (dgn->sigma2).resize(ctr->nRec);                  (dgn->sigma2).setZero();
  (dgn->kappa).resize(ctr->nRec);                   (dgn->kappa).setZero();
  (dgn->nu).resize(ctr->nRec);                      (dgn->nu).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec);        (dgn->tau).setZero();
  (dgn->muExp).resize(ctr->nExp, ctr->nRec);        (dgn->muExp).setZero();
  if (ctr->interaction > 0) {
    (dgn->muMix).resize(ctr->nMix, ctr->nRec);      (dgn->muMix).setZero();
    (dgn->mixInf).resize(ctr->nMix, ctr->nRec);     (dgn->mixInf).setZero();
    (dgn->mixCount).resize(ctr->nMix, ctr->nRec);   (dgn->muMix).setZero();
  } else {
    (dgn->muMix).resize(1, 1);                      (dgn->muMix).setZero();
    (dgn->mixInf).resize(1, 1);                     (dgn->mixInf).setZero();
    (dgn->mixCount).resize(1, 1);                   (dgn->mixCount).setZero();
  }
  (dgn->expProb).resize(ctr->nExp, ctr->nRec);      (dgn->expProb).setZero();
  (dgn->expCount).resize(ctr->nExp, ctr->nRec);     (dgn->expCount).setZero();
  (dgn->expInf).resize(ctr->nExp, ctr->nRec);       (dgn->expInf).setZero();
  (dgn->fhat).resize(ctr->n);                       (dgn->fhat).setZero();
  (dgn->termNodes).resize(ctr->nTrees, ctr->nRec);  (dgn->termNodes).setZero();
  (dgn->termNodes2).resize(ctr->nTrees, ctr->nRec); (dgn->termNodes2).setZero();
  (dgn->tree1Exp).resize(ctr->nTrees, ctr->nRec);   (dgn->tree1Exp).setZero();
  (dgn->tree2Exp).resize(ctr->nTrees, ctr->nRec);   (dgn->tree2Exp).setZero();

  // * Initial draws
  (ctr->fhat).resize(ctr->n);                     (ctr->fhat).setZero();
  ctr->R = ctr->Y;
  (ctr->gamma).resize(ctr->pZ);
  // Load initial params for faster convergence in binomial model
  if (ctr->binomial) {
    ctr->gamma = as<VectorXd>(model["initParams"]);
    ctr->Omega =rcpp_pgdraw(ctr->binomialSize, ctr->fhat + ctr->Z * ctr->gamma);
    ctr->Zw = ctr->Omega.asDiagonal() * ctr->Z;
    ctr->VgInv =   ctr->Z.transpose() * ctr->Zw;
    ctr->VgInv.diagonal().array() += 1 / 100000.0;
    ctr->Vg = ctr->VgInv.inverse();
    ctr->VgChol = ctr->Vg.llt().matrixL();
    // recalculate 'pseudo-Y' = kappa / omega, kappa = (y - n_b)/2
    ctr->Y = ctr->kappa.array() / ctr->Omega.array();
  }
  (ctr->totTermExp).resize(ctr->nExp);            (ctr->totTermExp).setZero();
  (ctr->sumTermT2Exp).resize(ctr->nExp);          (ctr->sumTermT2Exp).setZero();
  (ctr->muExp).resize(ctr->nExp); (ctr->muExp).setOnes();
  if (ctr->interaction) {
    (ctr->totTermMix).resize(ctr->nExp, ctr->nExp); (ctr->totTermMix).setZero();
    (ctr->sumTermT2Mix).resize(ctr->nExp, ctr->nExp); (ctr->sumTermT2Mix).setZero();
    (ctr->muMix).resize(ctr->nExp, ctr->nExp);    (ctr->muMix).setOnes();
    (ctr->mixInf).resize(ctr->nExp, ctr->nExp);   (ctr->mixInf).setZero();
    (ctr->mixCount).resize(ctr->nExp, ctr->nExp); (ctr->mixCount).setZero();
  }
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  ctr->nu = 1.0; // ! Need to define nu and sigma2 prior to ModelEst
  ctr->sigma2 = 1.0;
  tdlmModelEst(ctr);
  rHalfCauchyFC(&(ctr->nu), ctr->nTrees, 0.0);
  (ctr->tau).resize(ctr->nTrees);                 (ctr->tau).setOnes();
  for (t = 0; t < ctr->nTrees; ++t) {
    if (ctr->shrinkage > 1)
      rHalfCauchyFC(&(ctr->tau(t)), 0.0, 0.0);
  }
  ctr->nTerm.resize(ctr->nTrees);                 (ctr->nTerm).setOnes();
  ctr->nTerm2.resize(ctr->nTrees);                (ctr->nTerm2).setOnes();
  (ctr->Rmat).resize(ctr->n, ctr->nTrees);        (ctr->Rmat).setZero();

  // * Create progress meter
  progressMeter* prog = new progressMeter(ctr);

  // * Beginning of MCMC
  double sigmanu;
  std::size_t s;
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    Rcpp::checkUserInterrupt();
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }

    // * Update trees
    ctr->R += (ctr->Rmat).col(0); // Remove first tree est from R
    ctr->fhat.setZero();
    ctr->totTerm = 0;
    ctr->sumTermT2 = 0;
    ctr->totTermExp.setZero();
    ctr->sumTermT2Exp.setZero();
    ctr->expCount.setZero();
    ctr->mixCount.setZero();
    ctr->expInf.setZero();
    ctr->mixInf.setZero();
    if (ctr->interaction > 0) {
      (ctr->totTermMix).setZero();                (ctr->sumTermT2Mix).setZero();
    }

    for (t = 0; t < ctr->nTrees; ++t) {
      // Rcout << "\n" << t;
      tdlmmTreeMCMC(t, trees1[t], trees2[t], ctr, dgn, Exp);
      ctr->fhat += (ctr->Rmat).col(t);
      if (t < ctr->nTrees - 1)
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
    }
    // Rcout << " end trees";


    // * Pre-calculations for control and variance
    ctr->R = ctr->Y - ctr->fhat;
    ctr->sumTermT2 = (ctr->sumTermT2Exp).sum();
    ctr->totTerm = (ctr->totTermExp).sum();
    if(ctr->interaction) {
      ctr->sumTermT2 += (ctr->sumTermT2Mix).sum();
      ctr->totTerm += (ctr->totTermMix).sum();
    }

    // * Update model
    tdlmModelEst(ctr);

    // * Update variance parameters
    rHalfCauchyFC(&(ctr->nu), ctr->totTerm, ctr->sumTermT2 / ctr->sigma2);
    // if ((ctr->nu) != (ctr->nu)) 
    //   stop("\nNaN values (nu) occured during model run, rerun model.\n");
    sigmanu = ctr->sigma2 * ctr->nu;
    if ((ctr->shrinkage == 3) || (ctr->shrinkage == 1)) {
      for (int i = 0; i < ctr->nExp; ++i) {
        rHalfCauchyFC(&(ctr->muExp(i)), ctr->totTermExp(i),
                      ctr->sumTermT2Exp(i) / sigmanu);
        // if ((ctr->muExp)(i) != (ctr->muExp)(i)) 
        //   stop("\nNaN values (muExp) occured during model run, rerun model.\n");
        if (ctr->interaction) {
          for (int j = i; j < ctr->nExp; ++j) {
            if ((j > i) || (ctr->interaction == 2))
              rHalfCauchyFC(&(ctr->muMix(j, i)), ctr->totTermMix(j, i),
                            ctr->sumTermT2Mix(j, i) / sigmanu);
            // if ((ctr->muMix)(j,i) != (ctr->muMix)(j,i)) 
            //   stop("\nNaN values (muMix) occured during model run, rerun model.\n");
          } // end for loop updating interaction variances
        } // end if interactions
      } // end for loop updating exposure variances
    } // end if shrinkage == 3 or 1


    // * Update exposure selection probability
    if ((ctr->b > 1000) || (ctr->b > (0.5 * ctr->burn)))
      ctr->expProb = 
        rDirichlet(((ctr->expCount).array() + ctr->modKappa).matrix());
      
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

    // * Progress
    prog->printMark();
  } // End MCMC




  // * Setup data for return
  Eigen::MatrixXd DLM((dgn->DLMexp).size(), 7);
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
                            Named("fhat") = wrap(fhat),
                            Named("sigma2") = wrap(sigma2),
                            Named("nu") = wrap(nu),
                            Named("tau") = wrap(tau),
                            Named("termNodes") = wrap(termNodes),
                            Named("termNodes2") = wrap(termNodes2),
                            Named("tree1Exp") = wrap(tree1Exp),
                            Named("tree2Exp") = wrap(tree2Exp),
                            Named("expProb") = wrap(expProb),
                            Named("expInf") = wrap(expInf),
                            Named("expCount") = wrap(expCount),
                            Named("mixInf") = wrap(mixInf),
                            Named("mixCount") = wrap(mixCount),
                            Named("muExp") = wrap(muExp),
                            Named("muMix") = wrap(muMix),
                            Named("kappa") = wrap(kappa),
                            Named("treeAccept") = wrap(Accept)));
} // end tdlmm_Cpp






