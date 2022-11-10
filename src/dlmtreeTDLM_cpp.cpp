#include "RcppEigen.h"
#include "Node.h"
#include "NodeStruct.h"
#include "modDat.h"
#include "exposureDat.h"
#include "Fncs.h"
#include "modelCtr.h"
using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Lower;

double calcLogRatioTDLM(treeMHR mhr0, treeMHR mhr, double RtR, double RtZVgZtR,
                        dlmtreeCtr* ctr, double stepMhr, double treevar)
{
  // Rcout << stepMhr << " " << mhr.logVThetaChol << " " << mhr0.logVThetaChol <<
  // " " << mhr.beta << " " << mhr0.beta << " " << treevar << " " <<
  // mhr.totTerm << " " << mhr0.totTerm << "\n";
  if (ctr->binomial) {
    return(stepMhr + mhr.logVThetaChol - mhr0.logVThetaChol +
           0.5 * (mhr.beta - mhr0.beta) -
           log(treevar) * 0.5 * round(mhr.totTerm - mhr0.totTerm));
  } else {
    return(stepMhr + mhr.logVThetaChol - mhr0.logVThetaChol -
          (0.5 * (ctr->n + 1.0) *
            (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
            log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
          (log(treevar) * 0.5 * round(mhr.totTerm - mhr0.totTerm)));
  }
}

treeMHR dlmtreeNestedMHR(std::vector<Node*> modTerm, 
                         dlmtreeCtr* ctr, 
                         VectorXd ZtR, 
                         double treevar, 
                         bool updateNested)
{
  std::size_t s, s2;
  treeMHR out;
  int pXMod = int(modTerm.size());
  int totTerm = 0;
  std::vector<std::vector<Node*> > nestedTerm;
  for (s = 0; s < modTerm.size(); ++s) {
    nestedTerm.push_back(
      modTerm[s]->nodevals->nestedTree->listTerminal(updateNested));
    totTerm += nestedTerm[s].size();
  }
  MatrixXd ZtX, VgZtX, LInv;
  ZtX.resize(ctr->pZ, totTerm);
  VgZtX.resize(ctr->pZ, totTerm);
  
  if (modTerm.size() == 1) {

    // * Create design Xd, Z^tX, and VgZ^tX matrices
    MatrixXd Xtemp(ctr->n, totTerm);      Xtemp.setZero();
    for (s = 0; s < nestedTerm[0].size(); ++s) {
      // Rcout << s;
      Xtemp.col(s) = (nestedTerm[0][s]->nodevals)->X;
      
      if (ctr->binomial) {
        ZtX.col(s) = ctr->Zw.transpose() * (nestedTerm[0][s]->nodevals)->X;
        VgZtX.col(s) = ctr->Vg * ZtX.col(s);
        
      } else {
        ZtX.col(s) = (nestedTerm[0][s]->nodevals)->ZtX;
        VgZtX.col(s) = (nestedTerm[0][s]->nodevals)->VgZtX;
      }
    }

    // * Calculate covariance matrix V_theta      
    VectorXd XtVzInvR(totTerm);
    MatrixXd tempV(totTerm, totTerm);
    if (ctr->binomial) {
      tempV = Xtemp.transpose() * 
        (Xtemp.array().colwise() * ctr->Omega.array()).matrix();
      tempV.noalias() -= ZtX.transpose() * VgZtX;
      XtVzInvR = (Xtemp.array().colwise() * 
                  ctr->Omega.array()).matrix().transpose() * ctr->R;
      
    } else {
      tempV = Xtemp.transpose() * Xtemp;
      tempV.noalias() -= ZtX.transpose() * VgZtX;
      XtVzInvR = Xtemp.transpose() * ctr->R;
    }
    
    XtVzInvR.noalias() -= VgZtX.transpose() * ZtR;
    tempV.diagonal().array() += 1.0 / treevar;
    const MatrixXd VTheta = tempV.inverse();
    const MatrixXd VThetaChol = VTheta.llt().matrixL();
    const VectorXd ThetaHat = VTheta * XtVzInvR;
    
    out.draw = ThetaHat;
    out.draw.noalias() += 
      VThetaChol * as<VectorXd>(rnorm(totTerm, 0, sqrt(ctr->sigma2)));
    out.beta = ThetaHat.dot(XtVzInvR);
    out.logVThetaChol = VThetaChol.diagonal().array().log().sum();

    out.termT2 = (out.draw).dot(out.draw);
    out.totTerm = double(totTerm);
    out.nModTerm = 1.0;
    return(out);
    
  } // end if no modification

  std::vector<MatrixXd> X;
  MatrixXd Xtemp, Ztemp;
  VectorXd Rtemp, Otemp;
  MatrixXd XXiblock(totTerm, totTerm); XXiblock.setZero();
  VectorXd XtR(totTerm); XtR.setZero();
  int j, pX;
  int start = 0;
  
  // Rcout << "22";
  for (s = 0; s < modTerm.size(); ++s) {
    // Rcout << s;
    pX = nestedTerm[s].size();
    MatrixXd temp(ctr->n, pX); temp.setZero();
    X.push_back(temp);
    
    // create nested tree data matrices
    for (s2 = 0; s2 < nestedTerm[s].size(); ++s2) {
      X[s].col(s2) = nestedTerm[s][s2]->nodevals->X;
    } // end loop over nested tree nodes
    
    if ((modTerm[s]->nodevals->updateXmat) || ctr->binomial) {
      // Rcout << "u";
      Xtemp.resize(modTerm[s]->nodevals->idx.size(), pX); Xtemp.setZero();
      Ztemp.resize(modTerm[s]->nodevals->idx.size(), ctr->pZ); Ztemp.setZero();
      Rtemp.resize(modTerm[s]->nodevals->idx.size()); Rtemp.setZero();
      Otemp.resize(modTerm[s]->nodevals->idx.size()); Otemp.setZero();
      
      j = 0;
      for (int i : modTerm[s]->nodevals->idx) {
        Xtemp.row(j) = X[s].row(i);
        Ztemp.row(j) = ctr->Zw.row(i);
        Rtemp(j) = ctr->R(i);
        if (ctr->binomial)
          Otemp(j) = ctr->Omega(i);
        ++j;
      }
      
      if (ctr->binomial) {
        MatrixXd Xwtemp = (Xtemp.array().colwise() * Otemp.array()).matrix();
        modTerm[s]->nodevals->XtX = Xtemp.transpose() * Xwtemp;
        modTerm[s]->nodevals->ZtXmat = Ztemp.transpose() * Xtemp;
        modTerm[s]->nodevals->VgZtXmat = ctr->Vg * modTerm[s]->nodevals->ZtXmat;
        XtR.segment(start, pX) = Xwtemp.transpose() * Rtemp;
      } else {
        modTerm[s]->nodevals->XtX = Xtemp.transpose() * Xtemp;
        modTerm[s]->nodevals->ZtXmat = Ztemp.transpose() * Xtemp;
        modTerm[s]->nodevals->VgZtXmat = ctr->Vg * modTerm[s]->nodevals->ZtXmat;
        XtR.segment(start, pX) = Xtemp.transpose() * Rtemp;
        modTerm[s]->nodevals->updateXmat = 0;
      }
      
    } else { // reuse precalculated matrices
  
      for (int i : modTerm[s]->nodevals->idx)
        XtR.segment(start, pX).noalias() += X[s].row(i).transpose() * ctr->R(i);
      
    } // end update XtX and ZtX matrices    
    
    LInv.resize(pX, pX); LInv.setZero();
    LInv.diagonal().array() += 1.0 / treevar;
    if ((ctr->pZ < totTerm) || ctr->binomial) {
      XXiblock.block(start, start, pX, pX) =  
        (modTerm[s]->nodevals->XtX + LInv).inverse();
    } else {
      XXiblock.block(start, start, pX, pX) = modTerm[s]->nodevals->XtX + LInv;
    }
    ZtX.block(0, start, ctr->pZ, pX) = modTerm[s]->nodevals->ZtXmat;
    VgZtX.block(0, start, ctr->pZ, pX) = modTerm[s]->nodevals->VgZtXmat;
    
    start += pX;
  } // end loop over modifier nodes
  
  // Rcout << "e";
  MatrixXd VTheta(totTerm, totTerm);  VTheta.setZero();
  if ((ctr->pZ < totTerm) || ctr->binomial) {
    MatrixXd ZtXXi = ZtX * XXiblock;
    VTheta.triangularView<Lower>() = ZtXXi.transpose() *
      (ctr->VgInv - ZtXXi * ZtX.transpose()).inverse() * ZtXXi;
    VTheta.triangularView<Lower>() += XXiblock;
  } else {
    VTheta.triangularView<Lower>() = 
      (XXiblock - ZtX.transpose() * ctr->Vg * ZtX).inverse();
  }
  VectorXd XtVzInvR =   XtR - VgZtX.transpose() * ZtR;
  VectorXd ThetaHat =   VTheta.selfadjointView<Lower>() * XtVzInvR;
  MatrixXd VThetaChol = VTheta.selfadjointView<Lower>().llt().matrixL();
  // Rcout << ".";

  // Calculate fitted values
  out.draw = ThetaHat;
  out.draw.noalias() +=
    VThetaChol * as<VectorXd>(rnorm(totTerm, 0.0, sqrt(ctr->sigma2)));  
  // Rcout << ".";
  out.beta = ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
  out.termT2 = (out.draw).dot(out.draw);
  out.totTerm = double(totTerm); 
  out.nModTerm = double(pXMod);
  // Rcout << ".\n";
  return(out);
}




void dlmtreeTDLMTreeMCMC(int t, Node* modTree, NodeStruct* expNS,
                         dlmtreeCtr* ctr, dlmtreeLog *dgn,
                         modDat* Mod, exposureDat* Exp)
{
  // Rcout << ctr->sigma2 << " " << ctr->nu << " " << ctr->tau(t) << "\n";
  int step;
  int success =     0;
  double stepMhr =  0;
  double ratio =    0;
  double treevar =  ctr->nu * ctr->tau(t);
  VectorXd ZtR =    ctr->Zw.transpose() * ctr->R;
  double RtR =      0.0;
  double RtZVgZtR = 0.0;
  if (!(ctr->binomial)) {
    RtR =           ctr->R.dot(ctr->R);
    RtZVgZtR =      ZtR.dot(ctr->Vg * ZtR);
  }
  std::size_t s;
  std::vector<Node*> modTerm, dlmTerm, newDlmTerm, newModTerm;
  treeMHR mhr0, mhr;

  // -- List terminal nodes --
  modTerm = modTree->listTerminal();

  // -- Propose new modifier tree --
  switch (modTerm.size()) {
    case 1:  step = 0; break;
    case 2:  step = sampleInt(ctr->stepProbMod, 1 - ctr->stepProbMod[3]); break;
    default: step = sampleInt(ctr->stepProbMod, 1);
  }
  stepMhr = modProposeTree(modTree, Mod, ctr, step);
  success = modTree->isProposed();
  // Rcout << "1!";
  mhr0 = dlmtreeNestedMHR(modTerm, ctr, ZtR, treevar, 0);

  if (success) {
    newModTerm = modTree->listTerminal(1);
    
    // draw nested trees if grow or prune
    if (step < 2) {
      for (Node* tn : newModTerm) {
        if (tn->nodevals->nestedTree == 0) {
          tn->nodevals->updateXmat = 1;
          tn->nodevals->nestedTree = new Node(0, 1);
          tn->nodevals->nestedTree->nodestruct = expNS->clone();
          drawTree(tn->nodevals->nestedTree, tn->nodevals->nestedTree,
                   ctr->treePrior[0], ctr->treePrior[1]);
          for (Node* tn2 : tn->nodevals->nestedTree->listTerminal())
            Exp->updateNodeVals(tn2);
        }
      } // end loop over nested trees      
    } // end draw nested trees if grow or prune
    // Rcout << "2!";
    mhr =   dlmtreeNestedMHR(newModTerm, ctr, ZtR, treevar, 1);
    ratio = calcLogRatioTDLM(mhr0, mhr, RtR, RtZVgZtR, ctr, stepMhr, treevar);
    
    if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
      mhr0 = mhr;
      success = 2;
      modTree->accept();
      modTerm = modTree->listTerminal();
    }
  } // end modTree proposal
  if (success < 2)
    modTree->reject();
  
  // * Propose new nested tree at each modifier node
  for (Node* tn : modTerm) {
    dlmTerm = tn->nodevals->nestedTree->listTerminal();
    switch (dlmTerm.size()) {
      case 1: step = 0; break;
      default: step = sampleInt(ctr->stepProb, 1);
    }  
    stepMhr = tdlmProposeTree(tn->nodevals->nestedTree, Exp, ctr, step);
    success = tn->nodevals->nestedTree->isProposed();
    
    if (success) {
      MatrixXd XtX =       tn->nodevals->XtX;
      MatrixXd ZtXmat =    tn->nodevals->ZtXmat;
      MatrixXd VgZtXmat =  tn->nodevals->VgZtXmat;
      tn->nodevals->updateXmat =  1;
      // Rcout << "3!";
      mhr =   dlmtreeNestedMHR(modTerm, ctr, ZtR, treevar, 1);
      ratio = calcLogRatioTDLM(mhr0, mhr, RtR, RtZVgZtR, ctr, stepMhr, treevar);
      
      if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
        mhr0 = mhr;
        success = 2;
        tn->nodevals->nestedTree->accept();
        
      } else {
        tn->nodevals->XtX =       XtX;
        tn->nodevals->ZtXmat =    ZtXmat;
        tn->nodevals->VgZtXmat =  VgZtXmat;
      } // end MHR accept/reject
    } // end nested tree proposal
    
    if (success < 2)
      tn->nodevals->nestedTree->reject();
      
  } // end loop to update nested trees

  // -- Update variance and residuals --
  if (ctr->shrinkage) {
    rHalfCauchyFC(&(ctr->tau(t)), mhr0.totTerm,
                  mhr0.termT2 / (ctr->sigma2 * ctr->nu));
  }
  // ctr->Rmat.col(t) = mhr0.fitted;
  ctr->sumTermT2 +=   mhr0.termT2 / ctr->tau(t);
  ctr->totTerm +=     mhr0.totTerm;
  ctr->nTermMod(t) =  mhr0.nModTerm;
  ctr->nTerm(t) =     mhr0.totTerm;  
  
  // -- Count modifiers used in tree --
  VectorXd modCount = countMods(modTree, Mod);
  ctr->modCount +=    modCount;

  // -- Record --
  if (ctr->record > 0) {    
    for (int i = 0; i < modCount.size(); ++i) {
      if (modCount(i) > 0)
        ctr->modInf(i) += ctr->tau(t);
    }
  }

  // -- Update DLM partial estimate --
  std::string rule;
  std::vector<Node*> nested;
  VectorXd rec(9);
  MatrixXd Xmat;
  VectorXd draw;
  int drawIdx = 0;
  
  for (s = 0; s < modTerm.size(); ++s) {
    nested = modTerm[s]->nodevals->nestedTree->listTerminal();
    rule = modRuleStr(modTerm[s], Mod);
    Xmat.resize(ctr->n, nested.size());
    draw = mhr0.draw.segment(drawIdx, nested.size());
    drawIdx += nested.size();    
    
    for (std::size_t s2 = 0; s2 < nested.size(); ++s2) {
      Xmat.col(s2) = nested[s2]->nodevals->X;
      if (ctr->record > 0) {
        rec << ctr->record, t, s, s2, 
          (nested[s2]->nodestruct)->get(1),
          (nested[s2]->nodestruct)->get(2), 
          (nested[s2]->nodestruct)->get(3),
          (nested[s2]->nodestruct)->get(4), 
          draw(s2);
        dgn->termRule.push_back(rule);
        dgn->DLMexp.push_back(rec);
      } // end if record
    } // end loop over nested tree
    
    // calculate fitted values
    for (int i : modTerm[s]->nodevals->idx)
      ctr->Rmat(i, t) = Xmat.row(i) * draw;
  } // end loop over modifier nodes to update partial DLM est
} // end dlmtreeTDLMGaussian_TreeMCMC function


// [[Rcpp::export]]
Rcpp::List dlmtreeTDLM_cpp(const Rcpp::List model)
{
  int t;
  // ---- Set up general control variables ----
  dlmtreeCtr *ctr = new dlmtreeCtr;
  ctr->iter =         as<int>(model["nIter"]);
  ctr->burn =         as<int>(model["nBurn"]);
  ctr->thin =         as<int>(model["nThin"]);
  ctr->nRec =         floor(ctr->iter / ctr->thin);
  ctr->nTrees =       as<int>(model["nTrees"]);
  ctr->verbose =      as<bool>(model["verbose"]);
  ctr->diagnostics =  as<bool>(model["diagnostics"]);
  ctr->binomial =     as<bool>(model["binomial"]);
  ctr->stepProb =     as<std::vector<double> >(model["stepProbTDLM"]);
  ctr->stepProbMod =  as<std::vector<double> >(model["stepProbMod"]);
  ctr->treePrior =    as<std::vector<double> >(model["treePriorMod"]);
  ctr->treePriorMod = as<std::vector<double> >(model["treePriorTDLM"]);
  ctr->modZeta =      as<double>(model["zeta"]);
  ctr->modKappa =     100;
  ctr->shrinkage =    as<int>(model["shrinkage"]);

  // * Set up model data
  ctr->Y =      as<VectorXd>(model["Y"]);
  ctr->n =      ctr->Y.size();
  ctr->Z =      as<MatrixXd>(model["Z"]);
  ctr->Zw =     ctr->Z;
  ctr->pZ =     ctr->Z.cols();
  ctr->VgInv =  ctr->Z.transpose() * (ctr->Z);
  ctr->VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg =     ctr->VgInv.inverse();
  ctr->VgChol = ctr->Vg.llt().matrixL();
  
  // * Set up data for logistic model
  ctr->binomialSize.resize(ctr->n);           ctr->binomialSize.setZero();
  ctr->kappa.resize(ctr->n);                  ctr->kappa.setOnes();
  ctr->Omega.resize(ctr->n);                  ctr->Omega.setOnes();
  if (ctr->binomial) {
    ctr->binomialSize = as<VectorXd>(model["binomialSize"]);
    ctr->kappa = ctr->Y - 0.5 * (ctr->binomialSize);
    ctr->Y = ctr->kappa;
  }

  // * Setup modifier data
  modDat *Mod = new modDat(as<std::vector<int> >(model["modIsNum"]),
                           as<Rcpp::List>(model["modSplitIdx"]),
                           as<std::vector<int> >(model["fullIdx"]));
  NodeStruct *modNS;
  modNS = new ModStruct(Mod);
  ctr->pM = Mod->nMods;

  // * Setup exposure data
  ctr->XcenterIdx = 0;
  exposureDat *Exp;
  if (as<int>(model["nSplits"]) == 0) { // DLM
    if (ctr->binomial)
      Exp = new exposureDat(as<MatrixXd>(model["Tcalc"]));
    else
      Exp = new exposureDat(as<MatrixXd>(model["Tcalc"]),
                            ctr->Z, ctr->Vg);
  } else { // DLNM
    if (ctr->binomial)
      Exp = new exposureDat(as<MatrixXd>(model["X"]),
                            as<MatrixXd>(model["SE"]),
                            as<VectorXd>(model["Xsplits"]),
                            as<MatrixXd>(model["Xcalc"]),
                            as<MatrixXd>(model["Tcalc"]));
    else
      Exp = new exposureDat(as<MatrixXd>(model["X"]),
                            as<MatrixXd>(model["SE"]),
                            as<VectorXd>(model["Xsplits"]),
                            as<MatrixXd>(model["Xcalc"]),
                            as<MatrixXd>(model["Tcalc"]),
                            ctr->Z, ctr->Vg);
  }
  ctr->pX = Exp->pX;
  ctr->nSplits = Exp->nSplits;

  // * Pre-calculate single node tree matrices
  ctr->X1 = (Exp->Tcalc).col(ctr->pX - 1);
  ctr->ZtX1 = (ctr->Z).transpose() * (ctr->X1);
  ctr->VgZtX1 = (ctr->Vg).selfadjointView<Lower>() * (ctr->ZtX1);
  ctr->VTheta1Inv = (ctr->X1).dot(ctr->X1) - (ctr->ZtX1).dot(ctr->VgZtX1);

  // * Create trees
  NodeStruct *expNS;
  expNS = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                          as<VectorXd>(model["splitProb"]),
                          as<VectorXd>(model["timeProb"]));
  std::vector<Node*> modTrees;
  for (t = 0; t < ctr->nTrees; ++t) {
    // create modifier tree
    modTrees.push_back(new Node(0, 1));
    modTrees[t]->nodestruct = modNS->clone();
    Mod->updateNodeVals(modTrees[t]);
    
    // create nested tree
    modTrees[t]->nodevals->nestedTree = new Node(0, 1);
    modTrees[t]->nodevals->nestedTree->nodestruct = expNS->clone();
    drawTree(modTrees[t]->nodevals->nestedTree, 
             modTrees[t]->nodevals->nestedTree,
             ctr->treePrior[0], ctr->treePrior[1]);
    for (Node* tn : modTrees[t]->nodevals->nestedTree->listTerminal())
      Exp->updateNodeVals(tn);
  }
  delete modNS;
  ctr->nTerm.resize(ctr->nTrees);         ctr->nTerm.setOnes(); 
  ctr->nTermMod.resize(ctr->nTrees);      ctr->nTermMod.setOnes();
  ctr->Rmat.resize(ctr->n, ctr->nTrees);  ctr->Rmat.setZero();
  ctr->modCount.resize(ctr->pM);          ctr->modCount.setZero();
  ctr->modInf.resize(ctr->pM);            ctr->modInf.setZero();


  // * Setup model logs
  dlmtreeLog *dgn = new dlmtreeLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec);    (dgn->gamma).setZero();
  (dgn->sigma2).resize(ctr->nRec);            (dgn->sigma2).setZero();
  (dgn->nu).resize(ctr->nRec);                (dgn->nu).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec);  (dgn->tau).setZero();
  (dgn->fhat).resize(ctr->n);                 (dgn->fhat).setZero();
  (dgn->modProb).resize(ctr->pM, ctr->nRec);  (dgn->modProb).setZero();
  (dgn->modCount).resize(ctr->pM, ctr->nRec); (dgn->modCount).setZero();
  (dgn->modInf).resize(ctr->pM, ctr->nRec);   (dgn->modInf).setZero();
  (dgn->modKappa).resize(ctr->nRec);          (dgn->modKappa).setZero();
  (dgn->termNodesMod).resize(ctr->nTrees, ctr->nRec);
    (dgn->termNodesMod).setZero();
  (dgn->termNodesDLM).resize(ctr->nTrees, ctr->nRec);
    (dgn->termNodesDLM).setZero();
    

  // * Initial draws
  (ctr->fhat).resize(ctr->n);                 (ctr->fhat).setZero();
  ctr->R = ctr->Y;
  (ctr->gamma).resize(ctr->pZ);
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  ctr->nu = 1.0; // ! Need to define for first update of sigma2
  ctr->sigma2 = 1.0;
  tdlmModelEst(ctr);
  rHalfCauchyFC(&(ctr->nu), ctr->nTrees, 0.0);
  (ctr->tau).resize(ctr->nTrees); (ctr->tau).setOnes();
  if (ctr->shrinkage) {
    for (t = 0; t < ctr->nTrees; t++)
      rHalfCauchyFC(&(ctr->tau(t)), 0.0, 0.0);
  }
  
  // Create progress meter
  progressMeter* prog = new progressMeter(ctr);

  std::size_t s;
  
  // * Begin MCMC
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    ctr->record = 0;
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0))
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);

    // * Update trees 
    ctr->R += ctr->Rmat.col(0);
    ctr->fhat.setZero();
    ctr->totTerm = 0.0; 
    ctr->sumTermT2 = 0.0;
    ctr->modCount.setZero();
    ctr->modInf.setZero();
    
    for (t = 0; t < ctr->nTrees; t++) {
      dlmtreeTDLMTreeMCMC(t, modTrees[t], expNS, ctr, dgn, Mod, Exp);
      ctr->fhat += (ctr->Rmat).col(t);
      if (t < ctr->nTrees - 1)
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
    } // end update trees

    // * Update model
    ctr->R = ctr->Y - ctr->fhat;
    tdlmModelEst(ctr);
    rHalfCauchyFC(&(ctr->nu), ctr->totTerm, ctr->sumTermT2 / ctr->sigma2);
    
    // * Update modifier selection
    if ((ctr->b > 1000) || (ctr->b > (0.5 * ctr->burn))) {
      double beta = R::rbeta(ctr->modZeta, 1.0);
      double modKappaNew = beta * ctr->pM / (1 - beta);
      double mhrDir =
        logDirichletDensity(Mod->modProb,
                            (ctr->modCount.array() +
                             modKappaNew / ctr->pM).matrix()) -
        logDirichletDensity(Mod->modProb,
                            (ctr->modCount.array() +
                             ctr->modKappa / ctr->pM).matrix());
      if (log(R::runif(0, 1)) < mhrDir)
        ctr->modKappa = modKappaNew;

      Mod->modProb = rDirichlet((ctr->modCount.array() +
                                 ctr->modKappa / ctr->pM).matrix());
    } // end modifier selection

    // -- Record --
    if (ctr->record > 0) {
      (dgn->gamma).col(ctr->record - 1) = ctr->gamma;
      (dgn->sigma2)(ctr->record - 1) = ctr->sigma2;
      (dgn->nu)(ctr->record - 1) = ctr->nu;
      (dgn->tau).col(ctr->record - 1) = ctr->tau;
      (dgn->termNodesDLM).col(ctr->record - 1) = ctr->nTerm;
      (dgn->termNodesMod).col(ctr->record - 1) = ctr->nTermMod;
      (dgn->modProb).col(ctr->record - 1) = Mod->modProb;
      (dgn->modCount).col(ctr->record - 1) = ctr->modCount;
      (dgn->modInf).col(ctr->record - 1) = ctr->modInf / ctr->modInf.maxCoeff();
      (dgn->modKappa)(ctr->record - 1) = ctr->modKappa;
      dgn->fhat += ctr->fhat;
    } // end record

    // * Mark progress
    prog->printMark();
  } // end MCMC


  // * Prepare outout
  MatrixXd TreeStructs;
  TreeStructs.resize((dgn->DLMexp).size(), 9);
  for (s = 0; s < (dgn->DLMexp).size(); ++s)
    TreeStructs.row(s) = dgn->DLMexp[s];
  Rcpp::StringVector termRule(dgn->termRule.size());
  termRule = dgn->termRule;
  MatrixXd termNodesDLM = (dgn->termNodesDLM).transpose();

  VectorXd sigma2 = dgn->sigma2;
  VectorXd nu = dgn->nu;
  MatrixXd tau = (dgn->tau).transpose();
  VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  MatrixXd gamma = (dgn->gamma).transpose();

  MatrixXd termNodesMod = (dgn->termNodesMod).transpose();
  VectorXd kappa = dgn->modKappa;
  MatrixXd modProb = (dgn->modProb).transpose();
  MatrixXd modCount = (dgn->modCount).transpose();
  MatrixXd modInf = (dgn->modInf).transpose();

  MatrixXd modAccept((dgn->treeModAccept).size(), 5);
  MatrixXd dlmAccept((dgn->treeDLMAccept).size(), 5);

  delete prog;
  delete ctr;
  delete dgn;
  delete Exp;
  delete Mod;
  delete expNS;
  for (s = 0; s < modTrees.size(); ++s)
    delete modTrees[s];

  return(Rcpp::List::create(Named("TreeStructs") = wrap(TreeStructs),
                            Named("termRules") = wrap(termRule),
                            Named("termNodesDLM") = wrap(termNodesDLM),
                            Named("fhat") = wrap(fhat),
                            Named("sigma2") = wrap(sigma2),
                            Named("nu") = wrap(nu),
                            Named("tau") = wrap(tau),
                            Named("gamma") = wrap(gamma),
                            Named("termNodesMod") = wrap(termNodesMod),
                            Named("kappa") = wrap(kappa),
                            Named("modProb") = wrap(modProb),
                            Named("modCount") = wrap(modCount),
                            Named("modInf") = wrap(modInf),
                            Named("treeModAccept") = wrap(modAccept),
                            Named("treeDLMAccept") = wrap(dlmAccept)));

} // end dlmtreeTDLMGaussian