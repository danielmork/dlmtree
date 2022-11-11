#include "RcppEigen.h"
#include "Node.h"
#include "NodeStruct.h"
#include "modDat.h"
#include "exposureDat.h"
#include "Fncs.h"
#include "modelCtr.h"
using namespace Rcpp;



double calcLogRatio(treeMHR mhr0, treeMHR mhr, double RtR, double RtZVgZtR,
                    dlmtreeCtr* ctr, double stepMhr, double treevar)
{
  return(stepMhr +
         mhr.logVThetaChol - mhr0.logVThetaChol -
         (0.5 * (ctr->n + 1.0) *
           (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
           log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
         (log(treevar) * 0.5 * round(mhr.totTerm - mhr0.totTerm)));
}

treeMHR dlmtreeTDLMNested_MHR(std::vector<Node*> modTerm,
                   dlmtreeCtr* ctr, Eigen::VectorXd ZtR,
                   double treevar, bool updateNested)
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

  std::vector<Eigen::MatrixXd> X;
  Eigen::MatrixXd ZtX, VgZtX, LInv;
  ZtX.resize(ctr->pZ, totTerm); ZtX.setZero();
  VgZtX.resize(ctr->pZ, totTerm); VgZtX.setZero();
  Eigen::MatrixXd Xtemp, Ztemp;
  Eigen::VectorXd Rtemp;
  Eigen::MatrixXd XXiblock(totTerm, totTerm); XXiblock.setZero();
  Eigen::VectorXd XtR(totTerm); XtR.setZero();
  int j, pX;
  int start = 0;
  
  for (s = 0; s < modTerm.size(); ++s) {
    pX = nestedTerm[s].size();
    Eigen::MatrixXd temp(ctr->n, pX); temp.setZero();
    X.push_back(temp);
    
    // create nested tree data matrices
    for (s2 = 0; s2 < nestedTerm[s].size(); ++s2) {
      X[s].col(s2) = nestedTerm[s][s2]->nodevals->X;
    } // end loop over nested tree nodes
    
    if (modTerm[s]->nodevals->updateXmat) {                                  
      Xtemp.resize(modTerm[s]->nodevals->idx.size(), pX); 
      Xtemp.setZero();
      Ztemp.resize(modTerm[s]->nodevals->idx.size(), ctr->pZ); 
      Ztemp.setZero();
      Rtemp.resize(modTerm[s]->nodevals->idx.size()); 
      Rtemp.setZero();
      
      j = 0;
      for (int i : modTerm[s]->nodevals->idx) {
        Xtemp.row(j) = X[s].row(i);
        Ztemp.row(j) = ctr->Z.row(i);
        Rtemp(j) = ctr->R(i);
        ++j;
      }
      
      modTerm[s]->nodevals->XtX.resize(pX, pX);
      modTerm[s]->nodevals->XtX = Xtemp.transpose() * Xtemp;
      modTerm[s]->nodevals->ZtXmat.resize(ctr->pZ, pX);
      modTerm[s]->nodevals->ZtXmat = Ztemp.transpose() * Xtemp;
      modTerm[s]->nodevals->VgZtXmat.resize(ctr->pZ, pX);
      modTerm[s]->nodevals->VgZtXmat = ctr->Vg * modTerm[s]->nodevals->ZtXmat;
      modTerm[s]->nodevals->updateXmat = 0;
      XtR.segment(start, pX) = Xtemp.transpose() * Rtemp;
      
    } else { // reuse precalculated matrices
      
      for (int i : modTerm[s]->nodevals->idx) {
        XtR.segment(start, pX).noalias() += X[s].row(i).transpose() * ctr->R(i);
      }
      
    } // end update XtX and ZtX matrices    
    
    LInv.resize(pX, pX); LInv.setZero();
    LInv.diagonal().array() += 1.0 / treevar;
    XXiblock.block(start, start, pX, pX) =  
      (modTerm[s]->nodevals->XtX + LInv).inverse();
    ZtX.block(0, start, ctr->pZ, pX) = modTerm[s]->nodevals->ZtXmat;
    VgZtX.block(0, start, ctr->pZ, pX) = modTerm[s]->nodevals->VgZtXmat;
    
    start += pX;
  } // end loop over modifier nodes
  
  
  Eigen::MatrixXd ZtXXi = ZtX * XXiblock;
  Eigen::MatrixXd VTheta = XXiblock;
  VTheta.noalias() += ZtXXi.transpose() *
    (ctr->VgInv - ZtXXi * ZtX.transpose()).inverse() * ZtXXi;
  Eigen::VectorXd XtVzInvR = XtR;
  XtVzInvR.noalias() -= VgZtX.transpose() * ZtR;
  Eigen::VectorXd ThetaHat = VTheta * XtVzInvR;
  Eigen::MatrixXd VThetaChol = VTheta.llt().matrixL();

  // Calculate fitted values
  out.draw = ThetaHat;
  out.draw.noalias() +=
    VThetaChol * as<Eigen::VectorXd>(rnorm(totTerm, 0.0, sqrt(ctr->sigma2)));
  out.fitted.resize(ctr->n);
  Eigen::VectorXd drawTemp;
  start = 0;
  for (s = 0; s < modTerm.size(); ++s) {    
    pX = nestedTerm[s].size();
    drawTemp.resize(pX);    
    drawTemp = out.draw.segment(start, pX);
    
    for (int i : modTerm[s]->nodevals->idx)
      out.fitted(i) = X[s].row(i) * drawTemp;
      
    start += pX;
  } // end loop over modifier nodes to calculate dr
  
  out.beta = ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
  out.termT2 = (out.draw).dot(out.draw);
  out.totTerm = double(totTerm); 
  out.nModTerm = double(pXMod);
  return(out);
}

void dlmtreeTDLMNestedGaussian_TreeMCMC(int t, Node* modTree, NodeStruct* expNS,
                                   dlmtreeCtr* ctr, dlmtreeLog *dgn,
                                   modDat* Mod, exposureDat* Exp)
{
  int step;
  int success = 0;
  double stepMhr = 0;
  double ratio = 0;
  double treevar = (ctr->nu) * (ctr->tau)(t);
  std::size_t s;
  std::vector<Node*> modTerm, dlmTerm, newDlmTerm, newModTerm;
  Eigen::VectorXd ZtR = (ctr->Z).transpose() * (ctr->R);
  double RtR = (ctr->R).dot(ctr->R);
  double RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);
  treeMHR mhr0, mhr;

  // -- List terminal nodes --
  modTerm = modTree->listTerminal();

  // -- Propose new modifier tree --
  switch (modTerm.size()) {
    case 1: step = 0; break;
    case 2: step = sampleInt(ctr->stepProbMod, 1 - ctr->stepProbMod[3]); break;
    default: step = sampleInt(ctr->stepProbMod, 1);
  }
  stepMhr = modProposeTree(modTree, Mod, ctr, step);
  success = modTree->isProposed();
  mhr0 = dlmtreeTDLMNested_MHR(modTerm, ctr, ZtR, treevar, 0);
  
  if (success && (stepMhr == stepMhr)) {
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
      }
      
    } // end draw nested trees if grow or prune
    
    mhr = dlmtreeTDLMNested_MHR(newModTerm, ctr, ZtR, treevar, 0);
    ratio = calcLogRatio(mhr0, mhr, RtR, RtZVgZtR, ctr, stepMhr, treevar);
    
    if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
      mhr0 = mhr;
      success = 2;
      modTree->accept();
      modTerm = modTree->listTerminal();
    }
  } // end modTree proposal
  if (success < 2)
    modTree->reject();
  
  // -- Propose new nested tree at each modifier node --
  for (Node* tn : modTerm) {
    dlmTerm = tn->nodevals->nestedTree->listTerminal();
    switch (dlmTerm.size()) {
      case 1: step = 0; break;
      default: step = sampleInt(ctr->stepProb, 1);
    }  
    stepMhr = tdlmProposeTree(tn->nodevals->nestedTree, Exp, ctr, step);
    success = tn->nodevals->nestedTree->isProposed();
    
    if (success && (stepMhr == stepMhr)) {
      Eigen::MatrixXd XtX = tn->nodevals->XtX;
      Eigen::MatrixXd ZtXmat = tn->nodevals->ZtXmat;
      Eigen::MatrixXd VgZtXmat = tn->nodevals->VgZtXmat;
      tn->nodevals->updateXmat = 1;
      mhr = dlmtreeTDLMNested_MHR(modTerm, ctr, ZtR, treevar, 1);
      ratio = calcLogRatio(mhr0, mhr, RtR, RtZVgZtR, ctr, stepMhr, treevar);
      
      if ((log(R::runif(0, 1)) < ratio) && (ratio == ratio)) {
        mhr0 = mhr;
        success = 2;
        tn->nodevals->nestedTree->accept();
      } else {
        tn->nodevals->XtX = XtX;
        tn->nodevals->ZtXmat = ZtXmat;
        tn->nodevals->VgZtXmat = VgZtXmat;
      } // end MHR accept/reject
    } // end nested tree proposal
    if (success < 2)
      tn->nodevals->nestedTree->reject();
  } // end loop to update nested trees

  // -- Update variance and residuals --
  if (ctr->shrinkage)
    rHalfCauchyFC(&(ctr->tau(t)), mhr0.totTerm,
                  mhr0.termT2 / (ctr->sigma2 * ctr->nu));
    
  ctr->Rmat.col(t) = mhr0.fitted;
  ctr->sumTermT2 += mhr0.termT2 / (ctr->tau(t));
  ctr->totTerm += mhr0.totTerm;
  ctr->nTermMod(t) = mhr0.nModTerm;
  ctr->nTerm(t) = mhr0.totTerm;  
  
  // -- Count modifiers used in tree --
  Eigen::VectorXd modCount = countMods(modTree, Mod);
  ctr->modCount += modCount;

  // -- Record --
  if (ctr->record > 0) {    
    for (int i = 0; i < modCount.size(); ++i) {
      if (modCount(i) > 0)
        ctr->modInf(i) += ctr->tau(t);
    }

    // -- Update DLM partial estimate --
    std::string rule;
    Eigen::VectorXd rec(9);
    Eigen::VectorXd draw(ctr->pX);
    Node* tn;
    std::vector<Node*> nested;
    int drawIdx = 0;
    
    for (s = 0; s < modTerm.size(); ++s) {
      tn = modTerm[s];
      nested = tn->nodevals->nestedTree->listTerminal();
      rule = modRuleStr(tn, Mod);
      
      for (std::size_t s2 = 0; s2 < nested.size(); ++s2) {
        if (ctr->nSplits == 0) { // DLM
          for (int t2 = nested[s2]->nodestruct->get(3) - 1;
              t2 < nested[s2]->nodestruct->get(4); ++t2)
            draw(t2) = mhr0.draw(drawIdx);
        }
        
        rec << ctr->record, t, s, s2, 
          (nested[s2]->nodestruct)->get(1),
          (nested[s2]->nodestruct)->get(2), 
          (nested[s2]->nodestruct)->get(3),
          (nested[s2]->nodestruct)->get(4), 
          mhr0.draw(drawIdx);
        dgn->termRule.push_back(rule);
        dgn->DLMexp.push_back(rec);
        ++drawIdx;
      } // end loop over nested tree
    } // end loop over modifier nodes to update partial DLM est
  } // end record


} // end dlmtreeTDLMGaussian_TreeMCMC function

// [[Rcpp::export]]
Rcpp::List dlmtreeTDLMNestedGaussian(const Rcpp::List model)
{
  int t;
  // ---- Set up general control variables ----
  dlmtreeCtr *ctr = new dlmtreeCtr;
  ctr->iter = as<int>(model["nIter"]);
  ctr->burn = as<int>(model["nBurn"]);
  ctr->thin = as<int>(model["nThin"]);
  ctr->nRec = floor(ctr->iter / ctr->thin);
  ctr->nTrees = as<int>(model["nTrees"]);
  ctr->Y = as<Eigen::VectorXd>(model["Y"]);
  ctr->n = (ctr->Y).size();
  ctr->modZeta = as<double>(model["zeta"]);
  ctr->modKappa = 100;
  
  ctr->binomial = 0;
  ctr->verbose = bool (model["verbose"]);
  ctr->diagnostics = bool (model["diagnostics"]);
  ctr->stepProb = as<std::vector<double> >(model["stepProbTDLM"]);
  ctr->stepProbMod = as<std::vector<double> >(model["stepProbMod"]);
  ctr->treePrior = as<std::vector<double> >(model["treePriorMod"]);
  ctr->treePriorMod = as<std::vector<double> >(model["treePriorTDLM"]);
  ctr->shrinkage = as<int>(model["shrinkage"]);

  ctr->Z = as<Eigen::MatrixXd>(model["Z"]);
  ctr->Zw = ctr->Z;
  ctr->pZ = (ctr->Z).cols();
  ctr->VgInv = (ctr->Z).transpose() * (ctr->Z);
  ctr->VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg = ctr->VgInv.inverse();
  ctr->VgChol = (ctr->Vg).llt().matrixL();

  // ---- Setup modifier data ----
  modDat *Mod = new modDat(as<std::vector<int> >(model["modIsNum"]),
                           as<Rcpp::List>(model["modSplitIdx"]),
                           as<std::vector<int> >(model["fullIdx"]));
  NodeStruct *modNS;
  modNS = new ModStruct(Mod);
  ctr->pM = Mod->nMods;

  // * Setup exposure data
  ctr->XcenterIdx = 0;
  exposureDat *Exp;
  if (as<int>(model["nSplits"]) == 0) {
    Exp = new exposureDat(as<Eigen::MatrixXd>(model["Tcalc"]),
                          ctr->Z, ctr->Vg);
  } else {
    ctr->XcenterIdx = as<double>(model["XcenterIdx"]);
    Exp = new exposureDat(as<Eigen::MatrixXd>(model["X"]),
                          as<Eigen::MatrixXd>(model["SE"]),
                          as<Eigen::VectorXd>(model["Xsplits"]),
                          as<Eigen::MatrixXd>(model["Xcalc"]),
                          as<Eigen::MatrixXd>(model["Tcalc"]),
                          ctr->Z, ctr->Vg);
  }
  ctr->pX = Exp->pX;
  ctr->nSplits = Exp->nSplits;

  // * Pre-calculate single node tree matrices
  ctr->X1 = (Exp->Tcalc).col(ctr->pX - 1);
  ctr->ZtX1 = (ctr->Z).transpose() * (ctr->X1);
  ctr->VgZtX1 = (ctr->Vg).selfadjointView<Eigen::Lower>() * (ctr->ZtX1);
  ctr->VTheta1Inv = (ctr->X1).dot(ctr->X1) - (ctr->ZtX1).dot(ctr->VgZtX1);

  // * Create trees
  NodeStruct *expNS;
  expNS = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                          as<Eigen::VectorXd>(model["splitProb"]),
                          as<Eigen::VectorXd>(model["timeProb"]));
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


  // * Logs
  dlmtreeLog *dgn = new dlmtreeLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec); (dgn->gamma).setZero();
  (dgn->sigma2).resize(ctr->nRec); (dgn->sigma2).setZero();
  (dgn->nu).resize(ctr->nRec); (dgn->nu).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec); (dgn->tau).setZero();
  (dgn->fhat).resize(ctr->n); (dgn->fhat).setZero();
  (dgn->modProb).resize(ctr->pM, ctr->nRec); (dgn->modProb).setZero();
  (dgn->modCount).resize(ctr->pM, ctr->nRec); (dgn->modCount).setZero();
  (dgn->modInf).resize(ctr->pM, ctr->nRec); (dgn->modInf).setZero();
  (dgn->modKappa).resize(ctr->nRec); (dgn->modKappa).setZero();
  (dgn->termNodesMod).resize(ctr->nTrees, ctr->nRec);
    (dgn->termNodesMod).setZero();
  (dgn->termNodesDLM).resize(ctr->nTrees, ctr->nRec);
    (dgn->termNodesDLM).setZero();

  // * Initial draws
  ctr->fhat.resize(ctr->n);             ctr->fhat.setZero();
  ctr->tau.resize(ctr->nTrees);         ctr->tau.setOnes();
  ctr->nu = 1.0; // ! Need to define for first update of sigma2
  ctr->sigma2 = 1.0;
  ctr->R = ctr->Y;
  ctr->gamma.resize(ctr->pZ);
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  tdlmModelEst(ctr);
  rHalfCauchyFC(&(ctr->nu), ctr->nTrees, 0.0);
  if (ctr->shrinkage) {
    for (t = 0; t < ctr->nTrees; t++)
      rHalfCauchyFC(&(ctr->tau(t)), 0.0, 0.0);
  }

  // Create progress meter
  progressMeter* prog = new progressMeter(ctr);

  std::size_t s;
  // * Begin MCMC
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    Rcpp::checkUserInterrupt();
    ctr->record = 0;
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0))
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);

    // * Update trees
    ctr->R += (ctr->Rmat).col(0);
    (ctr->fhat).setZero();
    ctr->totTerm = 0.0; 
    ctr->sumTermT2 = 0.0;
    ctr->modCount.setZero();
    ctr->modInf.setZero();

    for (t = 0; t < ctr->nTrees; t++) {
      dlmtreeTDLMNestedGaussian_TreeMCMC(t, modTrees[t], expNS,
                                         ctr, dgn, Mod, Exp);
      ctr->fhat += (ctr->Rmat).col(t);
      if (t < ctr->nTrees - 1)
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
    }

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
      if (log(R::runif(0, 1)) < mhrDir) {
        ctr->modKappa = modKappaNew;
      }

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

    // -- Mark progress --
    prog->printMark();
  } // end MCMC


  // * Prepare outout
  Eigen::MatrixXd TreeStructs;
  TreeStructs.resize((dgn->DLMexp).size(), 9);
  for (s = 0; s < (dgn->DLMexp).size(); ++s)
    TreeStructs.row(s) = dgn->DLMexp[s];
  Rcpp::StringVector termRule(dgn->termRule.size());
  termRule = dgn->termRule;
  Eigen::MatrixXd termNodesDLM = (dgn->termNodesDLM).transpose();

  Eigen::VectorXd sigma2 = dgn->sigma2;
  Eigen::VectorXd nu = dgn->nu;
  Eigen::MatrixXd tau = (dgn->tau).transpose();
  Eigen::VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  Eigen::MatrixXd gamma = (dgn->gamma).transpose();

  Eigen::MatrixXd termNodesMod = (dgn->termNodesMod).transpose();
  Eigen::VectorXd kappa = dgn->modKappa;
  Eigen::MatrixXd modProb = (dgn->modProb).transpose();
  Eigen::MatrixXd modCount = (dgn->modCount).transpose();
  Eigen::MatrixXd modInf = (dgn->modInf).transpose();

  Eigen::MatrixXd modAccept((dgn->treeModAccept).size(), 5);
  Eigen::MatrixXd dlmAccept((dgn->treeDLMAccept).size(), 5);

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