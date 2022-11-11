#include "RcppEigen.h"
#include "Node.h"
#include "NodeStruct.h"
#include "modDat.h"
#include "exposureDat.h"
#include "Fncs.h"
#include "modelCtr.h"
using namespace Rcpp;


void dlmtreeTDLMFixed_Gaussian_TreeMCMC(
  int t, std::vector<Node*> fixedNodes,
  dlmtreeCtr* ctr, dlmtreeLog *dgn, exposureDat* Exp);

treeMHR dlmtreeTDLMFixedMHR(std::vector<Node*> fixedNodes,
                   dlmtreeCtr* ctr, 
                   Eigen::VectorXd ZtR,
                   double treevar,
                   double updateNested);

double calcLogRatioFixedTDLM(treeMHR mhr0, treeMHR mhr, double RtR, 
                             double RtZVgZtR, dlmtreeCtr* ctr, 
                             double stepMhr, double treevar);

// [[Rcpp::export]]
Rcpp::List dlmtreeTDLMFixedGaussian(const Rcpp::List model)
{
  int t;
  // ---- Set up general control variables ----
  dlmtreeCtr *ctr = new dlmtreeCtr;
  ctr->iter = as<int>(model["nIter"]);
  ctr->burn = as<int>(model["nBurn"]);
  ctr->thin = as<int>(model["nThin"]);
  ctr->nRec = floor(ctr->iter / ctr->thin);
  ctr->nTrees = as<int>(model["nTrees"]);
  ctr->verbose = bool (model["verbose"]);
  ctr->diagnostics = bool (model["diagnostics"]);
  ctr->binomial = 0;
  ctr->stepProb = as<std::vector<double> >(model["stepProbTDLM"]);
  ctr->treePrior = as<std::vector<double> >(model["treePriorMod"]);
  
  // ---- Load data ----
  ctr->Y = as<Eigen::VectorXd>(model["Y"]);
  ctr->n = (ctr->Y).size();
  ctr->Z = as<Eigen::MatrixXd>(model["Z"]);
  ctr->Zw = ctr->Z;
  ctr->pZ = (ctr->Z).cols();
  ctr->VgInv = (ctr->Z).transpose() * (ctr->Z);
  ctr->VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg = ctr->VgInv.inverse();
  ctr->VgChol = (ctr->Vg).llt().matrixL();

  // ---- Pre-calculate single node tree matrices ----
  ctr->XcenterIdx = 0;
  exposureDat *Exp;
  if (as<int>(model["nSplits"]) == 0) { // DLM
    Exp = new exposureDat(as<Eigen::MatrixXd>(model["Tcalc"]),
                          ctr->Z, ctr->Vg);
  } else { // DLNM
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
  ctr->X1 = (Exp->Tcalc).col(ctr->pX - 1);
  ctr->ZtX1 = (ctr->Z).transpose() * (ctr->X1);
  ctr->VgZtX1 = (ctr->Vg).selfadjointView<Eigen::Lower>() * (ctr->ZtX1);
  ctr->VTheta1Inv = (ctr->X1).dot(ctr->X1) - (ctr->ZtX1).dot(ctr->VgZtX1);

  // ---- Create trees ----
  NodeStruct *expNS;
  expNS = new DLNMStruct(0, ctr->nSplits + 1, 1, int(ctr->pX),
                          as<Eigen::VectorXd>(model["splitProb"]),
                          as<Eigen::VectorXd>(model["timeProb"]));

  std::vector<std::vector<Node*> > fixedTrees;
  std::vector<Node*> fixedNodes;
  Rcpp::List fixedIdx = as<Rcpp::List>(model["fixedIdx"]);
  Eigen::MatrixXd Xtemp, Ztemp;
  
  for (t = 0; t < ctr->nTrees; ++t) {
    fixedNodes.clear();
    for (int i = 0; i < fixedIdx.size(); ++i) {
      fixedNodes.push_back(new Node(0, 1));
      fixedNodes[i]->nodevals = new NodeVals(ctr->n);
      fixedNodes[i]->nodevals->idx = as<std::vector<int> >(fixedIdx[i]);
    
      // create nested tree
      fixedNodes[i]->nodevals->nestedTree = new Node(0, 1);
      fixedNodes[i]->nodevals->nestedTree->nodestruct = expNS->clone();
      drawTree(fixedNodes[i]->nodevals->nestedTree, 
               fixedNodes[i]->nodevals->nestedTree,
               ctr->treePrior[0], ctr->treePrior[1]);
      for (Node* tn : fixedNodes[i]->nodevals->nestedTree->listTerminal())
        Exp->updateNodeVals(tn);
    }
    fixedTrees.push_back(fixedNodes);
  }


  // ---- Logs ----
  dlmtreeLog *dgn = new dlmtreeLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec); (dgn->gamma).setZero();
  (dgn->sigma2).resize(ctr->nRec); (dgn->sigma2).setZero();
  (dgn->nu).resize(ctr->nRec); (dgn->nu).setZero();
  (dgn->phi).resize(ctr->nRec); (dgn->phi).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec); (dgn->tau).setZero();
  (dgn->fhat).resize(ctr->n); (dgn->fhat).setZero();

  // ---- Initial draws ----
  (ctr->fhat).resize(ctr->n); (ctr->fhat).setZero();
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
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  ctr->nu = 1.0; // Need to define for first update of sigma2
  ctr->sigma2 = 1.0;
  tdlmModelEst(ctr);
  double xiInv = R::rgamma(1, 0.5);
  ctr->nu = 1.0 / R::rgamma(0.5 * ctr->nTrees + 0.5, 1.0 / xiInv);
  (ctr->tau).resize(ctr->nTrees); (ctr->tau).setOnes();
  for (t = 0; t < ctr->nTrees; t++) {
    xiInv = R::rgamma(1, 0.5);
    (ctr->tau)(t) = 1.0 / R::rgamma(0.5, 1.0 / xiInv);
  }
  (ctr->Rmat).resize(ctr->n, ctr->nTrees); (ctr->Rmat).setZero();
  
  // Create progress meter
  progressMeter* prog = new progressMeter(ctr);

  std::size_t s;
  // ---- MCMC ----
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    Rcpp::checkUserInterrupt();
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }

    // -- Update trees --
    ctr->R += (ctr->Rmat).col(0);
    ctr->fhat.setZero();
    ctr->totTerm = 0.0; ctr->sumTermT2 = 0.0;
    // ctr->exDLM.setZero();
    ctr->phiMH = 0; ctr->phiMHNew = 0;
    for (t = 0; t < ctr->nTrees; t++) {
      dlmtreeTDLMFixed_Gaussian_TreeMCMC(t, fixedTrees[t], ctr, dgn, Exp);
      ctr->fhat += (ctr->Rmat).col(t);

      if (t < ctr->nTrees - 1)
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
    }

    // -- Update model --
    ctr->R = ctr->Y - ctr->fhat;
    tdlmModelEst(ctr);
    xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->nu)));
    ctr->nu = 1.0 / R::rgamma(0.5 * ctr->totTerm + 0.5,
                               1.0 / (0.5 * ctr->sumTermT2 / (ctr->sigma2) +
                                      xiInv));

    if ((ctr->sigma2 != ctr->sigma2) || (ctr->nu != ctr->nu)) {
      Rcout << "\n" << ctr->sigma2 << "\n" <<
        ctr->nu << "\n" << ctr->totTerm << "\n" << ctr->sumTermT2 << "\n" <<
        ctr->xiInvSigma2 << "\n" << xiInv << "\n" <<
        ctr->tau;
      stop("\nNaN values occured during model run, rerun model.\n");
    }
    
    
    
    // -- Record --
    if (ctr->record > 0) {
      (dgn->gamma).col(ctr->record - 1) = ctr->gamma;
      (dgn->sigma2)(ctr->record - 1) = ctr->sigma2;
      (dgn->nu)(ctr->record - 1) = ctr->nu;
      (dgn->tau).col(ctr->record - 1) = ctr->tau;
      (dgn->phi)(ctr->record - 1) = ctr->phi;
      dgn->fhat += ctr->fhat;
      // dlmtreeRecDLM(ctr, dgn);
    } // end record



    // -- Progress --
    prog->printMark();
  } // end MCMC


  // -- Prepare outout --  
  Eigen::MatrixXd TreeStructs((dgn->DLMexp).size(), 9);
  for (s = 0; s < (dgn->DLMexp).size(); ++s)
    TreeStructs.row(s) = dgn->DLMexp[s];

  Eigen::VectorXd sigma2 = dgn->sigma2;
  Eigen::VectorXd nu = dgn->nu;
  Eigen::MatrixXd tau = (dgn->tau).transpose();
  Eigen::VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  Eigen::MatrixXd gamma = (dgn->gamma).transpose();
  Eigen::VectorXd phi = dgn->phi;

  delete prog;
  delete ctr;
  delete dgn;
  for (s = 0; s < fixedTrees.size(); ++s) {
    for (t = 0; t < fixedTrees[s].size(); ++t)
      delete fixedTrees[s][t];
  }

  return(Rcpp::List::create(Named("TreeStructs") = wrap(TreeStructs),
                            Named("fhat") = wrap(fhat),
                            Named("sigma2") = wrap(sigma2),
                            Named("nu") = wrap(nu),
                            Named("tau") = wrap(tau),
                            Named("gamma") = wrap(gamma),
                            Named("phi") = wrap(phi)));

} // end dlmtreeTDLMGaussian



void dlmtreeTDLMFixed_Gaussian_TreeMCMC(
  int t, std::vector<Node*> fixedNodes,
  dlmtreeCtr* ctr, dlmtreeLog *dgn, exposureDat* Exp)
{
  int step;
  int success = 0;
  double stepMhr = 0;
  double ratio = 0;
  double treevar = (ctr->nu) * (ctr->tau)(t);
  std::size_t s;
  std::vector<Node*> dlmTerm;
  Eigen::VectorXd ZtR = (ctr->Z).transpose() * (ctr->R);
  treeMHR mhr0 = dlmtreeTDLMFixedMHR(fixedNodes, ctr, ZtR, treevar, 0);
  double RtR = (ctr->R).dot(ctr->R);
  double RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);
  treeMHR mhr;
  
  // -- Propose new nested tree at each modifier node --
  for (Node* tn : fixedNodes) {
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
      // if (ctr->nSplits == 0)
        mhr = dlmtreeTDLMFixedMHR(fixedNodes, ctr, ZtR, treevar, 1);
      // else
      //   mhr = dlmtreeTDLNMNested_MHR(modTerm, ctr, ZtR, treevar, 1);
      ratio = calcLogRatioFixedTDLM(mhr0, mhr, RtR, RtZVgZtR, ctr, 
                                    stepMhr, treevar);
      // Rcout << " ratioTT=" << ratio;
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
  double xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->tau)(t)));
  (ctr->tau)(t) = 1.0 / R::rgamma(0.5 * mhr0.draw.size() + 0.5,
                                  1.0 / ((0.5 * mhr0.termT2 /
                                          (ctr->sigma2 * ctr->nu)) + xiInv));
  ctr->Rmat.col(t) = mhr0.fitted;
  ctr->sumTermT2 += mhr0.termT2 / (ctr->tau(t));
  ctr->totTerm += static_cast<double>(mhr0.draw.size());
  
  // -- Record --
  if (ctr->record > 0) {
    
    // -- Update DLM partial estimate --
    std::string rule;
    Eigen::VectorXd rec(9);
    Eigen::VectorXd draw(ctr->pX);
    Node* tn;
    std::vector<Node*> nested;
    int drawIdx = 0;
    
    for (s = 0; s < fixedNodes.size(); ++s) {
      tn = fixedNodes[s];
      nested = tn->nodevals->nestedTree->listTerminal();
      
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
} // end dlmtreeTDLM_Gaussian_TreeMCMC function



// function to calculate part of MH ratio
treeMHR dlmtreeTDLMFixedMHR(std::vector<Node*> fixedNodes,
                   dlmtreeCtr* ctr, 
                   Eigen::VectorXd ZtR,
                   double treevar,
                   double updateNested)
{
  std::size_t s, s2;
  treeMHR out;
  int pXMod = int(fixedNodes.size());
  int totTerm = 0;
  std::vector<std::vector<Node*> > nestedTerm;
  for (s = 0; s < fixedNodes.size(); ++s) {
    nestedTerm.push_back(
      fixedNodes[s]->nodevals->nestedTree->listTerminal(updateNested));
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
  
  for (s = 0; s < fixedNodes.size(); ++s) {
    pX = nestedTerm[s].size();
    Eigen::MatrixXd temp(ctr->n, pX); temp.setZero();
    X.push_back(temp);
    
    // create nested tree data matrices
    for (s2 = 0; s2 < nestedTerm[s].size(); ++s2) {
      X[s].col(s2) = nestedTerm[s][s2]->nodevals->X;
    } // end loop over nested tree nodes
    
    if (fixedNodes[s]->nodevals->updateXmat) {                                  
      Xtemp.resize(fixedNodes[s]->nodevals->idx.size(), pX); 
      Xtemp.setZero();
      Ztemp.resize(fixedNodes[s]->nodevals->idx.size(), ctr->pZ); 
      Ztemp.setZero();
      Rtemp.resize(fixedNodes[s]->nodevals->idx.size()); 
      Rtemp.setZero();
      
      j = 0;
      for (int i : fixedNodes[s]->nodevals->idx) {
        Xtemp.row(j) = X[s].row(i);
        Ztemp.row(j) = ctr->Z.row(i);
        Rtemp(j) = ctr->R(i);
        ++j;
      }
      
      fixedNodes[s]->nodevals->XtX.resize(pX, pX);
      fixedNodes[s]->nodevals->XtX = Xtemp.transpose() * Xtemp;
      fixedNodes[s]->nodevals->ZtXmat.resize(ctr->pZ, pX);
      fixedNodes[s]->nodevals->ZtXmat = Ztemp.transpose() * Xtemp;
      fixedNodes[s]->nodevals->VgZtXmat.resize(ctr->pZ, pX);
      fixedNodes[s]->nodevals->VgZtXmat = ctr->Vg * fixedNodes[s]->nodevals->ZtXmat;
      fixedNodes[s]->nodevals->updateXmat = 0;
      XtR.segment(start, pX) = Xtemp.transpose() * Rtemp;
      
    } else { // reuse precalculated matrices
      
      for (int i : fixedNodes[s]->nodevals->idx) {
        XtR.segment(start, pX).noalias() += X[s].row(i).transpose() * ctr->R(i);
      }
      
    } // end update XtX and ZtX matrices    
    
    LInv.resize(pX, pX); LInv.setZero();
    LInv.diagonal().array() += 1.0 / treevar;
    XXiblock.block(start, start, pX, pX) =  
      (fixedNodes[s]->nodevals->XtX + LInv).inverse();
    ZtX.block(0, start, ctr->pZ, pX) = fixedNodes[s]->nodevals->ZtXmat;
    VgZtX.block(0, start, ctr->pZ, pX) = fixedNodes[s]->nodevals->VgZtXmat;
    
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
  for (s = 0; s < fixedNodes.size(); ++s) {    
    pX = nestedTerm[s].size();
    drawTemp.resize(pX);    
    drawTemp = out.draw.segment(start, pX);
    
    for (int i : fixedNodes[s]->nodevals->idx)
      out.fitted(i) = X[s].row(i) * drawTemp;
      
    start += pX;
  } // end loop over modifier nodes to calculate dr
  
  out.beta = ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();
  out.termT2 = (out.draw).dot(out.draw);
  out.totTerm = double(totTerm); 
  out.nModTerm = double(pXMod);
  return(out);
} // end dlmtreeMHR function

double calcLogRatioFixedTDLM(treeMHR mhr0, treeMHR mhr, double RtR, 
                             double RtZVgZtR, dlmtreeCtr* ctr, 
                             double stepMhr, double treevar)
{
  return(stepMhr +
         mhr.logVThetaChol - mhr0.logVThetaChol -
         (0.5 * (ctr->n + 1.0) *
           (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
           log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
         (log(treevar) * 0.5 * round(mhr.totTerm - mhr0.totTerm)));
}