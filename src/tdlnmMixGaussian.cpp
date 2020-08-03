#include <RcppEigen.h>
#include "tdlnmCtr.h"
#include "exposureDat.h"
#include "Fncs.h"
using namespace Rcpp;

void tdlnmMixGaussianTreeMCMC(int, Node*, Node*, tdlnmCtr*, tdlnmLog*,
                              std::vector<exposureDat*>);
double tdlnmProposeTree(Node *tree, exposureDat *Exp, tdlnmCtr *ctr, int step);

struct treeMixMHR {
public:
  Eigen::MatrixXd Xd;
  Eigen::MatrixXd tempV;
  Eigen::VectorXd draw1, draw2, drawMix, drawAll;
  double logVThetaChol, beta;
  double term1T2, term2T2, mixT2, nTerm1, nTerm2;
  int pXd;
};

treeMixMHR mixMHR(std::vector<Node*>, std::vector<Node*>, tdlnmCtr*,
                  Eigen::VectorXd, double, double, double, double,
                  Node*, bool);

// Rcpp::List tdlnmMixGaussian(Rcpp::List);

// [[Rcpp::export]]
Rcpp::List tdlnmMixGaussian(const Rcpp::List model)
{

  int i, j, k, t;
  // ---- Set up general control variables ----
  tdlnmCtr *ctr = new tdlnmCtr;
  ctr->iter = as<int>(model["nIter"]);
  ctr->burn = as<int>(model["nBurn"]);
  ctr->thin = as<int>(model["nThin"]);
  ctr->nRec = floor(ctr->iter / ctr->thin);
  ctr->nTrees = as<int>(model["nTrees"]);
  ctr->Y = as<Eigen::VectorXd>(model["Y"]);
  ctr->n = (ctr->Y).size();

  ctr->Z = as<Eigen::MatrixXd>(model["Z"]);
  ctr->pZ = (ctr->Z).cols();
  Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ);
  VgInv = (ctr->Z).transpose() * (ctr->Z);
  VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg = VgInv.inverse();
  VgInv.resize(0,0);
  ctr->VgChol = (ctr->Vg).llt().matrixL();

  // ---- Pre-calculate single node tree matrices ----
  std::vector<exposureDat*> Exp;
  Rcpp::List exp_dat = as<Rcpp::List>(model["X"]);
  ctr->nExp = exp_dat.size();
  for (i = 0; i < ctr->nExp; ++i) {
    Exp.push_back(new exposureDat(
        as<Eigen::MatrixXd>(as<Rcpp::List>(exp_dat[i])["Tcalc"]), ctr->Z, ctr->Vg));
  }
  ctr->pX = Exp[0]->pX;
  ctr->nSplits = 0;
  ctr->interaction = as<int>(model["interaction"]);
  ctr->nMix = 0;
  if (ctr->interaction) {
    ctr->nMix += int (ctr->nExp * (ctr->nExp - 1) / 2);
    if (ctr->interaction == 2) {
      ctr->nMix += ctr->nExp;
    }
  }
  // ctr->modZeta = 1;
  ctr->modKappa = double(ctr->nExp) / double(2 * ctr->nTrees);


  // ---- Logs ----;
  tdlnmLog *dgn = new tdlnmLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec); (dgn->gamma).setZero();
  (dgn->sigma2).resize(ctr->nRec); (dgn->sigma2).setZero();
  (dgn->nu).resize(ctr->nRec); (dgn->nu).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec); (dgn->tau).setZero();
  (dgn->muExp).resize(ctr->nExp, ctr->nRec); (dgn->muExp).setZero();
  if (ctr->interaction > 0) {
    (dgn->muMix).resize(ctr->nMix, ctr->nRec); (dgn->muMix).setZero();
  } else {
    (dgn->muMix).resize(1, 1); (dgn->muMix).setZero();
  }
  (dgn->expProb).resize(ctr->nExp, ctr->nRec); (dgn->expProb).setZero();
  (dgn->fhat).resize(ctr->n); (dgn->fhat).setZero();
  (dgn->termNodes).resize(ctr->nTrees, ctr->nRec); (dgn->termNodes).setZero();
  (dgn->termNodes2).resize(ctr->nTrees, ctr->nRec); (dgn->termNodes2).setZero();
  (dgn->tree1Exp).resize(ctr->nTrees, ctr->nRec); (dgn->tree1Exp).setZero();
  (dgn->tree2Exp).resize(ctr->nTrees, ctr->nRec); (dgn->tree2Exp).setZero();



  // ---- Initial draws ----
  (ctr->fhat).resize(ctr->n); (ctr->fhat).setZero();
  ctr->R = ctr->Y;
  (ctr->gamma).resize(ctr->pZ);
  (ctr->totTermExp).resize(ctr->nExp); (ctr->totTermExp).setZero();
  (ctr->sumTermT2Exp).resize(ctr->nExp); (ctr->sumTermT2Exp).setZero();
  (ctr->muExp).resize(ctr->nExp); (ctr->muExp).setOnes();
  if (ctr->interaction) {
    (ctr->totTermMix).resize(ctr->nExp, ctr->nExp); (ctr->totTermMix).setZero();
    (ctr->sumTermT2Mix).resize(ctr->nExp, ctr->nExp); (ctr->sumTermT2Mix).setZero();
    (ctr->muMix).resize(ctr->nExp, ctr->nExp); (ctr->muMix).setOnes();
  }
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  ctr->nu = 1; // Need to define for first update of sigma2
  ctr->sigma2 = 1;
  tdlnmModelEst(ctr);
  double xiInv = R::rgamma(1, 0.5);
  ctr->nu = 1 / R::rgamma(0.5, 1 / xiInv);
  (ctr->tau).resize(ctr->nTrees);
  for (t = 0; t < ctr->nTrees; t++) {
    xiInv = R::rgamma(1, 0.5);
    (ctr->tau)(t) = 1 / R::rgamma(0.5, 1 / xiInv);
  }
  ctr->nTerm.resize(ctr->nTrees); (ctr->nTerm).setOnes();
  ctr->nTerm2.resize(ctr->nTrees); (ctr->nTerm2).setOnes();
  (ctr->Rmat).resize(ctr->n, ctr->nTrees); (ctr->Rmat).setZero();





  // ---- Create trees ----
  (ctr->tree1Exp).resize(ctr->nTrees); (ctr->tree1Exp).setZero();
  (ctr->tree2Exp).resize(ctr->nTrees); (ctr->tree2Exp).setZero();
  ctr->expProb = as<Eigen::VectorXd>(model["expProb"]);
  (ctr->expCount).resize((ctr->expProb).size());
  std::vector<Node*> trees1;
  std::vector<Node*> trees2;
  NodeStruct *ns;
  ns = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                      as<Eigen::VectorXd>(model["splitProb"]),
                      as<Eigen::VectorXd>(model["timeProb"]));
  for (t = 0; t < ctr->nTrees; t++) {
    ctr->tree1Exp(t) = sampleInt(ctr->expProb);
    trees1.push_back(new Node(0, 1));
    trees1[t]->nodestruct = ns->clone();
    Exp[ctr->tree1Exp(t)]->updateNodeVals(trees1[t]);

    ctr->tree2Exp(t) = sampleInt(ctr->expProb);
    trees2.push_back(new Node(0, 1));
    trees2[t]->nodestruct = ns->clone();
    Exp[ctr->tree2Exp(t)]->updateNodeVals(trees2[t]);
  }
  delete ns;





  time_t startTime = time(NULL); double timediff;
  ctr->verbose = bool (model["verbose"]);
  ctr->diagnostics = bool (model["diagnostics"]);
  ctr->stepProb = as<std::vector<double>>(model["stepProb"]);
  ctr->treePrior = as<std::vector<double>>(model["treePrior"]);



  // if (ctr->sigma2 != ctr->sigma2) {
  //   int nAttempt = as<int>(model["try"]);
  //   if (nAttempt < 5) {
  //     model["try"] = nAttempt + 1;
  //     if (ctr->verbose)
  //       Rcout << "NaN values in initial draws, restarting model...attempt " <<
  //         nAttempt + 1 << "/5\n";
  //     std::this_thread::sleep_until(std::chrono::system_clock::now() +
  //       std::chrono::seconds(nAttempt + 1));
  //     return(tdlnmMixGaussian(model));
  //   } else {
  //     stop("Model run failed");
  //   }
  // }



  if (ctr->verbose)
    Rcout << "Running model '.' = 100 iterations\n";

  // ---- MCMC ----
  double sigmanu;
  std::size_t s;
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }


    // -- Update trees --
    ctr->R += (ctr->Rmat).col(0);
    (ctr->fhat).setZero();
    ctr->totTerm = 0;
    ctr->sumTermT2 = 0;
    (ctr->totTermExp).setZero();
    (ctr->sumTermT2Exp).setZero();
    (ctr->expCount).setZero();
    if (ctr->interaction > 0) {
      (ctr->totTermMix).setZero(); (ctr->sumTermT2Mix).setZero();
    }

    for (t = 0; t < ctr->nTrees; t++) {
      tdlnmMixGaussianTreeMCMC(t, trees1[t], trees2[t], ctr, dgn, Exp);
      ctr->fhat += (ctr->Rmat).col(t);
      ctr->R -= (ctr->Rmat).col(t);
      if (t < ctr->nTrees - 1) {
        ctr->R += (ctr->Rmat).col(t + 1);
      }
    }


    // -- Pre-calculations for control and variance --
    ctr->R = ctr->Y - ctr->fhat;
    ctr->sumTermT2 = (ctr->sumTermT2Exp).sum();
    ctr->totTerm = (ctr->totTermExp).sum();
    if(ctr->interaction) {
      ctr->sumTermT2 += (ctr->sumTermT2Mix).sum();
      ctr->totTerm += (ctr->totTermMix).sum();
    }



    // -- Update control --
    tdlnmModelEst(ctr);

    if ((ctr->sigma2 != ctr->sigma2)) {
      // int nAttempt = as<int>(model["try"]);
      // if ((nAttempt < 5) && (ctr->b < ctr->burn)) {
      //   model["try"] = nAttempt + 1;
      //   if (ctr->verbose)
      //     Rcout << "NaN values occured during model run, restarting model...attempt " <<
      //       nAttempt + 1 << "/5\n";
      //   std::this_thread::sleep_until(std::chrono::system_clock::now() +
      //     std::chrono::seconds(nAttempt + 1));
      //   return(tdlnmMixGaussian(model));
      // } else {
      stop("NaN values occured during model run, stopping.");
      // }
    }


    xiInv = R::rgamma(1, (ctr->nu) / (ctr->nu + 1.0));
    ctr->nu = 1.0 / R::rgamma(0.5 * ctr->totTerm + 0.5,
                            1.0 / (ctr->sumTermT2 / (2.0 * ctr->sigma2) + xiInv));
    sigmanu = ctr->sigma2 * ctr->nu;
    for (i = 0; i < ctr->nExp; ++i) {
      xiInv = R::rgamma(1, (ctr->muExp(i)) / (ctr->muExp(i) + 1.0));
      ctr->muExp(i) = 1.0 / R::rgamma(0.5 * ctr->totTermExp(i) + 0.5,
                 1.0 / (ctr->sumTermT2Exp(i) / (2.0 * sigmanu) + xiInv));
      if (ctr->interaction) {
        for (j = i; j < ctr->nExp; ++j) {
          if ((j > i) || (ctr->interaction == 2)) {
            xiInv = R::rgamma(1, (ctr->muMix(j, i)) / (ctr->muMix(j, i) + 1.0));
            ctr->muMix(j, i) = 1.0 / R::rgamma(0.5 * ctr->totTermMix(j, i) + 0.5,
                       1.0 / (ctr->sumTermT2Mix(j, i) / (2.0 * sigmanu) + xiInv));
            if (ctr->muMix(j, i) != ctr->muMix(j, i)) {
              stop("\nNaN values occured during model run, rerun model.\n");
            }
          }
        }
      }
    }


    // -- Update exposure selection probability
    if ((ctr->b > 1000) || (ctr->b > (0.5 * ctr->burn))) {
      // double beta = R::rbeta(ctr->modZeta, 1);
      // double newModKappa = beta * ctr->nExp / (1 - beta);
      //
      // double dirichRatio =
      //   logDirichletDensity(ctr->expProb,
      //                       ((ctr->expCount).array() +
      //                         (newModKappa / ctr->nExp)).matrix()) -
      //   logDirichletDensity(ctr->expProb,
      //                       ((ctr->expCount).array() +
      //                        (ctr->modKappa / ctr->nExp)).matrix());
      //
      // if (log(R::runif(0, 1)) < dirichRatio) {
      //   ctr->modKappa = newModKappa;
      // }

      ctr->expProb = rDirichlet(((ctr->expCount).array() + ctr->modKappa).matrix());
    }





    // -- Record --
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
      (dgn->muExp).col(ctr->record - 1) = ctr->muExp;
      if (ctr->interaction) {
        k = 0;
        for (i = 0; i < ctr->nExp; ++i) {
          for (j = i; j < ctr->nExp; ++j) {
            if ((j > i) || (ctr->interaction == 2)) {
              dgn->muMix(k, ctr->record - 1) = ctr->muMix(j, i);
              ++k;
            }
          }
        }
      }
    }

    // -- Progress --
    if (ctr->verbose) {
      if (ctr->b % 100 == 0)
        Rcout << ".";
      if ((ctr->b > ctr->burn) && ((ctr->b - ctr->burn) % 5000 == 0))
        Rcout << "5000\n";
      if (ctr->b == ctr->burn) {
        timediff = difftime(time(NULL), startTime);
        if (timediff > 3600) {
          Rprintf("\nBurn-in time: %.2g hours", round(100 * timediff / 3600) / 100);
        } else if (timediff > 60) {
          Rprintf("\nBurn-in time: %.2g minutes", round(100 * timediff / 60) / 100);
        } else {
          Rprintf("\nBurn-in time: %.2g seconds", round(100 * timediff) / 100);
        }

        timediff = timediff * ctr->iter / ctr->burn;
        if (timediff > 3600) {
          Rprintf("\nEstimated time to completion: %.2g hours\n",
                  round(100 * timediff / 3600) / 100);
        } else if (timediff > 60) {
          Rprintf("\nEstimated time to completion: %.2g minutes\n",
                  round(100 * timediff / 60) / 100);
        } else {
          Rprintf("\nEstimated time to completion: %.2g seconds\n", round(100 * timediff) / 100);
        }
      }
    }
  } // -- End MCMC --




  Eigen::MatrixXd DLM((dgn->DLMexp).size(), 7);
  i = 0;
  for (s = 0; s < (dgn->DLMexp).size(); ++s) {
    DLM.row(i) = dgn->DLMexp[s];
    ++i;
  }

  Eigen::VectorXd sigma2 = dgn->sigma2;
  Eigen::VectorXd nu = dgn->nu;
  Eigen::VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  Eigen::MatrixXd gamma = (dgn->gamma).transpose();
  Eigen::MatrixXd tau = (dgn->tau).transpose();
  Eigen::MatrixXd termNodes = (dgn->termNodes).transpose();
  Eigen::MatrixXd termNodes2 = (dgn->termNodes2).transpose();
  Eigen::MatrixXd tree1Exp = (dgn->tree1Exp).transpose();
  Eigen::MatrixXd tree2Exp = (dgn->tree2Exp).transpose();
  Eigen::MatrixXd expProb = (dgn->expProb).transpose();
  Eigen::MatrixXd muExp = (dgn->muExp).transpose();
  Eigen::MatrixXd muMix(1, 1); muMix.setZero();
  Eigen::MatrixXd MIX(0, 10); MIX.setZero();
  if (ctr->interaction) {
    muMix.resize((dgn->muMix).cols(), (dgn->muMix).rows());
    muMix = (dgn->muMix).transpose();
    MIX.resize((dgn->MIXexp).size(), 10);
    i = 0;
    for (s = 0; s < (dgn->MIXexp).size(); ++s) {
      MIX.row(i) = dgn->MIXexp[s];
      ++i;
    }
  }

  delete ctr;
  delete dgn;
  for (s = 0; s < Exp.size(); ++s) {
    delete Exp[s];
  }
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
                            Named("muExp") = wrap(muExp),
                            Named("muMix") = wrap(muMix)));
}






void tdlnmMixGaussianTreeMCMC(int t, Node *tree1, Node *tree2,
                              tdlnmCtr *ctr, tdlnmLog *dgn,
                              std::vector<exposureDat*> Exp)
{

  int m1, m2, newExp, success, step1, step2;
  double RtR, RtZVgZtR, stepMhr, ratio;
  double m1Var, m2Var, mixVar, newExpVar, newMixVar, treeVar;
  size_t s;
  std::vector<Node*> term1, term2, newTerm;
  Node* newTree = 0;
  treeMixMHR mhr0, mhr;

  term1 = tree1->listTerminal();
  term2 = tree2->listTerminal();
  treeVar = (ctr->nu) * (ctr->tau[t]);
  m1 = ctr->tree1Exp[t];
  m2 = ctr->tree2Exp[t];
  m1Var = ctr->muExp(m1);
  m2Var = ctr->muExp(m2);
  mixVar = 0;
  if ((ctr->interaction) && ((ctr->interaction == 2) || (m1 != m2))) {
    if (m1 <= m2) {
      mixVar = ctr->muMix(m2, m1);
    } else {
      mixVar = ctr->muMix(m1, m2);
    }
  }
  Eigen::VectorXd ZtR = (ctr->Z).transpose() * (ctr->R);
  RtR = (ctr->R).dot(ctr->R);
  RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);


  // -- Update tree 1 --
  newExp = m1;
  newExpVar = m1Var;
  newMixVar = mixVar;
  stepMhr = 0;
  success = 0;


  // List current tree terminal nodes
  step1 = sampleInt(ctr->stepProb, 1);
  if ((term1.size() == 1) && (step1 < 3)) {
    step1 = 0;
  }

  // propose update
  if (step1 < 3) {
    stepMhr = tdlnmProposeTree(tree1, Exp[m1], ctr, step1);
    success = tree1->isProposed();
    newTerm = tree1->listTerminal(1);

    // switch exposures
  } else {
    newExp = sampleInt(ctr->expProb);
    if (newExp != m1) {
      success = 1;
      stepMhr = log(ctr->expProb(m1)) - log(ctr->expProb(newExp));
      newExpVar = ctr->muExp(newExp);

      newTree = new Node(*tree1);
      newTree->setUpdate(1);
      newTerm = newTree->listTerminal();
      for (Node* nt : newTerm)
        Exp[newExp]->updateNodeVals(nt);


      if ((ctr->interaction) && ((ctr->interaction == 2) || (newExp != m2))) {
        if (newExp <= m2) {
          newMixVar = ctr->muMix(m2, newExp);
        } else {
          newMixVar = ctr->muMix(newExp, m2);
        }
      } else {
        newMixVar = 0;
      }
    }
  }


  // -- Tree 1 MHR --
  if ((tree1->nodevals->tempV).rows() == 0) {
    mhr0 = mixMHR(term1, term2, ctr, ZtR, treeVar, m1Var, m2Var, mixVar, tree1, 1);
  } else {
    mhr0 = mixMHR(term1, term2, ctr, ZtR, treeVar, m1Var, m2Var, mixVar, tree1, 0);
  }

  if (success) {
    mhr = mixMHR(newTerm, term2, ctr, ZtR, treeVar, newExpVar, m2Var, newMixVar,
                 tree1, 1);

    ratio =
      stepMhr +
      mhr.logVThetaChol - mhr0.logVThetaChol -
      (0.5 * (ctr->n + 1.0) *
      (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
      log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
      (0.5 * ((log(treeVar * newExpVar) * mhr.nTerm1) -
      (log(treeVar * m1Var) * mhr0.nTerm1)));

    if (newMixVar != 0) {
      ratio -= 0.5 * log(treeVar * newMixVar) * mhr.nTerm1 * mhr0.nTerm2;
    }
    if (mixVar != 0) {
      ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2;
    }

    if (log(R::runif(0, 1)) < ratio) {
      mhr0 = mhr;
      success = 2;

      if (step1 == 3) {
        m1 = newExp;
        m1Var = newExpVar;
        mixVar = newMixVar;
        tree1->replaceNodeVals(newTree);
      } else {
        tree1->accept();
      }
      (tree1->nodevals->tempV).resize(mhr0.pXd, mhr0.pXd);
      tree1->nodevals->tempV = mhr0.tempV;
      term1 = tree1->listTerminal();

    } else {
      tree1->reject();
    }

  } else {
    if (step1 < 3)
      tree1->reject();
  }

  if (newTree != 0)
    delete newTree;
  newTree = 0;



  // -- Update tree 2 --
  newExp = m2;
  newExpVar = m2Var;
  newMixVar = mixVar;
  stepMhr = 0;
  success = 0;


  // List current tree terminal nodes
  step2 = sampleInt(ctr->stepProb, 1);
  if ((term2.size() == 1) && (step2 < 3)) {
    step2 = 0;
  }


  // propose update
  if (step2 < 3) {
    stepMhr = tdlnmProposeTree(tree2, Exp[m2], ctr, step2);
    success = tree2->isProposed();
    newTerm = tree2->listTerminal(1);

    // switch exposures
  } else {
    newExp = sampleInt(ctr->expProb);
    if (newExp != m2) {
      success = 1;
      stepMhr = log(ctr->expProb(m2)) - log(ctr->expProb(newExp));
      newExpVar = ctr->muExp(newExp);

      newTree = new Node(*tree2);
      newTree->setUpdate(1);
      newTerm = newTree->listTerminal();
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
    mhr = mixMHR(term1, newTerm, ctr, ZtR, treeVar, m1Var, newExpVar, newMixVar,
                 tree1, 1);


    ratio =
      stepMhr +
      mhr.logVThetaChol - mhr0.logVThetaChol -
      (0.5 * (ctr->n + 1.0) *
      (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
      log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
      (0.5 * ((log(treeVar * newExpVar) * mhr.nTerm2) -
      (log(treeVar * m2Var) * mhr0.nTerm2)));

    if (newMixVar != 0) {
      ratio -= 0.5 * log(treeVar * newMixVar) * mhr0.nTerm1 * mhr.nTerm2;
    }
    if (mixVar != 0) {
      ratio += 0.5 * log(treeVar * mixVar) * mhr0.nTerm1 * mhr0.nTerm2;
    }

    if (log(R::runif(0, 1)) < ratio) {

      mhr0 = mhr;
      success = 2;

      if (step2 == 3) {
        m2 = newExp;
        m2Var = newExpVar;
        mixVar = newMixVar;
        tree2->replaceNodeVals(newTree);

      } else {
        tree2->accept();
      }
      (tree1->nodevals->tempV).resize(mhr0.pXd, mhr0.pXd);
      tree1->nodevals->tempV = mhr0.tempV;
      term2 = tree2->listTerminal();

    } else {
      tree2->reject();
    }

  } else {
    tree2->reject();
  }

  if (newTree != 0)
    delete newTree;
  newTree = 0;


  // -- Update variance and residuals --
  double xiInv = R::rgamma(1, (ctr->tau[t]) / (ctr->tau[t] + 1.0));
  double tauT2 = mhr0.term1T2 / m1Var + mhr0.term2T2 / m2Var;
  int totTerm = mhr0.nTerm1 + mhr0.nTerm2;
  if (mixVar != 0) {
    tauT2 += mhr0.mixT2 / mixVar;
    totTerm += mhr0.nTerm1 * mhr0.nTerm2;
  }
  ctr->tau[t] = 1.0 / R::rgamma(0.5 * totTerm + 0.5,
                              1.0 / (tauT2 / (2.0 * ctr->sigma2 * ctr->nu) + xiInv));
  (ctr->nTerm)(t) = mhr0.nTerm1;
  (ctr->nTerm2)(t) = mhr0.nTerm2;
  (ctr->tree1Exp)(t) = m1;
  (ctr->tree2Exp)(t) = m2;
  (ctr->expCount)(m1)++;
  (ctr->expCount)(m2)++;
  (ctr->totTermExp)(m1) += mhr0.nTerm1;
  (ctr->totTermExp)(m2) += mhr0.nTerm2;
  (ctr->sumTermT2Exp)(m1) += mhr0.term1T2 / ctr->tau[t];
  (ctr->sumTermT2Exp)(m2) += mhr0.term2T2 / ctr->tau[t];
  if (mixVar != 0) {
    if (m1 <= m2) {
      (ctr->totTermMix)(m2, m1) += mhr0.nTerm1 * mhr0.nTerm2;
      (ctr->sumTermT2Mix)(m2, m1) += mhr0.mixT2 / ctr->tau[t];
    } else {
      (ctr->totTermMix)(m1, m2) += mhr0.nTerm1 * mhr0.nTerm2;
      (ctr->sumTermT2Mix)(m1, m2) += mhr0.mixT2 / ctr->tau[t];
    }
  }

  (ctr->Rmat).col(t) = mhr0.Xd * mhr0.drawAll;


  // Record
  if (ctr->record > 0) {
    Eigen::VectorXd rec(7);
    Eigen::VectorXd mix(10);
    rec << ctr->record, t, 0, 0, 0, 0, 0;
    mix << ctr->record, t, 0, 0, 0, 0, 0, 0, 0, 0;
    int k = 0;
    for(int i = 0; i < mhr0.nTerm1; ++i) {
      rec[2] = m1;
      rec[3] = (term1[i]->nodestruct)->get(3);
      rec[4] = (term1[i]->nodestruct)->get(4);
      rec[5] = mhr0.draw1(i);
      (dgn->DLMexp).push_back(rec);
      for (int j = 0; j < mhr0.nTerm2; ++j) {
        if (i == 0) {
          rec[2] = m2;
          rec[3] = (term2[j]->nodestruct)->get(3);
          rec[4] = (term2[j]->nodestruct)->get(4);
          rec[5] = mhr0.draw2(j);
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
}


treeMixMHR mixMHR(std::vector<Node*> nodes1, std::vector<Node*> nodes2,
                  tdlnmCtr *ctr, Eigen::VectorXd ZtR,
                  double treeVar, double m1Var, double m2Var, double mixVar,
                  Node* tree, bool newTree)
{
  treeMixMHR out;
  int pX1 = nodes1.size();
  int pX2 = nodes2.size();
  int pXd = pX1 + pX2;
  int interaction = 0;
  if (mixVar != 0) {
    pXd += pX1 * pX2;
    interaction = 1;
  }

  out.Xd.resize(ctr->n, pXd); out.Xd.setZero();
  Eigen::MatrixXd ZtX(ctr->pZ, pXd); ZtX.setZero();
  Eigen::VectorXd diagVar(pXd); diagVar.setZero();

  int i, j, k;
  for (i = 0; i < pX1; ++i) {
    out.Xd.col(i) = (nodes1[i]->nodevals)->X;
    ZtX.col(i) = (nodes1[i]->nodevals)->ZtX;
    diagVar(i) = 1.0 / (m1Var * treeVar);
  }

  for (j = 0; j < pX2; ++j) {
    k = pX1 + j;
    out.Xd.col(k) = (nodes2[j]->nodevals)->X;
    ZtX.col(k) = (nodes2[j]->nodevals)->ZtX;
    diagVar(k) = 1.0 / (m2Var * treeVar);
  }

  if (interaction) {
    for (i = 0; i < pX1; ++i) {
      for (j = 0; j < pX2; ++j) {
        k = pX1 + pX2 + i * pX2 + j;
        out.Xd.col(k) =
          (((nodes1[i]->nodevals)->X).array() *
          ((nodes2[j]->nodevals)->X).array()).matrix();
        ZtX.col(k) = (ctr->Z).transpose() * out.Xd.col(k);
        diagVar(k) = 1.0 / (mixVar * treeVar);
      }
    }
  }

  // calculate MHR
  const Eigen::MatrixXd VgZtX = ctr->Vg * ZtX;
  Eigen::MatrixXd tempV(pXd, pXd);
  if (newTree) {
    tempV = out.Xd.transpose() * out.Xd - ZtX.transpose() * VgZtX;
    out.tempV = tempV;
  } else {
    tempV = tree->nodevals->tempV;
  }
  tempV.diagonal() += diagVar;

  Eigen::MatrixXd VTheta(pXd, pXd);
  VTheta.triangularView<Eigen::Lower>() =
    tempV.selfadjointView<Eigen::Lower>().llt().solve(
        Eigen::MatrixXd::Identity(pXd, pXd));
  Eigen::VectorXd XtVzInvR = out.Xd.transpose() * ctr->R;
  XtVzInvR.noalias() -= VgZtX.transpose() * ZtR;
  const Eigen::VectorXd ThetaHat =
    VTheta.selfadjointView<Eigen::Lower>() * XtVzInvR;
  const Eigen::MatrixXd VThetaChol =
    VTheta.selfadjointView<Eigen::Lower>().llt().matrixL();


  Eigen::VectorXd ThetaDraw = ThetaHat;
  ThetaDraw.noalias() += VThetaChol * as<Eigen::VectorXd>(rnorm(pXd, 0, sqrt(ctr->sigma2)));

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

