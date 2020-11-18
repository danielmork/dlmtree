#include <RcppEigen.h>
#include "tdlmCtr.h"
#include "exposureDat.h"
#include "Fncs.h"
#include "omp.h"
using namespace Rcpp;

void tdlmmBinomialTreeMCMC(int, Node*, Node*, tdlmCtr*, tdlmLog*,
                              std::vector<exposureDat*>);
double tdlmProposeTree(Node *tree, exposureDat *Exp, tdlmCtr *ctr, int step);

struct treeMixBinomMHR {
public:
  Eigen::MatrixXd Xd;
  Eigen::MatrixXd tempV;
  Eigen::VectorXd draw1, draw2, drawMix, drawAll;
  double logVThetaChol, beta;
  double term1T2, term2T2, mixT2, nTerm1, nTerm2;
  int pXd;
};

treeMixBinomMHR mixBinomMHR(std::vector<Node*>, std::vector<Node*>, tdlmCtr*,
                  Eigen::VectorXd, double, double, double, double);


// [[Rcpp::export]]
Rcpp::List tdlmmBinomial(const Rcpp::List model)
{

  int i, j, k, t;
  // ---- Set up general control variables ----
  tdlmCtr *ctr = new tdlmCtr;
  ctr->iter = as<int>(model["nIter"]);
  ctr->burn = as<int>(model["nBurn"]);
  ctr->thin = as<int>(model["nThin"]);
  // ctr->threads = model["threads"];
  // omp_set_num_threads(ctr->threads);
  // Eigen::setNbThreads(ctr->threads);
  ctr->nRec = floor(ctr->iter / ctr->thin);
  ctr->nTrees = as<int>(model["nTrees"]);
  ctr->Y = as<Eigen::VectorXd>(model["Y"]);
  ctr->n = (ctr->Y).size();
  ctr->binomialSize = as<Eigen::VectorXd>(model["binomialSize"]);
  ctr->kappa = ctr->Y - 0.5 * (ctr->binomialSize);
  ctr->Lambda = ctr->kappa;
  (ctr->Omega).resize(ctr->n); (ctr->Omega).setOnes();

  ctr->Z = as<Eigen::MatrixXd>(model["Z"]);
  ctr->Zw = ctr->Z;
  ctr->pZ = (ctr->Z).cols();
  Eigen::MatrixXd VgInv = (ctr->Z).transpose() * (ctr->Z);
  VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg = VgInv.inverse();
  VgInv.resize(0,0);

  // ---- Pre-calculate single node tree matrices ----
  std::vector<exposureDat*> Exp;
  Rcpp::List exp_dat = as<Rcpp::List>(model["X"]);
  ctr->nExp = exp_dat.size();
  for (i = 0; i < ctr->nExp; ++i) {
    Exp.push_back(new exposureDat(
        as<Eigen::MatrixXd>(as<Rcpp::List>(exp_dat[i])["Tcalc"])));
  }
  ctr->pX = Exp[0]->pX;
  ctr->nSplits = 0;
  ctr->interaction = as<int>(model["interaction"]);
  ctr->nMix = 0;
  if (ctr->interaction) {
    ctr->nMix += int (ctr->nExp * (ctr->nExp - 1.0) / 2.0);
    if (ctr->interaction == 2) {
      ctr->nMix += ctr->nExp;
    }
  }
  // ctr->modZeta = 1;
  ctr->modKappa = 1;//as<double>(model["mix.prior"]);
  //std::max(ctr->nExp, ctr->nMix) / double(ctr->nTrees);


  // ---- Logs ----;
  tdlmLog *dgn = new tdlmLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec); (dgn->gamma).setZero();
  (dgn->kappa).resize(ctr->nRec); (dgn->kappa).setZero();
  (dgn->nu).resize(ctr->nRec); (dgn->nu).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec); (dgn->tau).setZero();
  (dgn->muExp).resize(ctr->nExp, ctr->nRec); (dgn->muExp).setZero();
  if (ctr->interaction > 0) {
    (dgn->muMix).resize(ctr->nMix, ctr->nRec); (dgn->muMix).setZero();
    (dgn->mixInf).resize(ctr->nMix, ctr->nRec); (dgn->mixInf).setZero();
    (dgn->mixCount).resize(ctr->nMix, ctr->nRec); (dgn->muMix).setZero();
  } else {
    (dgn->muMix).resize(1, 1); (dgn->muMix).setZero();
    (dgn->mixInf).resize(1, 1); (dgn->mixInf).setZero();
    (dgn->mixCount).resize(1, 1); (dgn->muMix).setZero();
  }
  (dgn->expProb).resize(ctr->nExp, ctr->nRec); (dgn->expProb).setZero();
  (dgn->expInf).resize(ctr->nExp, ctr->nRec); (dgn->expInf).setZero();
  (dgn->expCount).resize(ctr->nExp, ctr->nRec); (dgn->expCount).setZero();
  (dgn->fhat).resize(ctr->n); (dgn->fhat).setZero();
  (dgn->termNodes).resize(ctr->nTrees, ctr->nRec); (dgn->termNodes).setZero();
  (dgn->termNodes2).resize(ctr->nTrees, ctr->nRec); (dgn->termNodes2).setZero();
  (dgn->tree1Exp).resize(ctr->nTrees, ctr->nRec); (dgn->tree1Exp).setZero();
  (dgn->tree2Exp).resize(ctr->nTrees, ctr->nRec); (dgn->tree2Exp).setZero();



  // ---- Initial draws ----
  (ctr->fhat).resize(ctr->n); (ctr->fhat).setZero();
  ctr->R = ctr->Lambda;
  (ctr->gamma).resize(ctr->pZ);
  (ctr->totTermExp).resize(ctr->nExp); (ctr->totTermExp).setZero();
  (ctr->sumTermT2Exp).resize(ctr->nExp); (ctr->sumTermT2Exp).setZero();
  (ctr->muExp).resize(ctr->nExp); (ctr->muExp).setOnes();
  if (ctr->interaction) {
    (ctr->totTermMix).resize(ctr->nExp, ctr->nExp); (ctr->totTermMix).setZero();
    (ctr->sumTermT2Mix).resize(ctr->nExp, ctr->nExp); (ctr->sumTermT2Mix).setZero();
    (ctr->muMix).resize(ctr->nExp, ctr->nExp); (ctr->muMix).setOnes();
    (ctr->mixInf).resize(ctr->nExp, ctr->nExp); (ctr->mixInf).setZero();
    (ctr->mixCount).resize(ctr->nExp, ctr->nExp); (ctr->mixCount).setZero();
  }
  ctr->totTerm = 0;
  ctr->sumTermT2 = 0;
  ctr->nu = 1; // Need to define for first update of nu
  tdlmModelEstBinomial(ctr);
  double xiInv = R::rgamma(1, 0.5);
  ctr->nu = 1.0 / R::rgamma(0.5, 1.0 / xiInv);
  (ctr->tau).resize(ctr->nTrees); (ctr->tau).setOnes();
  for (t = 0; t < ctr->nTrees; t++) {
    xiInv = R::rgamma(1, 0.5);
    (ctr->tau)(t) = 1.0 / R::rgamma(0.5, 1.0 / xiInv);
  }
  ctr->nTerm.resize(ctr->nTrees); (ctr->nTerm).setOnes();
  ctr->nTerm2.resize(ctr->nTrees); (ctr->nTerm2).setOnes();
  (ctr->Rmat).resize(ctr->n, ctr->nTrees); (ctr->Rmat).setZero();





  // ---- Create trees ----
  (ctr->tree1Exp).resize(ctr->nTrees); (ctr->tree1Exp).setZero();
  (ctr->tree2Exp).resize(ctr->nTrees); (ctr->tree2Exp).setZero();
  ctr->expProb = as<Eigen::VectorXd>(model["expProb"]);
  (ctr->expCount).resize((ctr->expProb).size());
  (ctr->expInf).resize((ctr->expProb).size());
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




  if (ctr->verbose)
    Rcout << "Burn-in % complete \n" <<
      "[0--------25--------50--------75--------100]\n '";
  double burnProgInc = (ctr->burn / 42.0);
  double burnProgMark = burnProgInc;
  double iterProgInc = (ctr->iter / 42.0);
  double iterProgMark = double(ctr->burn) + iterProgInc;

  // ---- MCMC ----
  std::size_t s;
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }


    // -- Update trees --
    ctr->R = ctr->Lambda - ctr->fhat + (ctr->Rmat).col(0);
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
      (ctr->totTermMix).setZero(); (ctr->sumTermT2Mix).setZero();
    }

    for (t = 0; t < ctr->nTrees; t++) {

      tdlmmBinomialTreeMCMC(t, trees1[t], trees2[t], ctr, dgn, Exp);
      ctr->fhat += (ctr->Rmat).col(t);

      if (t < ctr->nTrees - 1) {
        ctr->R += (ctr->Rmat).col(t + 1) - (ctr->Rmat).col(t);
      }
    }


    // -- Pre-calculations for control and variance --
    ctr->sumTermT2 = (ctr->sumTermT2Exp).sum();
    ctr->totTerm = (ctr->totTermExp).sum();
    if(ctr->interaction) {
      ctr->sumTermT2 += (ctr->sumTermT2Mix).sum();
      ctr->totTerm += (ctr->totTermMix).sum();
    }



    // -- Update control --
    ctr->R = ctr->Lambda - ctr->fhat;
    tdlmModelEstBinomial(ctr);


    xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->nu)));
    ctr->nu = 1.0 / R::rgamma(0.5 * ctr->totTerm + 0.5,
                              1.0 / (0.5 * ctr->sumTermT2 + xiInv));

    for (i = 0; i < ctr->nExp; ++i) {
      xiInv = R::rgamma(1, (ctr->muExp(i)) / (ctr->muExp(i) + 1.0));
      ctr->muExp(i) = 1.0 / R::rgamma(0.5 * ctr->totTermExp(i) + 0.5,
                 1.0 / (0.5 * ctr->sumTermT2Exp(i) / (ctr->nu) + xiInv));
      if (ctr->interaction) {
        for (j = i; j < ctr->nExp; ++j) {
          if ((j > i) || (ctr->interaction == 2)) {
            xiInv = R::rgamma(1, (ctr->muMix(j, i)) / (ctr->muMix(j, i) + 1.0));
            ctr->muMix(j, i) = 1.0 / R::rgamma(0.5 * ctr->totTermMix(j, i) + 0.5,
                       1.0 / (0.5 * ctr->sumTermT2Mix(j, i) / (ctr->nu) + xiInv));
            if (ctr->muMix(j, i) != ctr->muMix(j, i)) {
              stop("\nNaN values occured during model run, rerun model.\n");
            }
          }
        }
      }
    }


    // -- Update exposure selection probability
    if ((ctr->b > 1000) || (ctr->b > (0.5 * ctr->burn))) {
      // double modKappaNew = exp(log(ctr->modKappa) + R::rnorm(0, 0.5));
      double modKappaNew = R::rgamma(1.0, ctr->nTrees/4.0);
      double mhrDir =
        logDirichletDensity(ctr->expProb,
                            ((ctr->expCount).array() + modKappaNew).matrix()) -
        // R::dgamma(modKappaNew, 0.5, ctr->nTrees/2.0, true) -
        logDirichletDensity(ctr->expProb,
                            ((ctr->expCount).array() + ctr->modKappa).matrix());
        // R::dgamma(ctr->modKappa, 0.5, ctr->nTrees/2.0, true);

      if (log(R::runif(0, 1)) < mhrDir) {
        ctr->modKappa = modKappaNew;
      }

      ctr->expProb = rDirichlet(((ctr->expCount).array() + ctr->modKappa).matrix());
    }





    // -- Record --
    if (ctr->record > 0) {
      dgn->fhat += ctr->fhat;
      (dgn->gamma).col(ctr->record - 1) = ctr->gamma;
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
        k = 0;
        for (i = 0; i < ctr->nExp; ++i) {
          for (j = i; j < ctr->nExp; ++j) {
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

    // -- Progress --
    if (ctr->verbose) {
      if (ctr->b < ctr->burn) {
        if (ctr->b >= burnProgMark) {
          Rcout << "'";
          burnProgMark += burnProgInc;
        }
      } else {
        if (ctr->b >= iterProgMark) {
          Rcout << "'";
          iterProgMark += iterProgInc;
        }
        if (ctr->b == (ctr->burn + ctr->iter)) {
          Rcout << "\n";
        }
      }

      // if (ctr->b % 100 == 0)
      //   Rcout << ".";
      // if ((ctr->b > ctr->burn) && ((ctr->b - ctr->burn) % 5000 == 0))
      //   Rcout << "5000\n";
      if (ctr->b == ctr->burn) {
        timediff = difftime(time(NULL), startTime);
        // if (timediff > 3600) {
        //   Rprintf("\nBurn-in time: %.2g hours", round(100 * timediff / 3600) / 100);
        // } else if (timediff > 60) {
        //   Rprintf("\nBurn-in time: %.2g minutes", round(100 * timediff / 60) / 100);
        // } else {
        //   Rprintf("\nBurn-in time: %.2g seconds", round(100 * timediff) / 100);
        // }

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
    }
  } // -- End MCMC --




  Eigen::MatrixXd DLM((dgn->DLMexp).size(), 7);

  for (s = 0; s < (dgn->DLMexp).size(); ++s) {
    DLM.row(s) = dgn->DLMexp[s];
  }


  Eigen::VectorXd nu = dgn->nu;
  Eigen::VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  Eigen::VectorXd kappa = dgn->kappa;
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
  Eigen::MatrixXd mixCount = (dgn->mixCount).transpose();
  Eigen::MatrixXd muExp = (dgn->muExp).transpose();
  Eigen::MatrixXd muMix(1, 1); muMix.setZero();
  Eigen::MatrixXd MIX(0, 10); MIX.setZero();
  if (ctr->interaction) {
    muMix.resize((dgn->muMix).cols(), (dgn->muMix).rows());
    muMix = (dgn->muMix).transpose();
    MIX.resize((dgn->MIXexp).size(), 10);
    for (s = 0; s < (dgn->MIXexp).size(); ++s) {
      MIX.row(s) = dgn->MIXexp[s];
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
                            Named("kappa") = wrap(kappa),
                            Named("muExp") = wrap(muExp),
                            Named("muMix") = wrap(muMix)));
}






void tdlmmBinomialTreeMCMC(int t, Node *tree1, Node *tree2,
                              tdlmCtr *ctr, tdlmLog *dgn,
                              std::vector<exposureDat*> Exp)
{

  int m1, m2, newExp, success, step1, step2;
  double stepMhr, ratio;
  double m1Var, m2Var, mixVar, newExpVar, newMixVar, treeVar;
  size_t s;
  std::vector<Node*> term1, term2, newTerm;
  Node* newTree = 0;
  treeMixBinomMHR mhr0, mhr;

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
  Eigen::VectorXd ZwtR = (ctr->Zw).transpose() * (ctr->R);


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
    stepMhr = tdlmProposeTree(tree1, Exp[m1], ctr, step1);
    success = tree1->isProposed();
    newTerm = tree1->listTerminal(1);

    // switch exposures
  } else {
    newExp = sampleInt(ctr->expProb);
    if (newExp != m1) {
      success = 1;
      // stepMhr = log(ctr->expProb(m1)) - log(ctr->expProb(newExp));
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
    mhr0 = mixBinomMHR(term1, term2, ctr, ZwtR, treeVar, m1Var, m2Var, mixVar);
  } else {
    mhr0 = mixBinomMHR(term1, term2, ctr, ZwtR, treeVar, m1Var, m2Var, mixVar);
  }

  if (success) {
    mhr = mixBinomMHR(newTerm, term2, ctr, ZwtR, treeVar, newExpVar, m2Var, newMixVar);

    ratio =
      stepMhr +
      mhr.logVThetaChol - mhr0.logVThetaChol +
      0.5 * (mhr.beta - mhr0.beta) -
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
    stepMhr = tdlmProposeTree(tree2, Exp[m2], ctr, step2);
    success = tree2->isProposed();
    newTerm = tree2->listTerminal(1);

    // switch exposures
  } else {
    newExp = sampleInt(ctr->expProb);
    if (newExp != m2) {
      success = 1;
      // stepMhr = log(ctr->expProb(m2)) - log(ctr->expProb(newExp));
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
    mhr = mixBinomMHR(term1, newTerm, ctr, ZwtR, treeVar, m1Var, newExpVar, newMixVar);


    ratio =
      stepMhr +
      mhr.logVThetaChol - mhr0.logVThetaChol +
      0.5 * (mhr.beta - mhr0.beta) -
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
  (ctr->tau)(t) = 1.0 / R::rgamma(0.5 * totTerm + 0.5,
   1.0 / (0.5 * tauT2 / (ctr->nu) + xiInv));
  (ctr->nTerm)(t) = mhr0.nTerm1;
  (ctr->nTerm2)(t) = mhr0.nTerm2;
  (ctr->tree1Exp)(t) = m1;
  (ctr->tree2Exp)(t) = m2;
  (ctr->expCount)(m1)++;
  (ctr->expCount)(m2)++;
  (ctr->expInf)(m1) += (ctr->tau)(t);
  (ctr->expInf)(m2) += (ctr->tau)(t);
  (ctr->totTermExp)(m1) += mhr0.nTerm1;
  (ctr->totTermExp)(m2) += mhr0.nTerm2;
  (ctr->sumTermT2Exp)(m1) += mhr0.term1T2 / (ctr->tau)(t);
  (ctr->sumTermT2Exp)(m2) += mhr0.term2T2 / (ctr->tau)(t);
  if (mixVar != 0) {
    if (m1 <= m2) {
      (ctr->totTermMix)(m2, m1) += mhr0.nTerm1 * mhr0.nTerm2;
      (ctr->sumTermT2Mix)(m2, m1) += mhr0.mixT2 / ctr->tau[t];
      (ctr->mixInf)(m2, m1) += ((ctr->tau)(t));
      (ctr->mixCount)(m2, m1)++;
    } else {
      (ctr->totTermMix)(m1, m2) += mhr0.nTerm1 * mhr0.nTerm2;
      (ctr->sumTermT2Mix)(m1, m2) += mhr0.mixT2 / ctr->tau[t];
      (ctr->mixInf)(m1, m2) += ((ctr->tau)(t));
      (ctr->mixCount)(m1, m2)++;
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


treeMixBinomMHR mixBinomMHR(std::vector<Node*> nodes1, std::vector<Node*> nodes2,
                  tdlmCtr *ctr, Eigen::VectorXd ZwtR,
                  double treeVar, double m1Var, double m2Var, double mixVar)
{
  treeMixBinomMHR out;
  int pX1 = nodes1.size();
  int pX2 = nodes2.size();
  int pXd = pX1 + pX2;
  int interaction = 0;
  if (mixVar != 0) {
    pXd += pX1 * pX2;
    interaction = 1;
  }

  out.Xd.resize(ctr->n, pXd); out.Xd.setZero();
  Eigen::VectorXd diagVar(pXd); diagVar.setZero();

  int i, j, k;
  for (i = 0; i < pX1; ++i) {
    out.Xd.col(i) = (nodes1[i]->nodevals)->X;
    diagVar(i) = 1.0 / (m1Var * treeVar);
  }

  for (j = 0; j < pX2; ++j) {
    k = pX1 + j;
    out.Xd.col(k) = (nodes2[j]->nodevals)->X;
    diagVar(k) = 1.0 / (m2Var * treeVar);
  }

  if (interaction) {
    for (i = 0; i < pX1; ++i) {
      for (j = 0; j < pX2; ++j) {
        k = pX1 + pX2 + i * pX2 + j;
        out.Xd.col(k) =
          (((nodes1[i]->nodevals)->X).array() *
          ((nodes2[j]->nodevals)->X).array()).matrix();
        diagVar(k) = 1.0 / (mixVar * treeVar);
      }
    }
  }

  // calculate MHR
  const Eigen::MatrixXd ZwtX = ctr->Zw.transpose() * out.Xd;
  const Eigen::MatrixXd VgZwtX = ctr->Vg * ZwtX;
  const Eigen::MatrixXd Xdw = (ctr->Omega).asDiagonal() * out.Xd;
  Eigen::MatrixXd tempV(pXd, pXd);
    tempV.triangularView<Eigen::Lower>() =
      Xdw.transpose() * out.Xd - ZwtX.transpose() * VgZwtX;
    tempV.diagonal() += diagVar;

  Eigen::MatrixXd VTheta(pXd, pXd);
    VTheta.triangularView<Eigen::Lower>() =
      tempV.selfadjointView<Eigen::Lower>().llt().solve(
          Eigen::MatrixXd::Identity(pXd, pXd));

  Eigen::VectorXd XtVzInvR = Xdw.transpose() * ctr->R;
    XtVzInvR.noalias() -= VgZwtX.transpose() * ZwtR;
  const Eigen::VectorXd ThetaHat =
    VTheta.selfadjointView<Eigen::Lower>() * XtVzInvR;
  const Eigen::MatrixXd VThetaChol =
    VTheta.selfadjointView<Eigen::Lower>().llt().matrixL();

  out.drawAll = ThetaHat;
  out.drawAll.noalias() += VThetaChol * as<Eigen::VectorXd>(rnorm(pXd, 0, 1));
  out.draw1 = out.drawAll.head(pX1);
  out.term1T2 = (out.draw1).dot(out.draw1);
  out.nTerm1 = double(pX1);
  out.draw2 = out.drawAll.segment(pX1, pX2);
  out.term2T2 = (out.draw2).dot(out.draw2);
  out.nTerm2 = double(pX2);
  if (interaction) {
    out.drawMix = out.drawAll.tail(pXd - pX1 - pX2);
    out.mixT2 = (out.drawMix).dot(out.drawMix);
  }
  out.beta = ThetaHat.dot(XtVzInvR);
  out.logVThetaChol = VThetaChol.diagonal().array().log().sum();

  out.pXd = pXd;
  return(out);

}

