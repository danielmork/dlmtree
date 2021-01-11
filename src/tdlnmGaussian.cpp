#include <RcppEigen.h>
#include "tdlmCtr.h"
#include "exposureDat.h"
#include "Fncs.h"
using namespace Rcpp;

void tdlnmGaussianTreeMCMC(int, Node*, tdlmCtr*, tdlmLog*, exposureDat*);
double tdlmProposeTree(Node *tree, exposureDat *Exp, tdlmCtr *ctr, int step);

struct treeMHR {
public:
  Eigen::VectorXd draw;
  Eigen::MatrixXd tempV;
  Eigen::MatrixXd Xd;
  double logVThetaChol, beta, nTerm, termT2;
};

treeMHR dlnmMHR(std::vector<Node*>, tdlmCtr*, Eigen::VectorXd, double, Node*, bool);


// [[Rcpp::export]]
Rcpp::List tdlnmGaussian(const Rcpp::List model)
{
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

  ctr->Z = as<Eigen::MatrixXd>(model["Z"]);
  ctr->pZ = (ctr->Z).cols();
  Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ);
  VgInv = (ctr->Z).transpose() * (ctr->Z);
  VgInv.diagonal().array() += 1.0 / 100000.0;
  ctr->Vg = VgInv.inverse();
  VgInv.resize(0,0);
  ctr->VgChol = (ctr->Vg).llt().matrixL();

  // ---- Pre-calculate single node tree matrices ----
  exposureDat *Exp;
  if (as<int>(model["nSplits"]) == 0) {
    Exp = new exposureDat(as<Eigen::MatrixXd>(model["Tcalc"]),
                          ctr->Z, ctr->Vg);
  } else {
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
  int t;
  std::vector<Node*> trees;
  NodeStruct *ns;
  ns = new DLNMStruct(0, ctr->nSplits + 1, 1, int (ctr->pX),
                      as<Eigen::VectorXd>(model["splitProb"]),
                      as<Eigen::VectorXd>(model["timeProb"]));
  for (t = 0; t < ctr->nTrees; t++) {
    trees.push_back(new Node(0, 1));
    trees[t]->nodestruct = ns->clone();
    Exp->updateNodeVals(trees[t]);
  }
  delete ns;



  // ---- Logs ----
  tdlmLog *dgn = new tdlmLog;
  (dgn->gamma).resize(ctr->pZ, ctr->nRec); (dgn->gamma).setZero();
  (dgn->sigma2).resize(ctr->nRec); (dgn->sigma2).setZero();
  (dgn->nu).resize(ctr->nRec); (dgn->nu).setZero();
  (dgn->tau).resize(ctr->nTrees, ctr->nRec); (dgn->tau).setZero();
  (dgn->fhat).resize(ctr->n); (dgn->fhat).setZero();
  (dgn->termNodes).resize(ctr->nTrees, ctr->nRec); (dgn->termNodes).setZero();



  // ---- Initial draws ----
  (ctr->fhat).resize(ctr->n); (ctr->fhat).setZero();
  ctr->R = ctr->Y;
  (ctr->gamma).resize(ctr->pZ);
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
  ctr->nTerm.resize(ctr->nTrees);
  (ctr->nTerm).array() = 1;
  (ctr->Rmat).resize(ctr->n, ctr->nTrees); (ctr->Rmat).setZero();



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

  std::size_t s;
  // ---- MCMC ----
  for (ctr->b = 1; ctr->b <= (ctr->iter + ctr->burn); (ctr->b)++) {
    if ((ctr->b > ctr->burn) && (((ctr->b - ctr->burn) % ctr->thin) == 0)) {
      ctr->record = floor((ctr->b - ctr->burn) / ctr->thin);
    } else {
      ctr->record = 0;
    }

    // -- Update trees --
    ctr->R += (ctr->Rmat).col(0);
    (ctr->fhat).setZero();
    ctr->totTerm = 0.0; ctr->sumTermT2 = 0.0;

    for (t = 0; t < ctr->nTrees; t++) {

      tdlnmGaussianTreeMCMC(t, trees[t], ctr, dgn, Exp);
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
      (dgn->termNodes).col(ctr->record - 1) = ctr->nTerm;
      dgn->fhat += ctr->fhat;
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
    }


  }

  Eigen::MatrixXd DLM((dgn->DLMexp).size(), 8);
  for (s = 0; s < (dgn->DLMexp).size(); ++s)
    DLM.row(s) = dgn->DLMexp[s];
  Eigen::VectorXd sigma2 = dgn->sigma2;
  Eigen::VectorXd nu = dgn->nu;
  Eigen::VectorXd fhat = (dgn->fhat).array() / ctr->nRec;
  Eigen::MatrixXd gamma = (dgn->gamma).transpose();
  Eigen::MatrixXd tau = (dgn->tau).transpose();
  Eigen::MatrixXd termNodes = (dgn->termNodes).transpose();

  Eigen::MatrixXd Accept((dgn->TreeAccept).size(), 5);
  for (s = 0; s < (dgn->TreeAccept).size(); ++s)
    Accept.row(s) = dgn->TreeAccept[s];

  delete ctr;
  delete dgn;
  delete Exp;

  for (s = 0; s < trees.size(); ++s) {
    delete trees[s];
  }

  return(Rcpp::List::create(Named("DLM") = wrap(DLM),
                            Named("fhat") = wrap(fhat),
                            Named("sigma2") = wrap(sigma2),
                            Named("nu") = wrap(nu),
                            Named("tau") = wrap(tau),
                            Named("termNodes") = wrap(termNodes),
                            Named("gamma") = wrap(gamma),
                            Named("treeAccept") = wrap(Accept)));
}




void tdlnmGaussianTreeMCMC(int t, Node *tree, tdlmCtr *ctr, tdlmLog *dgn, exposureDat *Exp)
{
  int step;
  int success = 0;
  double stepMhr = 0;
  double ratio = 0;
  double treevar = (ctr->nu) * (ctr->tau)(t);
  double RtR, RtZVgZtR, xiInv;
  std::size_t s;
  std::vector<Node*> dlnmTerm, newDlnmTerm;
  treeMHR mhr0, mhr;


  // List current tree terminal nodes
  dlnmTerm = tree->listTerminal();
  // calculate part of MHR and draw new node-specific effects
  Eigen::VectorXd ZtR = (ctr->Z).transpose() * (ctr->R);
  mhr0 = dlnmMHR(dlnmTerm, ctr, ZtR, treevar, tree, 0);

  if (dlnmTerm.size() > 1) {
    step = sampleInt(ctr->stepProb, 1);
  } else { // if single terminal node, grow is only option
    step = 0;
  }

  // propose update
  stepMhr = tdlmProposeTree(tree, Exp, ctr, step);
  success = tree->isProposed();


  if (success) {
    // calculate new tree part of MHR and draw node effects
    newDlnmTerm = tree->listTerminal(1);
    mhr = dlnmMHR(newDlnmTerm, ctr, ZtR, treevar, tree, 1);

    // combine mhr parts into log-MH ratio
    RtR = (ctr->R).dot(ctr->R);
    RtZVgZtR = ZtR.dot((ctr->Vg).selfadjointView<Eigen::Lower>() * ZtR);
    ratio =
      stepMhr +
      mhr.logVThetaChol - mhr0.logVThetaChol -
      (0.5 * (ctr->n + 1.0) *
       (log(0.5 * (RtR - RtZVgZtR - mhr.beta) + ctr->xiInvSigma2) -
        log(0.5 * (RtR - RtZVgZtR - mhr0.beta) + ctr->xiInvSigma2))) -
      (log(treevar) * 0.5 * (mhr.nTerm - mhr0.nTerm));

    if (log(R::runif(0, 1)) < ratio) {
      mhr0 = mhr;
      success = 2;
      tree->accept();
      (tree->nodevals->tempV).resize(dlnmTerm.size(), dlnmTerm.size());
      tree->nodevals->tempV = mhr.tempV;
      dlnmTerm = tree->listTerminal();

    } else {
      tree->reject();
    }

  } else {
    tree->reject();
  }


  // Update variance and residuals
  xiInv = R::rgamma(1, 1.0 / (1.0 + 1.0 / (ctr->tau)(t)));
  (ctr->tau)(t) = 1.0 / R::rgamma(0.5 * mhr0.nTerm + 0.5,
                              1.0 / ((0.5 * mhr0.termT2 /
                                     (ctr->sigma2 * ctr->nu)) + xiInv));

  if ((ctr->tau)(t) != (ctr->tau)(t)) { // stop if nan values
    Rcout << ctr->sigma2 << "\n" <<
      ctr->nu << "\n" << mhr0.nTerm << "\n" << mhr0.termT2 << "\n" <<
      xiInv << "\n" <<
      t << "\n" << ctr->tau;
    stop("\nNaN values occured during model run, rerun model.\n");
  }

  (ctr->nTerm)(t) = mhr0.nTerm;
  ctr->totTerm += mhr0.nTerm;
  ctr->sumTermT2 += mhr0.termT2 / ((ctr->tau)(t));
  ctr->Rmat.col(t) = mhr0.Xd * mhr0.draw;


  // Record
  if (ctr->record > 0) {
    Eigen::VectorXd rec(8);
    rec << ctr->record, t, (dlnmTerm[0]->nodestruct)->get(1),
    (dlnmTerm[0]->nodestruct)->get(2), (dlnmTerm[0]->nodestruct)->get(3),
    (dlnmTerm[0]->nodestruct)->get(4), mhr0.draw(0), 0;
    (dgn->DLMexp).push_back(rec);

    for(s = 1; s < dlnmTerm.size(); ++s) {
      rec[2] = (dlnmTerm[s]->nodestruct)->get(1);
      rec[3] = (dlnmTerm[s]->nodestruct)->get(2);
      rec[4] = (dlnmTerm[s]->nodestruct)->get(3);
      rec[5] = (dlnmTerm[s]->nodestruct)->get(4);
      rec[6] = mhr0.draw(s);
      (dgn->DLMexp).push_back(rec);
    }

    if (ctr->diagnostics) {
      Eigen::VectorXd acc(5);
      acc << step, success, dlnmTerm.size(), stepMhr, ratio;
      (dgn->TreeAccept).push_back(acc);
    }
  }
}


treeMHR dlnmMHR(std::vector<Node*> nodes, tdlmCtr *ctr,
                Eigen::VectorXd ZtR, double var, Node* tree, bool newTree)
{
  treeMHR out;
  std::size_t s;
  int pX = int(nodes.size());


  // single terminal node, use precalculated values
  if (pX == 1) {
    double VTheta = 1.0 / (ctr->VTheta1Inv + 1.0 / var);
    double XtVzInvR = (ctr->X1).dot(ctr->R) - (ctr->VgZtX1).dot(ZtR);
    double ThetaHat = VTheta * XtVzInvR;
    double VThetaChol = sqrt(VTheta);
    (out.Xd).resize(ctr->n, 1);
    (out.Xd).col(0) = ctr->X1;
    (out.draw).resize(1);
    (out.draw)(0) = VThetaChol * R::rnorm(0, sqrt(ctr->sigma2)) + ThetaHat;
    out.beta = ThetaHat * XtVzInvR;
    out.logVThetaChol = log(VThetaChol);


  // pX > 1
  } else {
    Eigen::MatrixXd tempV(pX, pX); tempV.setZero();
    // Eigen::VectorXd XtVzInvR(pX);
    (out.Xd).resize(ctr->n, pX);
    Eigen::MatrixXd ZtX(ctr->pZ, pX); ZtX.setZero();
    Eigen::MatrixXd VgZtX(ctr->pZ, pX); VgZtX.setZero();

    for (s = 0; s < nodes.size(); ++s) {
      // XtVzInvR(i) = ((nodes[i]->nodevals)->X).dot(ctr->R);
      // XtVzInvR(i) -= ((nodes[i]->nodevals)->VgZtX).transpose() * ZtR;
      (out.Xd).col(s) = (nodes[s]->nodevals)->X;
      ZtX.col(s) = (nodes[s]->nodevals)->ZtX;
      VgZtX.col(s) = (nodes[s]->nodevals)->VgZtX;

      // if (newTree) {
      //   tempV(i, i) = (((nodes[i]->nodevals)->X).array().square()).matrix().sum();
      //   tempV(i, i) -= ((nodes[i]->nodevals)->ZtX).dot((nodes[i]->nodevals)->VgZtX);
      //   for (int j = (i + 1); j < pX; j++) {
      //     tempV(j, i) = ((nodes[i]->nodevals)->X).dot((nodes[j]->nodevals)->X);
      //     tempV(j, i) -= ((nodes[i]->nodevals)->ZtX).dot((nodes[j]->nodevals)->VgZtX);
      //     tempV(i, j) = tempV(j, i);
      //   }
      // }

    }

    if (newTree) {
      // tempV.triangularView<Eigen::Lower>() =
      //   out.Xd.transpose() * out.Xd - ZtX.transpose() * VgZtX;
      tempV = out.Xd.transpose() * out.Xd - ZtX.transpose() * VgZtX;
      out.tempV = tempV;
    } else {
      tempV = tree->nodevals->tempV;
    }

    tempV.diagonal().array() += 1.0 / var;
    const Eigen::MatrixXd VTheta = tempV.inverse();//(pX, pX);
    // VTheta.triangularView<Eigen::Lower>() = tempV.inverse();

    const Eigen::VectorXd XtVzInvR =
      out.Xd.transpose() * ctr->R - VgZtX.transpose() * ZtR;
    const Eigen::MatrixXd VThetaChol = VTheta.llt().matrixL();
    const Eigen::VectorXd ThetaHat = VTheta * XtVzInvR;

    out.draw = ThetaHat;
    (out.draw).noalias() += VThetaChol *
      as<Eigen::VectorXd>(rnorm(pX, 0, sqrt(ctr->sigma2)));
    out.beta = ThetaHat.dot(XtVzInvR);
    out.logVThetaChol = VThetaChol.diagonal().array().log().sum();

  }

  out.termT2 = (out.draw).dot(out.draw);
  out.nTerm = double(pX);
  return(out);
}

