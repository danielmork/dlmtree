#include "RcppEigen.h"
#include "tdlnmCtr.h"
#include "exposureDat.h"
#include "Fncs.h"
using namespace Rcpp;

void tdlnmModelEst(tdlnmCtr *ctr)
{
  const Eigen::VectorXd ZR = (ctr->Z).transpose() * (ctr->R);
  const Eigen::VectorXd gammaHat = ctr->Vg * ZR;
  ctr->xiInvSigma2 = R::rgamma(1, 1.0 / (1.0 + 1.0 / ctr->sigma2));
  ctr->sigma2 = 1.0 / R::rgamma(0.5 * ((double) ctr->n + ctr->totTerm) + 0.5,
                              1.0 / (0.5 *
                                     ((ctr->R).dot(ctr->R) -
                                      ZR.dot(gammaHat) +
                                      ctr->sumTermT2 / ctr->nu) +
                              ctr->xiInvSigma2));
  ctr->gamma = gammaHat;
  (ctr->gamma).noalias() += ctr->VgChol *
    as<Eigen::VectorXd>(rnorm(ctr->pZ, 0, sqrt(ctr->sigma2)));

  if ((ctr->sigma2 != ctr->sigma2)) {
    Rcout << "\nSigma2:" << ctr->sigma2 << "\nNu:" <<
        ctr->nu << "\nRtR:" << (ctr->R).dot(ctr->R) << "\nRtZGZtR:" <<
        ZR.dot(gammaHat) << "\nnTerm:" <<
        ctr->totTerm << "\nsumT2:" << ctr->sumTermT2 << "\nxiInvSigma2:" <<
        ctr->xiInvSigma2 << "\nTau:" << ctr->tau;
    stop("\nNaN values occured during model run, rerun model.\n");
  }
}


void tdlnmModelEstBinomial(tdlnmCtr *ctr)
{
  int i;
  ctr->gamma = ctr->Vg.transpose() * (ctr->Zw).transpose() * ctr->R;
  (ctr->gamma).noalias() += ctr->VgChol * as<Eigen::VectorXd>(rnorm(ctr->pZ, 0, 1));

  Eigen::VectorXd psi = ctr->fhat;
  for (i = 0; i < ctr->pZ; ++i)
    psi.noalias() += (ctr->Z).col(i) * (ctr->gamma)(i);

  ctr->Omega = rcpp_pgdraw(ctr->binomialSize, psi);
  ctr->Zw = ctr->Z;
  Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ);
  for (i = 0; i < ctr->pZ; ++i) {
    (ctr->Zw).col(i).array() *= (ctr->Omega).array();
  }
  VgInv.triangularView<Eigen::Lower>() = (ctr->Z).transpose() * ctr->Zw;
  VgInv.diagonal().array() += 1/10000;
  ctr->Vg.triangularView<Eigen::Lower>() =
    VgInv.selfadjointView<Eigen::Lower>().llt().solve(
        Eigen::MatrixXd::Identity(ctr->pZ, ctr->pZ));
  ctr->VgChol = (ctr->Vg).selfadjointView<Eigen::Lower>().llt().matrixL();
  ctr->Lambda = (ctr->kappa).array() / (ctr->Omega).array();
}



double tdlnmProposeTree(Node* tree, exposureDat* Exp, tdlnmCtr* ctr, int step)
{
  int no = 0;
  double stepMhr = 0;
  std::vector<Node*> dlnmTerm, tempNodes;


  // List current tree terminal nodes
  dlnmTerm = tree->listTerminal();


  // Grow
  if (step == 0) {
    // Rcout << "G";
    no = (std::size_t) floor(R::runif(0, dlnmTerm.size())); // select node to grow

    if (dlnmTerm[no]->grow()) { // propose new split
      double nGen2 = double(tree->nGen2());
      if (dlnmTerm[no]->depth == 0) { // depth == 0
        ++nGen2;
      } else {
        if (!(dlnmTerm[no]->parent->isGen2())) {
          ++nGen2;
        }
      }
      //   if (!(dlnmTerm[no]->parent)->isGen2())
      //     ++nGen2;
      // }
      stepMhr = log((double)tree->nTerminal()) - log(nGen2) +
        2 * logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth + 1, 1) +
        logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth, 0) -
        logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth, 1);

      Exp->updateNodeVals((dlnmTerm[no]->proposed)->c1); // update node values
      // newDlnmTerm = tree->listTerminal(1); // list proposed terminal nodes
    }


    // Prune
  } else if (step == 1) {
    // Rcout << "P";
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
    // Rcout << "C";
    tempNodes = tree->listInternal();
    no = floor(R::runif(0, tempNodes.size())); // select internal nodes to change

    if (tempNodes[no]->change()) { // propose new split

      // newDlnmTerm = tree->listTerminal(1); // update proposed terminal nodes
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
}
