#include <RcppEigen.h>
using namespace Rcpp;

struct tdlnmCtr {
public:
  bool verbose, diagnostics;
  int n, pZ, pX, nRec, nSplits, nTrees;
  int b, iter, thin, burn, record, threads;
  double sigma2, xiInvSigma2, nu, VTheta1Inv, totTerm, sumTermT2;
  std::vector<double> stepProb, treePrior;
  Eigen::VectorXd Y;
  Eigen::VectorXd fhat;
  Eigen::VectorXd X1;
  Eigen::VectorXd ZtX1;
  Eigen::VectorXd VgZtX1;
  Eigen::VectorXd tau;
  Eigen::VectorXd nTerm;
  Eigen::VectorXd R;
  Eigen::VectorXd gamma;
  Eigen::MatrixXd Rmat;
  Eigen::MatrixXd Z;
  Eigen::MatrixXd Vg;
  Eigen::MatrixXd VgChol;

  // Binomial
  Eigen::VectorXd kappa;
  Eigen::VectorXd binomialSize;
  Eigen::VectorXd Lambda;
  Eigen::VectorXd Omega;
  Eigen::MatrixXd Zw;

  // Mixtures
  int interaction, nExp, nMix;
  double modZeta, modKappa;
  Eigen::VectorXd expProb;
  Eigen::VectorXd expCount;
  Eigen::VectorXd nTerm2;
  Eigen::VectorXd tree1Exp;
  Eigen::VectorXd tree2Exp;
  Eigen::VectorXd totTermExp;
  Eigen::MatrixXd totTermMix;
  Eigen::VectorXd sumTermT2Exp;
  Eigen::MatrixXd sumTermT2Mix;
  Eigen::VectorXd muExp;
  Eigen::MatrixXd muMix;
};


struct tdlnmLog {
public:
  std::vector<Eigen::VectorXd> DLMexp;
  std::vector<Eigen::VectorXd> TreeAccept;
  Eigen::MatrixXd gamma;
  Eigen::VectorXd sigma2;
  Eigen::VectorXd nu;
  Eigen::MatrixXd tau;
  Eigen::VectorXd fhat;
  Eigen::MatrixXd termNodes;

  // Mixtures
  std::vector<Eigen::VectorXd> MIXexp;
  Eigen::MatrixXd termNodes2;
  Eigen::MatrixXd expProb;
  Eigen::MatrixXd tree1Exp;
  Eigen::MatrixXd tree2Exp;
  Eigen::MatrixXd muExp;
  Eigen::MatrixXd muMix;
};

void tdlnmModelEst(tdlnmCtr *ctr);

// Binomial model
Eigen::VectorXd rcpp_pgdraw(Eigen::VectorXd, Eigen::VectorXd);
void tdlnmModelEstBinomial(tdlnmCtr *ctr);
