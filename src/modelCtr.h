#include <RcppEigen.h>
using namespace Rcpp;

struct modelCtr {
public:
  bool verbose, diagnostics;
  int n, pZ, pX, nRec, nSplits, nTrees;
  int b, iter, thin, burn, record, threads, shrinkage;
  double sigma2, xiInvSigma2, nu, VTheta1Inv, totTerm, sumTermT2;
  double modKappa, modZeta;
  std::vector<double> stepProb, treePrior;
  Eigen::VectorXd Y;
  Eigen::MatrixXd Z;
  Eigen::VectorXd R;
  Eigen::MatrixXd Rmat;
  Eigen::MatrixXd Vg;
  Eigen::MatrixXd VgInv;
  Eigen::MatrixXd VgChol;
  Eigen::VectorXd X1;
  Eigen::VectorXd ZtX1;
  Eigen::VectorXd VgZtX1;
  Eigen::VectorXd gamma;
  Eigen::VectorXd fhat;
  Eigen::VectorXd tau;
  
  // Binomial
  bool binomial;
  Eigen::VectorXd Omega;
  Eigen::MatrixXd Zw;
  Eigen::VectorXd kappa;
  Eigen::VectorXd binomialSize;
  Eigen::VectorXd Lambda;
};

struct tdlmCtr : modelCtr {
public:
  Eigen::VectorXd nTerm;

  // Mixtures
  int interaction, nExp, nMix;
  double modZeta, modKappa;
  Eigen::VectorXd expProb;
  Eigen::VectorXd expCount;
  Eigen::MatrixXd mixCount;
  Eigen::VectorXd expInf;
  Eigen::MatrixXd mixInf;
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


struct tdlmLog {
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
  Eigen::VectorXd kappa;
  Eigen::MatrixXd termNodes2;
  Eigen::MatrixXd expCount;
  Eigen::MatrixXd mixCount;
  Eigen::MatrixXd expProb;
  Eigen::MatrixXd expInf;
  Eigen::MatrixXd mixInf;
  Eigen::MatrixXd tree1Exp;
  Eigen::MatrixXd tree2Exp;
  Eigen::MatrixXd muExp;
  Eigen::MatrixXd muMix;
};

struct dlmtreeCtr : modelCtr {
public:
  int pM;
  double XcenterIdx;
  std::vector<double> stepProbMod, treePriorMod;
  Eigen::VectorXd nTerm;
  Eigen::VectorXd nTermMod;
  Eigen::MatrixXd exDLM;
  Eigen::VectorXd modCount;
  Eigen::VectorXd modInf;
  
  // Gaussian Process
  Eigen::MatrixXd X;
  Eigen::MatrixXd XtXall;
  Eigen::MatrixXd ZtXall;
  Eigen::MatrixXd VgZtXall;
  Eigen::MatrixXd VThetaInvall;
  Eigen::MatrixXd DistMat;
  Eigen::MatrixXd LambdaInv;
  Eigen::MatrixXd LambdaInvNew;
  double phi, phiNew, phiMH, phiMHNew;
  double logLambdaDet, logLambdaDetNew;
  int covarType;

  // Mixtures
  // int interaction, nExp, nMix;
  // double modZeta, modKappa;
  // Eigen::VectorXd expProb;
  // Eigen::VectorXd expCount;
  // Eigen::MatrixXd mixCount;
  // Eigen::VectorXd expInf;
  // Eigen::MatrixXd mixInf;
  // Eigen::VectorXd nTerm2;
  // Eigen::VectorXd tree1Exp;
  // Eigen::VectorXd tree2Exp;
  // Eigen::VectorXd totTermExp;
  // Eigen::MatrixXd totTermMix;
  // Eigen::VectorXd sumTermT2Exp;
  // Eigen::MatrixXd sumTermT2Mix;
  // Eigen::VectorXd muExp;
  // Eigen::MatrixXd muMix;
};


struct dlmtreeLog {
public:
  // General model logs
  Eigen::MatrixXd gamma;
  Eigen::VectorXd sigma2;
  Eigen::VectorXd nu;
  Eigen::MatrixXd tau;
  Eigen::VectorXd fhat;  
  
  // Modifier tree logs
  std::vector<Eigen::VectorXd> treeModAccept;
  Eigen::MatrixXd termNodesMod;
  Eigen::VectorXd modKappa;
  Eigen::MatrixXd modProb;
  Eigen::MatrixXd modCount;
  Eigen::MatrixXd modInf;
    
  // DLM tree logs
  std::vector<Eigen::VectorXd> treeDLMAccept;
  Eigen::MatrixXd termNodesDLM;
  
  // DLM and cumulative effect estimates
  Eigen::MatrixXd exDLM;
  Eigen::MatrixXd ex2DLM;
  Eigen::VectorXd cumDLM;
  Eigen::VectorXd cum2DLM;
  std::vector<std::string> termRule;
  std::vector<Eigen::VectorXd> DLMexp;
  
  // GP
  Eigen::VectorXd phi;

  // Mixtures
  // std::vector<Eigen::VectorXd> MIXexp;
  // Eigen::VectorXd kappa;
  // Eigen::MatrixXd termNodes2;
  // Eigen::MatrixXd expCount;
  // Eigen::MatrixXd mixCount;
  // Eigen::MatrixXd expProb;
  // Eigen::MatrixXd expInf;
  // Eigen::MatrixXd mixInf;
  // Eigen::MatrixXd tree1Exp;
  // Eigen::MatrixXd tree2Exp;
  // Eigen::MatrixXd muExp;
  // Eigen::MatrixXd muMix;
};

void tdlmModelEst(modelCtr *ctr);
// Binomial model
Eigen::VectorXd rcpp_pgdraw(Eigen::VectorXd, Eigen::VectorXd);
// void tdlmModelEstBinomial(modelCtr *ctr);

void dlmtreeRecDLM(dlmtreeCtr* ctr, dlmtreeLog* dgn);
class Node;
class exposureDat;
class modDat;
double tdlmProposeTree(Node* tree, exposureDat* Exp, modelCtr* ctr, int step);
double modProposeTree(Node* tree, modDat* Mod, dlmtreeCtr* ctr, int step);
std::string modRuleStr(Node* n, modDat* Mod);
Eigen::VectorXd countMods(Node* tree, modDat* Mod);
void drawTree(Node* tree, Node* n, double alpha, double beta);
void updateGPMats(Node* n, dlmtreeCtr* ctr);


struct treeMHR {
public:
  Eigen::VectorXd draw;
  Eigen::VectorXd draw1, draw2, drawMix, drawAll;
  Eigen::VectorXd fitted;
  Eigen::MatrixXd tempV;
  Eigen::MatrixXd Xd;
  double logVThetaChol, beta, termT2;
  double nNodes, nModTerm, nDlmTerm, totTerm, nTerm;
  double term1T2, term2T2, mixT2, nTerm1, nTerm2;
  int pXd;
};


class progressMeter {
public:
  progressMeter(modelCtr* ctr);
  ~progressMeter();
  modelCtr* ctr;
  double burnProgMark, burnProgInc, iterProgMark, iterProgInc, timediff;
  time_t startTime;
  void printMark();
};