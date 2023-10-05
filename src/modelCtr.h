#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * @brief Data container for model control variables. Passed as pointer throughout model functions.
 * 
 */
struct modelCtr {
public:
  bool verbose, diagnostics, debug;
  int n, pZ, pX, nRec, nSplits, nTrees;
  int b, iter, thin, burn, record, threads, shrinkage;
  double sigma2, xiInvSigma2, nu, VTheta1Inv, totTerm, sumTermT2;
  double modKappa, modZeta;
  std::vector<double> stepProb, treePrior, treePrior2;
  VectorXd Y;
  MatrixXd Z;
  MatrixXd Xd;
  VectorXd R;
  MatrixXd Rmat;
  MatrixXd Vg;
  MatrixXd VgInv;
  MatrixXd VgChol;
  VectorXd X1;
  VectorXd ZtX1;
  VectorXd VgZtX1;
  VectorXd gamma;
  VectorXd fhat;
  VectorXd tau;
  
  // Binomial
  bool binomial;
  VectorXd Omega;
  MatrixXd Zw;
  VectorXd kappa;
  VectorXd binomialSize;
  VectorXd Lambda;
};

struct tdlmCtr : modelCtr {
public:
  VectorXd nTerm;

  // Mixtures
  int interaction, nExp, nMix;
  double modZeta, modKappa;
  VectorXd expProb;
  VectorXd expCount;
  MatrixXd mixCount;
  VectorXd expInf;
  MatrixXd mixInf;
  VectorXd nTerm2;
  VectorXd tree1Exp;
  VectorXd tree2Exp;
  VectorXd totTermExp;
  MatrixXd totTermMix;
  VectorXd sumTermT2Exp;
  MatrixXd sumTermT2Mix;
  VectorXd muExp;
  MatrixXd muMix;
};


struct tdlmLog {
public:
  std::vector<VectorXd> DLMexp;
  std::vector<VectorXd> TreeAccept;
  MatrixXd gamma;
  VectorXd sigma2;
  VectorXd nu;
  MatrixXd tau;
  VectorXd fhat;
  VectorXd fhat2;
  MatrixXd termNodes;

  // Mixtures
  std::vector<VectorXd> MIXexp;
  VectorXd kappa;
  MatrixXd termNodes2;
  MatrixXd expCount;
  MatrixXd mixCount;
  MatrixXd expProb;
  MatrixXd expInf;
  MatrixXd mixInf;
  MatrixXd tree1Exp;
  MatrixXd tree2Exp;
  MatrixXd muExp;
  MatrixXd muMix;
};

struct dlmtreeCtr : modelCtr {
public:
  int pM;
  double XcenterIdx;
  std::vector<double> stepProbMod, treePriorMod;
  VectorXd nTerm;     // HDLM
  VectorXd nTermDLM;     // TDLM
  VectorXd nTermDLM1;    // TDLMM
  VectorXd nTermDLM2;    // TDLMM
  VectorXd nTermMod;
  MatrixXd exDLM;
  VectorXd modCount;
  VectorXd modInf;
  VectorXd kappa;
  
  // Gaussian Process
  MatrixXd X;
  MatrixXd XtXall;
  MatrixXd ZtXall;
  MatrixXd VgZtXall;
  MatrixXd VThetaInvall;
  MatrixXd DistMat;
  MatrixXd LambdaInv;
  MatrixXd LambdaInvNew;
  double phi, phiNew, phiMH, phiMHNew;
  double logLambdaDet, logLambdaDetNew;
  int covarType;

  // Mixtures
  int interaction, nExp, nMix;
  double modZeta, modKappa, mixKappa;

  VectorXd expProb;
  VectorXd expCount;
  MatrixXd mixCount;
  VectorXd expInf;
  MatrixXd mixInf;
  VectorXd dlmTree1Exp;
  VectorXd dlmTree2Exp;
  VectorXd totTermExp;
  MatrixXd totTermMix; 
  VectorXd sumTermT2Exp;
  MatrixXd sumTermT2Mix;
  VectorXd muExp;
  MatrixXd muMix;
};


struct dlmtreeLog {
public:
  // General model logs
  MatrixXd gamma;
  VectorXd sigma2;
  VectorXd nu;
  MatrixXd tau;
  VectorXd fhat;
  VectorXd totTerm;
  
  // Modifier tree logs
  std::vector<VectorXd> treeModAccept;
  MatrixXd termNodesMod;
  VectorXd modKappa;
  MatrixXd modProb;
  MatrixXd modCount;
  MatrixXd modInf;
    
  // DLM tree logs
  std::vector<VectorXd> treeDLMAccept;
  MatrixXd termNodesDLM;    // TDLM
  MatrixXd termNodesDLM1;   // TDLMM
  MatrixXd termNodesDLM2;   // TDLMM
  VectorXd kappa;
  VectorXd mixKappa;

  // Mixtures
  std::vector<VectorXd> MIXexp;
  MatrixXd termNodes1;
  MatrixXd termNodes2;
  MatrixXd expCount;
  MatrixXd mixCount;
  MatrixXd expProb;
  MatrixXd expInf;
  MatrixXd mixInf;
  MatrixXd dlmTree1Exp;
  MatrixXd dlmTree2Exp;
  MatrixXd muExp;
  MatrixXd muMix;
  
  // DLM and cumulative effect estimates
  MatrixXd exDLM;
  MatrixXd ex2DLM;
  VectorXd cumDLM;
  VectorXd cum2DLM;
  std::vector<VectorXd> DLMexp;
  std::vector<std::string> termRule;
  std::vector<std::string> termRuleMIX;
  
  // GP
  VectorXd phi;
};



class Node;
class exposureDat;
class modDat;
void tdlmModelEst(modelCtr *ctr);
VectorXd rcpp_pgdraw(VectorXd b, VectorXd c);
double tdlmProposeTree(Node* tree, exposureDat* Exp, modelCtr* ctr, int step);
double modProposeTree(Node* tree, modDat* Mod, dlmtreeCtr* ctr, int step);
std::string modRuleStr(Node* n, modDat* Mod);
VectorXd countMods(Node* tree, modDat* Mod);
void drawTree(Node* tree, Node* n, double alpha, double beta);
void updateGPMats(Node* n, dlmtreeCtr* ctr);
// void dlmtreeRecDLM(dlmtreeCtr* ctr, dlmtreeLog* dlmtreeLog);


struct treeMHR {
public:
  VectorXd draw;
  VectorXd draw1, draw2, drawMix, drawAll;
  VectorXd fitted;
  MatrixXd tempV;
  MatrixXd Xd, Dtrans;
  double logVThetaChol, beta, termT2, cdf;
  double nNodes, nModTerm, nDlmTerm, nDlmTerm1, nDlmTerm2, totTerm, nTerm;
  double term1T2, term2T2, mixT2, nTerm1, nTerm2;
  double pXd;
  double m1Var, m2Var;
  double nTermMix;
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


VectorXd rtmvnorm(VectorXd mu, MatrixXd sigma);
double zeroToInfNormCDF(VectorXd mu, MatrixXd sigma);