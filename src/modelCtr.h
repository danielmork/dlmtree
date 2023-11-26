#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Eigen;

/**
 * @brief Data container for model control variables. Passed as pointer throughout model functions.
 * 
 */
struct modelCtr {
public:
  bool verbose, diagnostics, debug;
  int n, pZ, pZ1, pX, nRec, nSplits, nTrees;
  int b, iter, thin, burn, record, threads, shrinkage;
  double sigma2, xiInvSigma2, nu, VTheta1Inv, totTerm, sumTermT2;
  double modKappa, modZeta;
  std::vector<double> stepProb, treePrior, treePrior2;
  VectorXd Y0;         // Fixed response
  MatrixXd Z;          // Design matrix for fixed effect
  VectorXd R;          // Partial residual (Also, Y - fhat): Store the current one -> update the next one
  MatrixXd Rmat;       // Each column is partial residual
  MatrixXd Vg;         // V_gamma
  MatrixXd VgInv;      // V_gamma inverse
  MatrixXd VgChol;     // V_gamma cholesky decomposition
  VectorXd X1;        
  VectorXd ZtX1;       // Z transponse * X1
  VectorXd VgZtX1;     // V_gamma * Z transpose * X1 -> Maybe for tdlnm?
  VectorXd gamma;
  VectorXd fhat;
  VectorXd tau;

  // Monotone
  VectorXd zirtGamma0;      // confounding coefficients
  VectorXd zirtGamma;
  MatrixXd zirtSigma;
  bool zirtUpdateSigma;
  VectorXd zirtSplitCounts;
  VectorXd timeSplitProb0;       // fitted values
  VectorXd timeSplitProbs;
  VectorXd timeSplitCounts;
  double timeKappa;
  bool updateTimeKappa;        // tau for IG
  
  // Binomial ----------------------------------------------
  bool binomial;
  VectorXd Omega;        // The latent variable for Polya-Gamma
  MatrixXd Zw;           // Z * Omega
  MatrixXd Zw1;          // Z * Omega1 (binom for Zinb)
  VectorXd kappa;        // Kappa = y_i - n_i/2
  VectorXd Ystar;        // Ystar which is updated every iteration of MCMC, also called z1 in ZINB
  VectorXd binomialSize; // n from Binomial (n, p)
  VectorXd Lambda;       

  // ZINB & NB --------------------------------------------------
  bool zinb; // Indicator boolean for ZINB

  // Binary component of ZINB (labelled with 1)
  MatrixXd Z1;         // Design matrix for fixed effect
  VectorXd b1;         // coefficients
  MatrixXd Vg1;        // V_gamma
  MatrixXd VgInv1;     // V_gamma inverse
  MatrixXd VgChol1;    // V_gamma cholesky decomposition
  VectorXd omega1;     // Polya-gamma latent variable (nx1)
  VectorXd z1;
  MatrixXd Sigma1;     // Var-Cov for MCMC update (nxn)
  VectorXd Mu1;        // Mean for MCMC update (nx1)
  
  // Count component(Negative Binomial) of ZINB (labelled with 2)
  int nStar;                  // Number of At-risk observations (w = 1)
  int yZeroN;                 // Number of zeros in the data
  VectorXd b2;         // coefficients
  std::vector<int> yZeroIdx;  // A vector of indices where y = 0
  VectorXd omega2;     // Polya-gamma latent variable (n x 1)
  MatrixXd Sigma2;     // Var-Cov for MCMC update (nStar x nStar)
  VectorXd Mu2;        // Mean for MCMC update (nStar x 1)
  MatrixXd Zstar;      // Z with Z with non At-risk individuals zeroed out
  VectorXd z2;         // Ystar

  // Dispersion parameter component of ZINB
  int r;                      // dispersion parameter of negative binomial
  VectorXd rVec;       // a vector of dispersion parameter: rep(r, n)
  double MHratio;             // Metropolis-Hasting ratio

  // Updating at-risk component
  VectorXd w;          // At-risk latent variable
  std::vector<int> NBidx;     // Vector containing non-zero y indices
  
  VectorXd Ytemp;      // Fixed response
  MatrixXd Ztemp;      // Design matrix for fixed effect
  VectorXd Rtemp;      // Partial residual (Also, Y - fhat): Store the current one -> update the next one
  MatrixXd Rmat_temp;  // Each column is partial residual
  
  VectorXd ones;       // Vector of ones
};

struct tdlmCtr : modelCtr { // tdlmCtr: Child class of modelCtr
public:
  VectorXd nTerm;

  // Mixtures
  int interaction, nExp, nMix;
  double modZeta, modKappa;
  VectorXd expProb;        // Probability to choose exposure
  VectorXd expCount;       // How many trees use a certain exposure: length of exposures we have
  MatrixXd mixCount;       // Same thing with tree pairs (triangle)
  VectorXd expInf;
  MatrixXd mixInf;
  VectorXd nTerm2;
  VectorXd tree1Exp;       // Exposure of tree1 of #A tree pairs
  VectorXd tree2Exp;       // Exposure of tree2 of #A tree pairs
  VectorXd totTermExp;     // tree terminal nodes related to exposure
  MatrixXd totTermMix;     // tree1 with 3 x tree2 with 4 -> (1, 2) = 12: row, column represents exposures (triangle)
  VectorXd sumTermT2Exp;   // Terminal node effect squared
  MatrixXd sumTermT2Mix;   // Same with mixture
  VectorXd muExp;          // Exposure specific-variance parameter
  MatrixXd muMix;          // Same with mixture
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

  // Monotone
  MatrixXd timeProbs;
  MatrixXd zirtSplitCounts;
  MatrixXd zirtGamma;

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

  // ZINB
  MatrixXd b1;
  MatrixXd b2;
  VectorXd r;
  MatrixXd wMat;
};

struct dlmtreeCtr : modelCtr {
public:
  int pM;
  double XcenterIdx;
  std::vector<double> stepProbMod, treePriorMod;
  VectorXd nTerm;
  VectorXd nTermMod;
  MatrixXd exDLM;
  VectorXd modCount;
  VectorXd modInf;
  
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
  // int interaction, nExp, nMix;
  // double modZeta, modKappa;
  // VectorXd expProb;
  // VectorXd expCount;
  // MatrixXd mixCount;
  // VectorXd expInf;
  // MatrixXd mixInf;
  // VectorXd nTerm2;
  // VectorXd tree1Exp;
  // VectorXd tree2Exp;
  // VectorXd totTermExp;
  // MatrixXd totTermMix;
  // VectorXd sumTermT2Exp;
  // MatrixXd sumTermT2Mix;
  // VectorXd muExp;
  // MatrixXd muMix;
};


struct dlmtreeLog {
public:
  // General model logs
  MatrixXd gamma;
  VectorXd sigma2;
  VectorXd nu;
  MatrixXd tau;
  VectorXd fhat;  
  
  // Modifier tree logs
  std::vector<VectorXd> treeModAccept;
  MatrixXd termNodesMod;
  VectorXd modKappa;
  MatrixXd modProb;
  MatrixXd modCount;
  MatrixXd modInf;
    
  // DLM tree logs
  std::vector<VectorXd> treeDLMAccept;
  MatrixXd termNodesDLM;
  
  // DLM and cumulative effect estimates
  MatrixXd exDLM;
  MatrixXd ex2DLM;
  VectorXd cumDLM;
  VectorXd cum2DLM;
  std::vector<std::string> termRule;
  std::vector<VectorXd> DLMexp;
  
  // GP
  VectorXd phi;

  // ZINB
  MatrixXd b1; // Binary coefficient
  MatrixXd b2; // Count coefficient
  VectorXd r; // dispersion parameter
  MatrixXd wMat;

  // Mixtures
  // std::vector<VectorXd> MIXexp;
  // VectorXd kappa;
  // MatrixXd termNodes2;
  // MatrixXd expCount;
  // MatrixXd mixCount;
  // MatrixXd expProb;
  // MatrixXd expInf;
  // MatrixXd mixInf;
  // MatrixXd tree1Exp;
  // MatrixXd tree2Exp;
  // MatrixXd muExp;
  // MatrixXd muMix;
};


void dlmtreeRecDLM(dlmtreeCtr* ctr, dlmtreeLog* dgn);
class Node;
class exposureDat;
class modDat;
class NodeStruct;
void tdlmModelEst(modelCtr *ctr);
double rcpp_pgdraw(double, double); // Binomial / ZINB model
VectorXd rcpp_pgdraw(VectorXd b, VectorXd c);
double tdlmProposeTree(Node* tree, exposureDat* Exp = 0, 
                       modelCtr* ctr = 0, int step = 0,
                       double depth = 0.0);
double modProposeTree(Node* tree, modDat* Mod, dlmtreeCtr* ctr, int step);
std::string modRuleStr(Node* n, modDat* Mod);
VectorXd countMods(Node* tree, modDat* Mod);
VectorXd countTimeSplits(Node* tree, modelCtr* ctr);
void drawTree(Node* tree, Node* n, double alpha, double beta, 
              double depth = 0.0);
void drawZirt(Node* eta, tdlmCtr* ctr, NodeStruct* nsX);
double zeroInflatedTreeMHR(VectorXd timeProbs, std::vector<Node*> trees,
                           int t, double newProb);
void updateGPMats(Node* n, dlmtreeCtr* ctr);
// void dlmtreeRecDLM(dlmtreeCtr* ctr, dlmtreeLog* dgn);
void updateTimeSplitProbs(std::vector<Node*> trees, modelCtr* ctr);
int updateZirtSigma(std::vector<Node*> trees, modelCtr* ctr, 
 int curCov, std::vector<MatrixXd> zirtSigmaInv, 
 std::vector<double> zirtSigmaDet);
 void updateZirtGamma(std::vector<Node*> trees, modelCtr* ctr);


struct treeMHR {
public:
  VectorXd draw;
  VectorXd draw1, draw2, drawMix, drawAll;
  VectorXd fitted;
  MatrixXd tempV;
  MatrixXd Xd, Dtrans;
  double logVThetaChol, beta, termT2, cdf;
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


VectorXd rtmvnorm(VectorXd mu, MatrixXd sigma, int iter = 3);
double zeroToInfNormCDF(VectorXd mu, MatrixXd sigma);