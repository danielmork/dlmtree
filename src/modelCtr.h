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
  Eigen::VectorXd Y0; // Fixed response
  Eigen::MatrixXd Z; // Design matrix for fixed effect
  Eigen::VectorXd R; // Partial residual (Also, Y - fhat): Store the current one -> update the next one
  Eigen::MatrixXd Rmat; // Each column is partial residual -> 
  Eigen::MatrixXd Vg; // V_gamma
  Eigen::MatrixXd VgInv; // V_gamma inverse
  Eigen::MatrixXd VgChol; // V_gamma cholesky decomposition
  Eigen::VectorXd X1; // X1? -> Maybe for tdlnm?
  Eigen::VectorXd ZtX1; // Z transponse * X1
  Eigen::VectorXd VgZtX1; // V_gamma * Z transpose * X1 -> Maybe for tdlnm?
  Eigen::VectorXd gamma; // confounding coefficients
  Eigen::VectorXd fhat; // fitted values
  Eigen::VectorXd tau; // tau for IG
  
  // Binomial ----------------------------------------------
  bool binomial;
  Eigen::VectorXd Omega; // The latent variable for Polya-Gamma
  Eigen::MatrixXd Zw; // Z * Omega
  Eigen::VectorXd kappa; // Kappa = y_i - n_i/2
  Eigen::VectorXd Ystar; // Ystar which is updated every iteration of MCMC, also called z1 in ZINB
  Eigen::VectorXd binomialSize; // n from Binomial (n, p)
  Eigen::VectorXd Lambda; // Kappa / Omega

  // ZINB & NB --------------------------------------------------
  bool zinb; // Indicator boolean for ZINB

  // Binary component of ZINB (labelled with 1)
  Eigen::VectorXd b1;       // coefficients
  Eigen::MatrixXd Vg1;       // V_gamma
  Eigen::MatrixXd VgInv1;    // V_gamma inverse
  Eigen::MatrixXd VgChol1;   // V_gamma cholesky decomposition
  Eigen::VectorXd omega1;   // Polya-gamma latent variable (nx1)
  Eigen::MatrixXd Omega1;   // diagonal matrix with omega1 (nxn)
  Eigen::VectorXd z1;
  Eigen::MatrixXd Sigma1;   // Var-Cov for MCMC update (nxn)
  Eigen::VectorXd Mu1;      // Mean for MCMC update (nx1)

  // Spatial component of ZINB
  Eigen::MatrixXd Vp;       // V_gamma
  Eigen::MatrixXd VpInv;    // V_spPhi inverse
  Eigen::MatrixXd VpChol;   // V_spPhi cholesky decomposition
  Eigen::VectorXd omegaPhi;
  Eigen::MatrixXd OmegaPhi;
  Eigen::VectorXd zPhi;

  // Count component(Negative Binomial) of ZINB (labelled with 2)
  int nStar;                  // Number of At-risk observations (w = 1)
  int yZeroN;                 // Number of zeros in the data
  Eigen::VectorXd b2;         // coefficients
  std::vector<int> yZeroIdx;  // A vector of indices where y = 0
  Eigen::VectorXd omega2;     // Polya-gamma latent variable (n x 1)
  Eigen::VectorXd omega2star; // Polya-gamma latent variable (n* x 1)
  Eigen::MatrixXd Omega2;     // diagonal matrix with omega2 (n* x n*)
  Eigen::MatrixXd Sigma2;     // Var-Cov for MCMC update (nStar x nStar)
  Eigen::VectorXd Mu2;        // Mean for MCMC update (nStar x 1)
  Eigen::MatrixXd Zstar;      // Z with Z with non At-risk individuals zeroed out
  Eigen::VectorXd z2;         // Ystar

  // Dispersion parameter component of ZINB
  int r;                      // dispersion parameter of negative binomial
  Eigen::VectorXd rVec;       // a vector of dispersion parameter: rep(r, n)
  double MHratio;             // Metropolis-Hasting ratio

  // Spatial component
  // rho update
  bool spatial; 
  Eigen::VectorXd spNodes1;
  Eigen::VectorXd spNodes2;
  // phi update
  int spN;
  Eigen::MatrixXd areaD;
  Eigen::MatrixXd areaW;
  Eigen::MatrixXd areaQ;
  Eigen::MatrixXd areaQinv;
  Eigen::MatrixXd areaA;
  Eigen::VectorXd spPhi;
  double rho;                 // Spatial correlation
  double spTau;               // Spatial precision

  // Updating at-risk component
  Eigen::VectorXd w;            // At-risk latent variable
  std::vector<int> atRiskIdx;   // Vector containing non-zero y indices

  // Useful vector
  Eigen::VectorXd ones;
};

struct tdlmCtr : modelCtr { // tdlmCtr: Child class of modelCtr
public:
  Eigen::VectorXd nTerm;

  // Mixtures
  int interaction, nExp, nMix;
  double modZeta, modKappa;
  Eigen::VectorXd expProb; // Probability to choose exposure
  Eigen::VectorXd expCount; // How many trees use a certain exposure: length of exposures we have
  Eigen::MatrixXd mixCount; // Same thing with tree pairs (triangle)
  Eigen::VectorXd expInf;
  Eigen::MatrixXd mixInf;
  Eigen::VectorXd nTerm2;
  Eigen::VectorXd tree1Exp; // Exposure of tree1 of #A tree pairs
  Eigen::VectorXd tree2Exp; // Exposure of tree2 of #A tree pairs
  Eigen::VectorXd totTermExp; // tree terminal nodes related to exposure
  Eigen::MatrixXd totTermMix; // tree1 with 3 x tree2 with 4 -> (1, 2) = 12: row, column represents exposures (triangle)
  Eigen::VectorXd sumTermT2Exp; // Terminal node effect squared
  Eigen::MatrixXd sumTermT2Mix; // Same with mixture
  Eigen::VectorXd muExp; // Exposure specific-variance parameter
  Eigen::MatrixXd muMix; // Same with mixture
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

  // ZINB
  Eigen::MatrixXd b1;
  Eigen::MatrixXd b2;
  Eigen::VectorXd r;
  Eigen::MatrixXd wMat;
  Eigen::MatrixXd spPhi;
  Eigen::VectorXd rho;
  Eigen::VectorXd spTau;
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

  // ZINB
  Eigen::MatrixXd b1; // Binary coefficient
  Eigen::MatrixXd b2; // Count coefficient
  Eigen::VectorXd r; // dispersion parameter
  Eigen::MatrixXd wMat;
  Eigen::MatrixXd spPhi;
  Eigen::VectorXd rho;
  Eigen::VectorXd spTau;

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

// ----- Binomial Model -----
void tdlmModelEst(modelCtr *ctr);
// Binomial model
Eigen::VectorXd rcpp_pgdraw(Eigen::VectorXd, Eigen::VectorXd); // Multiple Draw from Polya-Gamma (Vector)
double rcpp_pgdraw(double, double); // Single Draw from Polya-Gamma
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