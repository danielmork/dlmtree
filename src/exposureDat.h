#include <RcppEigen.h>
using Eigen::VectorXd;
using Eigen::MatrixXd;
class Node;

class exposureDat {
public:
  int n, nSplits, pX, pZ;
  bool preset, se, lowmem;
  exposureDat(MatrixXd Tcalc_in); // Binomial DLM
  exposureDat(MatrixXd Tcalc_in, MatrixXd Z_in,
              MatrixXd Vg_in); // Gaussian DLM
  exposureDat(MatrixXd X_in, MatrixXd SE_in, VectorXd Xsplits_in,
              MatrixXd Xcalc_in, MatrixXd Tcalc_in, 
              bool lowmem_in = 0); // Binomial DLNM
  exposureDat(MatrixXd X_in, MatrixXd SE_in, VectorXd Xsplits_in,
              MatrixXd Xcalc_in, MatrixXd Tcalc_in, MatrixXd Z_in,
              MatrixXd Vg_in, bool lowmem_in = 0); // Gaussian DLNM
  ~exposureDat();

  MatrixXd X;
  MatrixXd Z;
  MatrixXd Vg;
  MatrixXd SE;
  VectorXd Xsplits;

  MatrixXd Xcalc;
  MatrixXd ZtXcalc;
  MatrixXd VgZtXcalc;

  MatrixXd Tcalc;
  MatrixXd ZtTcalc;
  MatrixXd VgZtTcalc;
  
  std::vector<MatrixXd> Xsave;
  std::vector<MatrixXd> ZtXsave;
  std::vector<MatrixXd> VgZtXsave;

  void updateNodeVals(Node*);
};
