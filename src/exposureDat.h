#include <RcppEigen.h>
class Node;

class exposureDat {
public:
  int n, nSplits, pX, pZ;
  bool preset, se;
  exposureDat(Eigen::MatrixXd); // Binomial DLM
  exposureDat(Eigen::MatrixXd, Eigen::MatrixXd,
              Eigen::MatrixXd); // Gaussian DLM
  exposureDat(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd,
              Eigen::MatrixXd, Eigen::MatrixXd); // Binomial DLNM
  exposureDat(Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd,
              Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd,
              Eigen::MatrixXd); // Gaussian DLNM
  ~exposureDat();

  Eigen::MatrixXd X;
  Eigen::MatrixXd Z;
  Eigen::MatrixXd Vg;
  Eigen::MatrixXd SE;
  Eigen::VectorXd Xsplits;

  Eigen::MatrixXd Xcalc;
  Eigen::MatrixXd ZtXcalc;
  Eigen::MatrixXd VgZtXcalc;

  Eigen::MatrixXd Tcalc;
  Eigen::MatrixXd ZtTcalc;
  Eigen::MatrixXd VgZtTcalc;

  void updateNodeVals(Node*);
};
