/**
 * @file modelEst.cpp
 * @author Daniel Mork (danielmork.github.io)
 * @brief Functions for model estimation to accompany Bayesian treed distributed lag methods
 * @version 1.0
 * @date 2021-02-24
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include "RcppEigen.h"
#include "modelCtr.h"
#include "exposureDat.h"
#include "modDat.h"
#include "Fncs.h"
#include "Node.h"
#include "NodeStruct.h"
#include <random>
#include <iostream>
using namespace Rcpp;

/**
 * @brief function to update fixed effect coef., sigma^2, polya-gamma
 * 
 * @param ctr model control data
 */
void tdlmModelEst(modelCtr *ctr)
{ 
  if(!(ctr->zinb)){ // If not ZINB, go to either Binomial or Gaussian
    // Update the fixed coefficient - Calculate the mean of the gamma full conditional 
    const Eigen::VectorXd ZR = ctr->Zw.transpose() * ctr->R; // R = Y - f(DLM effects)
    ctr->gamma = ctr->Vg * ZR; // For mean vector of gamma's normal distribution: V_gamma * Zt * (Y - f)
    
    // * Update sigma^2 and xi_sigma2
    if (!(ctr->binomial)) {
      rHalfCauchyFC(&(ctr->sigma2), (double)ctr->n + (double)ctr->totTerm, 
                    ctr->R.dot(ctr->R) - ZR.dot(ctr->gamma) + ctr->sumTermT2 / ctr->nu, &(ctr->xiInvSigma2));
      if ((ctr->sigma2 != ctr->sigma2)) // ! stop if infinte or nan variance
        stop("\nNaN values (sigma) occured during model run, rerun model.\n");
    }

    // * Draw fixed effect coefficients' variance
    ctr->gamma.noalias() += ctr->VgChol * as<Eigen::VectorXd>(rnorm(ctr->pZ, 0, sqrt(ctr->sigma2))); 
      // square-root of (V_gamma * N(0, sigma^2))
      // -> chol(Vg) * sigma2 -> Covariance matrix -> Add to the fixed coefficient
      // cholesky for inverse is more efficient than naive inverse functions
      
    // * Update polya gamma vars
    if (ctr->binomial) {
      // Rcout << "Logistic Model Estimation \n";
      // DLM + Z * Gamma (Equation 24) for PG's second parameter
      Eigen::VectorXd psi = ctr->fhat; // f part
      psi.noalias() += ctr->Z * ctr->gamma; // Z_gamma + f
      
      // Latent variable, Omega
      // Omega|psi ~ Polya-Gamma(n_i, psi_i) Appendix B.5 (Equation 26)
      // -> returns a vector of omega
      ctr->Omega = rcpp_pgdraw(ctr->binomialSize, psi); // Second line from the supplemental
      
      // asDiagonal() creates a Matrix with the input vector as diagonal elements
      // Omega.asDiagonal() is the big Omega below equation 29
      ctr->Zw = ctr->Omega.asDiagonal() * ctr->Z; // Zw = Omega * Z for binomial-specific

      // Constructing V_gamma Inverse 
      // Note V_gamma inverse and V_gamma are both symmetric
      Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ); 
      VgInv.triangularView<Eigen::Lower>() = ctr->Z.transpose() * ctr->Zw; // Lower triangle = Zt * Omega * Z 
      VgInv.diagonal().array() += 1 / 100000.0; // ZT*Omega*Z  + c*I where c = 100000
      VgInv.triangularView<Eigen::Upper>() = VgInv.transpose().eval(); // Upper triangle = transpose of the lower triangle

      // Constructing V_gamma = Inverse of V_gamma Inverse
      ctr->Vg.triangularView<Eigen::Lower>() = VgInv.inverse();
      ctr->Vg.triangularView<Eigen::Upper>() = ctr->Vg.transpose().eval();
      
      // Update the V_gamma cholesky using LLT Decomposition, Lower triangular part of matrix L
      ctr->VgChol = ctr->Vg.llt().matrixL(); 
      ctr->Ystar = (ctr->kappa).array() / ctr->Omega.array(); 
      // This is Y^* in the Dan's document and z in ZINB

      ctr->R = ctr->Ystar - ctr->fhat; // Update Ystar - f
    }

  } else { // ZINB
    // *** Step 1: Update the latent At-Risk Indicators ***
    // Step 1-1: Calculate eta1 for logit1, eta2 (binary) & logit2 & nu2 (count)
    // Rcout << "Step 1: Updating the latent indicator variable \n";
    Eigen::VectorXd eta1 = (ctr->Z) * (ctr->b1);              // Z * \gamma
    eta1.noalias() += ctr->areaA * ctr->spPhi;                // Z * \gamma + spatial random effect
    Eigen::VectorXd eta2 = (ctr->fhat);                       // dlm effect
    eta2.noalias() += ctr->Z * ctr->b2;                       // dlm effect + fixed effect

    // Rcout << "here2 \n";
    Eigen::VectorXd logit1 = 1 / (1 + exp(-(eta1).array()));
    Eigen::VectorXd logit2 = 1 / (1 + exp(-(eta2).array()));  // psi 
    Eigen::VectorXd nu2 = 1 - logit2.array();

    // Rcout << "Step 1: Resetting the latent indicator variable \n";
    ctr->w.setOnes(); // Reset w
    // Rcout << "here3 \n";
    // Sampling w
    for(int z = 0; z < (ctr->yZeroN); z++){ // For all y = 0,
      // Find the index of y = 0,
      int idx = (ctr->yZeroIdx)[z]; 

      // Bernoulli probability
      long double prob = log(logit1[idx]) + (ctr->r) * log((nu2)[idx])
                        - log(1 - logit1[idx] + logit1[idx] * pow((nu2)[idx], (ctr->r))); 

      // Update the index with a probability with either 0 or 1
      (ctr->w)[idx] = R::rbinom(1, exp(prob)); 

      // Rcout << (ctr->w)[idx];
      // Rcout << "\n";
      // Rcout << "---------- \n";

      if(isnan(ctr->w[idx])){
        (ctr->w)[idx] = 0;
      }
    }
    // Rcout << "here4 \n";
    //ctr->w.setOnes(); // Reset w (### Comment this or not to determine ZINB and NB)
    
    // Update the number of at-risk individuals
    ctr->nStar = (ctr->w).sum();

    // Update the indices of at-risk
    ctr->atRiskIdx.clear();  // Clear the vector
    for(int j = 0; j < (ctr->n); j++){
      if((ctr->w)[j] == 1){
        ctr->atRiskIdx.push_back(j); // Save the index if w = 1
      }
    }

    // Check the atRiskIdx vector
    if(ctr->atRiskIdx.size() != (ctr->nStar)){
      stop("The number of at-risk observations don't match");
    }

    // Check for NaN of the latent variable w.
    if(ctr->w != ctr->w){
      stop("NaN value occured for the latent variable, w");
    }

    // *** Step 2: Binary component aka Update b1 ***
    // 2-0: Store the previous b1 before update
    Eigen::VectorXd b1Prev = ctr->b1;
    // 2-1: Update Omega1 ~ PG(1, eta1)
    // Rcout << "Step 2-1: Sampling Omega1 \n";
    ctr->omega1 = rcpp_pgdraw(ctr->ones, eta1);     // Sample PG
    ctr->Omega1 = (ctr->omega1).asDiagonal();       // Omega1 matrix

    // 2-2: Update Vg1 and z1
    // Vg1
    // Rcout << "Step 2-2: Updating Vg1 \n";
    Eigen::MatrixXd VgInv1(ctr->pZ, ctr->pZ); 
    VgInv1.triangularView<Eigen::Lower>() = ctr->Z.transpose() * (ctr->Omega1) * (ctr->Z); // Zt * Omega * Z : (pxp)
    VgInv1.diagonal().array() += 1 / 100.0; // ZT*Omega*Z  + c*I(prior) where c = 100
    VgInv1.triangularView<Eigen::Upper>() = VgInv1.transpose().eval(); 
    ctr->Vg1.triangularView<Eigen::Lower>() = VgInv1.inverse();
    ctr->Vg1.triangularView<Eigen::Upper>() = ctr->Vg1.transpose().eval();   
    ctr->VgChol1 = ctr->Vg1.llt().matrixL();

    // z1
    // Rcout << "Step 2-2: Updating z1 \n";
    ctr->z1 = ((ctr->w).array() - 0.5).array() / ctr->omega1.array(); 

    // 2-3: Sample b1
    // Rcout << "Step 2-3: Sampling b1 \n";
    // Mean of b1
    Eigen::VectorXd R1 = ctr->z1 - (ctr->areaA * ctr->spPhi);
    const Eigen::VectorXd ZR1 = ctr->Z.transpose() * (ctr->Omega1) * (R1); // (px1)(nxn)(nx1) = (px1)(no DLM effects)
    ctr->b1 = ctr->Vg1 * ZR1; // Mean (pxp)(px1)

    // Variance of b1 using cholesky
    ctr->b1.noalias() += ctr->VgChol1 * as<Eigen::VectorXd>(rnorm(ctr->pZ, 0, sqrt(ctr->sigma2)));

    // *** Extra Step 2: Spatial random effect ***
    if(ctr->spatial){ // For spatial analysis, add a random effect
      // [Update spTau with full conditional (Strictly positive)] --------------------------------------------------------
      double spTauP = ctr->spTau + R::rnorm(0, 1);
      while(spTauP < 0){ // If negative, re-sample
        spTauP = ctr->spTau + R::rnorm(0, 1);
      }

      // 2. Calcuate the full conditional of spTau : p(spTau) * p(spPhi|spTau) 
      double spTauMHratio = 0;

      // Add prior probability: Prior is set as Gamma(3, 3)
      spTauMHratio += R::dgamma(spTauP, 3, 3, true);      // Proposed (Numerator)
      spTauMHratio -= R::dgamma(ctr->spTau, 3, 3, true);  // Current (Denominator)

      for (int spNodeIndex = 0; spNodeIndex < (ctr->spNodes1).size(); spNodeIndex++){
        // Converting from R indexing to C++ index
        int spNode1 = (ctr->spNodes1)[spNodeIndex] - 1;
        int spNode2 = (ctr->spNodes2)[spNodeIndex] - 1;

        // Numerator (Proposed)
        spTauMHratio += - (spTauP) * 0.5 * (pow((ctr->spPhi)[spNode1], 2) - 2 * (ctr->rho) * (ctr->spPhi)[spNode1] * (ctr->spPhi)[spNode2] + pow((ctr->spPhi)[spNode2], 2)); 
        // Denominator (Current)
        spTauMHratio -= - (ctr->spTau) * 0.5 * (pow((ctr->spPhi)[spNode1], 2) - 2 * (ctr->rho) * (ctr->spPhi)[spNode1] * (ctr->spPhi)[spNode2] + pow((ctr->spPhi)[spNode2], 2));
      }

      spTauMHratio = std::min(1.0, exp(spTauMHratio));

      // 3. Accept/Reject
      // Proposed probability is higher -> accept with 1
      if(R::runif(0, 1) < spTauMHratio){
        ctr->spTau = spTauP;
      } // spTau update end

      // ctr->spTau = 1; // Fix spTau as 1 for now.

      // [Update rho with full conditional] ----------------------------------------------------------------
      // Metropolis Hasting for rho 
      double rhoP = ctr->rho + R::rnorm(0, 0.1); // Proposal: q(rho*|rho) ~ N(rho, 0.1) -> Transition in MH ratio cancels out
      while(rhoP < 0 || rhoP > 1){               // rho must be between 0 and 1
        rhoP = ctr->rho + R::rnorm(0, 0.1);
      }
      
      // 2. Calcuate the full conditional of rho : p(rho) * p(spPhi|rho)
      // FC of rho = log(dnorm(0, 0.2)) + log(prob of phi1^2 - 2 * rho * phi1 * phi2 + phi2^2)
      double rhoMHratio = 0;

      // Add prior probability in log: Prior is set as Unif(0, 1) = Beta(1, 1) as rho is between -1 and 1
      rhoMHratio += log(0.5);        // Proposed
      rhoMHratio -= log(0.5);        // Current

      for (int spNodeIndex = 0; spNodeIndex < (ctr->spNodes1).size(); spNodeIndex++){
        // Converting from R indexing to C++ index
        int spNode1 = (ctr->spNodes1)[spNodeIndex] - 1;
        int spNode2 = (ctr->spNodes2)[spNodeIndex] - 1;

        // Proposed
        rhoMHratio += - (ctr->spTau) * 0.5 * (pow(ctr->spPhi[spNode1], 2) - 2 * (rhoP) * ctr->spPhi[spNode1] * ctr->spPhi[spNode2] + pow(ctr->spPhi[spNode2], 2));
        // Current
        rhoMHratio -= - (ctr->spTau) * 0.5 * (pow(ctr->spPhi[spNode1], 2) - 2 * (ctr->rho) * ctr->spPhi[spNode1] * ctr->spPhi[spNode2] + pow(ctr->spPhi[spNode2], 2));
      }

      rhoMHratio = std::min(1.0, exp(rhoMHratio));

      // 3. Accept/Reject
      // Proposed probability is higher -> accept with 1
      if(R::runif(0, 1) < rhoMHratio){
        ctr->rho = rhoP;
      } // rho update end

      // [With updated rho and spTau, update Qinv] ------------------------------------------------
      ctr->areaQ = (ctr->spTau) * (ctr->areaD - (ctr->rho) * ctr->areaW);
      //ctr->areaQinv = (ctr->areaQ).inverse();

      // Rcout << "Vp and zPhi \n";
      // 2-2: Update Vp and zPhi

      // zPhi
      // Rcout << "Step 2-2: Updating z1 \n";
      // Recalculate eta1 with the updated b1
      // eta1.setZero();
      // eta1 = (ctr->Z) * (ctr->b1);                        // Z * \gamma(updated)
      // eta1.noalias() += ctr->areaA * ctr->spPhi;          // Z * \gamma(updated) + spatial random effect

      // Resample PG variable with the updated eta1
      ctr->omegaPhi = ctr->omega1;
      // ctr->omegaPhi = rcpp_pgdraw(ctr->ones, eta1);       // Sample PG
      ctr->OmegaPhi = (ctr->omegaPhi).asDiagonal();       // OmegaPhi matrix
      ctr->zPhi = ((ctr->w).array() - 0.5).array() / ctr->omegaPhi.array(); 

      // Vp
      // Rcout << "Step 2-2: Updating Vg1 \n";
      Eigen::MatrixXd VpInv(ctr->spN, ctr->spN); 
      // At * Omega * A : (spN x spN)
      VpInv.triangularView<Eigen::Lower>() = ctr->areaA.transpose() * (ctr->OmegaPhi) * (ctr->areaA);
      VpInv += ctr->areaQ; // AT*OmegaPhi*A  + CAR prior inverse
      VpInv.triangularView<Eigen::Upper>() = VpInv.transpose().eval(); 
      ctr->Vp.triangularView<Eigen::Lower>() = VpInv.inverse();
      ctr->Vp.triangularView<Eigen::Upper>() = ctr->Vp.transpose().eval();   
      ctr->VpChol = ctr->Vp.llt().matrixL();

      // 2-3: Sample phi
      //Rcout << "Step 2-3: Sampling spPhi \n";
      // Mean of spPhi
      Eigen::VectorXd Rphi = ctr->zPhi - ((ctr->Z) * (b1Prev));
      const Eigen::VectorXd ARphi = ctr->areaA.transpose() * (ctr->OmegaPhi) * (Rphi); // (px1)(nxn)(nx1) = (px1)(no DLM effects)
      ctr->spPhi = ctr->Vp * ARphi; // Mean (spN x spN)(spN x 1)

      // Variance of spPhi using cholesky
      ctr->spPhi.noalias() += ctr->VpChol * as<Eigen::VectorXd>(rnorm(ctr->spN, 0, sqrt(ctr->sigma2)));
    } // Spatial finish

    // Rcout << ctr->rho;
    // Rcout << "\n";
    // Rcout << ctr->spatial;
    // Rcout << "\n";

    // *** Step 3: Negative Binomial (count) component ***
    // 3-1: Update Omega2 ~ PG(y + r, psi2) with dlm effect
    // Rcout << "Step 3-1: Sampling Omega2 \n";
    // Start omega2 vector with ones and sample only for the at-risk (No need to sample PG(y + r, 0))
    (ctr->omega2).setOnes(); // Can set anything as it will be multiplied to zero
    for(int k = 0; k < (ctr->nStar); k++){
      int tmp = (ctr->atRiskIdx)[k];
      (ctr->omega2)[tmp] = rcpp_pgdraw((ctr->Y0)[tmp] + ctr->r, eta2[tmp]);
    }
    ctr->Omega2 = (ctr->omega2).asDiagonal();                 // Omega2 matrix: (nxn)

    // 3-2: Update Zstar, Zw, Vg, z2, R
    // Compute Zstar (Z with non At-risk individuals zeroed out)
    ctr->Zstar = (ctr->Z).array().colwise() * (ctr->w).array();

    // Zw
    // Rcout << "Step 3-2: Updating Zw \n";
    ctr->Zw = ctr->Omega2 * ctr->Zstar;  // (nxn) x (nxp) = (nxp) with only with at-risk

    // Vg (pxp)
    // Rcout << "Step 3-2: Updating Vg \n";
    Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ); // (pxp)
    VgInv.triangularView<Eigen::Lower>() = ctr->Zstar.transpose() * (ctr->Zw); // Zt * Omega * Z: (pxn)(nxn)(nxp)
    VgInv.diagonal().array() += 1 / 100.0; // Zt*Omega*Z  + c*I where c = 100000
    VgInv.triangularView<Eigen::Upper>() = VgInv.transpose().eval(); 
    ctr->Vg.triangularView<Eigen::Lower>() = VgInv.inverse();
    ctr->Vg.triangularView<Eigen::Upper>() = ctr->Vg.transpose().eval();   
    ctr->VgChol = ctr->Vg.llt().matrixL();

    // z2 (= Ystar) (nx1)
    // Rcout << "Step 3-2: Updating z2 \n";
    ctr->z2 = (ctr->Y0 - ctr->rVec).array() / (2*(ctr->omega2).array()).array();
    ctr->Ystar = (ctr->z2).array() * (ctr->w).array();

    // R (partial residual) (nx1)
    // Rcout << "Step 3-2: Updating R \n";
    Eigen::VectorXd fhatStar = ctr->fhat.array() * ctr->w.array(); 

    ctr->R = ctr->Ystar - fhatStar; // (nx1) - (nx1) = (nx1)

    // 3-3: Sample b2
    // Rcout << "Step 3-3: Sampling b2 \n";
    const Eigen::VectorXd ZR = ctr->Zw.transpose() * (ctr->R); // R = Y - f(DLM effects) : (pxn) x (nx1) = (px1)
    ctr->b2 = ctr->Vg * ZR; // Mean: (pxp) x (px1) = (px1)
    ctr->b2.noalias() += ctr->VgChol * as<Eigen::VectorXd>(rnorm(ctr->pZ, 0, sqrt(ctr->sigma2))); // Variance

    // *** Step 4: Update the dispersion parameter, r ***
    // Rcout << "Step 4: Updating dispersion parameter, r \n";
    // Propose a new r with a random walk
    int rC = ctr->r;
    int rP;
    if(R::runif(0, 1) < 0.5){ // Walk to left
      rP = rC - 1;
    } else { // Walk to right
      rP = rC + 1;
    }

    // Support range = 1 ~ 10
    if(rP == 0){
      rP++;
    }

    if(rP == 11){
      rP--;
    }

    // Reset & Compute MH Ratio
    // Rcout << "Step 3: Resetting MHratio \n";
    ctr->MHratio = 0;

    // Rcout << "Step 3: MHratio calculation \n";
    for(int q = 0; q < (ctr->nStar); q++){
      int index_aR = (ctr->atRiskIdx)[q];
      (ctr->MHratio) += R::dnbinom((ctr->Y0)[index_aR], rP, nu2[index_aR], true);
      (ctr->MHratio) -= R::dnbinom((ctr->Y0)[index_aR], rC, nu2[index_aR], true);
    }

    //Rcout << "MH ratio: " << exp(ctr->MHratio) << "\n";
    ctr->MHratio = std::min(1.0, exp(ctr->MHratio));
    
    // Accept/Reject 
    if(R::runif(0, 1) < ctr->MHratio){
      ctr->r = rP;
      ctr->rVec = (ctr->ones).array() * (ctr->r);
    }
    
    // End ZINB
  }
}   // end tdlmModelEst function


// void tdlmModelEstBinomial(modelCtr *ctr)
// {
//   ctr->gamma = ctr->Vg.transpose() * ctr->Zw.transpose() * ctr->R;
//   ctr->VgChol = ctr->Vg.llt().matrixL();
//   ctr->gamma.noalias() += 
//     ctr->VgChol * as<Eigen::VectorXd>(rnorm(ctr->pZ, 0, 1));

//   Eigen::VectorXd psi = ctr->fhat;
//   psi.noalias() += ctr->Z * ctr->gamma;
//   ctr->Omega = rcpp_pgdraw(ctr->binomialSize, psi);
//   ctr->Zw = ctr->Omega.asDiagonal() * ctr->Z;
  
//   Eigen::MatrixXd VgInv(ctr->pZ, ctr->pZ);
//   VgInv.triangularView<Eigen::Lower>() = ctr->Z.transpose() * ctr->Zw;
//   VgInv.diagonal().array() += 1 / 100000.0;
//   VgInv.triangularView<Eigen::Upper>() = VgInv.transpose().eval();
//   ctr->Vg.triangularView<Eigen::Lower>() = VgInv.inverse();
//   ctr->Vg.triangularView<Eigen::Upper>() = ctr->Vg.transpose().eval();
//   ctr->Lambda = (ctr->kappa).array() / ctr->Omega.array();
// }

/**
 * @brief Construct a new progress Meter::progress Meter object
 * 
 * @param c model control data to track meter progress
 */
progressMeter::progressMeter(modelCtr* c)
{
  ctr = c;
  startTime = time(NULL);
  if (ctr->verbose)
    Rcout << "Burn-in % complete \n" <<
      "[0--------25--------50--------75--------100]\n '";
  burnProgInc = (ctr->burn / 42.0);
  burnProgMark = burnProgInc;
  iterProgInc = (ctr->iter / 42.0);
  iterProgMark = double(ctr->burn) + iterProgInc;
}
progressMeter::~progressMeter()
{
  ctr = 0;
}
/**
 * @brief print progress mark "'" at given intervals thorughout model run 
 */
void progressMeter::printMark()
{
  if (ctr->verbose) {
    if (ctr->b > ctr->burn) {
      if (ctr->b >= iterProgMark) {
        Rcout << "'"; 
        iterProgMark += iterProgInc;
      }
    } else {
      if (ctr->b >= burnProgMark) {
        Rcout << "'"; 
        burnProgMark += burnProgInc;
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
}



double tdlmProposeTree(Node* tree, exposureDat* Exp, modelCtr* ctr, int step)
{
  int no = 0;
  double stepMhr = 0;
  std::vector<Node*> dlnmTerm, tempNodes;

  // List current tree terminal nodes
  dlnmTerm = tree->listTerminal();

  // Grow
  if (step == 0) {
    // select node to grow
    no = (std::size_t) floor(R::runif(0, dlnmTerm.size())); // Uniform selection among terminal nodes

    if (dlnmTerm[no]->grow()) { // propose new split
      double nGen2 = double(tree->nGen2());
      if (dlnmTerm[no]->depth == 0) { // If selected terminal node is the root node,
        ++nGen2;
      } else {
        if (!(dlnmTerm[no]->parent->isGen2())) {
          ++nGen2;
        }
      }

      // Compute MH ratio for growing
      stepMhr = log((double)tree->nTerminal()) - log(nGen2) +
        2 * logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth + 1, 1) +
        logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth, 0) -
        logPSplit((ctr->treePrior)[0], (ctr->treePrior)[1], dlnmTerm[no]->depth, 1);

      Exp->updateNodeVals((dlnmTerm[no]->proposed)->c1); // update node values
      // newDlnmTerm = tree->listTerminal(1); // list proposed terminal nodes
    }


  // Prune
  } else if (step == 1) {
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
    tempNodes = tree->listInternal();
    no = floor(R::runif(0, tempNodes.size())); // select internal nodes to change 
    if (tempNodes[no]->change()) { // propose new split
      for (Node* tn : tempNodes[no]->proposed->listTerminal())
        Exp->updateNodeVals(tn);
      // Rcout << "!";
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
} // end tdlmProposeTree function


double modProposeTree(Node* tree, modDat* Mod, dlmtreeCtr* ctr, int step)
{
  int no = 0;
  double stepMhr = 0;
  std::vector<Node*> modTerm, tempNodes;


  // List current tree terminal nodes
  modTerm = tree->listTerminal();
  
  // Grow
  if (step == 0) {
    // select node to grow
    no = (std::size_t) floor(R::runif(0, modTerm.size())); 

    if (modTerm[no]->grow()) { // propose new split
      double nGen2 = double(tree->nGen2());
      if (modTerm[no]->depth == 0) { // depth == 0
        ++nGen2;
      } else {
        if (!(modTerm[no]->parent->isGen2())) {
          ++nGen2;
        }
      }
      
      stepMhr = log((double)tree->nTerminal()) - log(nGen2) +
        2 * logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                      modTerm[no]->depth + 1, 1) +
        logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                  modTerm[no]->depth, 0) -
        logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                  modTerm[no]->depth, 1);

      Mod->updateNodeVals((modTerm[no]->proposed)->c1); // update node values
      if ((modTerm[no]->proposed)->c1->nodevals->idx.size() == 0 ||
          (modTerm[no]->proposed)->c2->nodevals->idx.size() == 0) {
        tree->reject();
        return(0);
      } // end reject if empty
    } // end grow proposal

  // Prune
  } else if (step == 1) {
    tempNodes = tree->listGen2();
    no = floor(R::runif(0, tempNodes.size())); // select gen2 node to prune

    stepMhr = log((double)tree->nGen2()) - log((double)tree->nTerminal() - 1.0) -
      2 * logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                    tempNodes[no]->depth + 1, 1) -
      logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                tempNodes[no]->depth, 0) +
      logPSplit((ctr->treePriorMod)[0], (ctr->treePriorMod)[1],
                tempNodes[no]->depth, 1);

    tempNodes[no]->prune(); // prune nodes
    if (tempNodes[no]->nodevals->nestedTree != 0)
      delete tempNodes[no]->nodevals->nestedTree;
    tempNodes[no]->nodevals->nestedTree = 0;

  // Change
  } else if (step == 2) {
    tempNodes = tree->listInternal();
    no = floor(R::runif(0, tempNodes.size())); // select internal nodes to change

    if (tempNodes[no]->change()) { // propose new split
      for (Node* tn : tempNodes[no]->proposed->listTerminal()) {
        Mod->updateNodeVals(tn);
        if (tn->nodevals->idx.size() == 0) {
          tree->reject();
          return(0);
        } // end reject if empty
      }

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

  // Swap
  } else {
    tempNodes = tree->listInternal();
    no = floor(R::runif(0, tempNodes.size() - 1)); // dont select top node
    if (tempNodes[no]->parent->swap(tempNodes[no])) {
      for (Node* tn : tempNodes[no]->parent->proposed->listTerminal()) {
        Mod->updateNodeVals(tn);
        if (tn->nodevals->idx.size() == 0) {
          tree->reject();
          return(0);
        } // end reject if empty
      }

      for (Node* tn : tempNodes[no]->parent->listInternal())
        stepMhr -= (tn->nodestruct)->logPRule();
      for (Node* tn : tempNodes[no]->parent->proposed->listInternal())
        stepMhr += (tn->nodestruct)->logPRule();

    }
  }

  return(stepMhr);
} // end modProposeTree function


void dlmtreeRecDLM(dlmtreeCtr* ctr, dlmtreeLog* dgn)
{
  double sumDLM;
  for (int i = 0; i < ctr->n; ++i) {
    dgn->exDLM.col(i) += ctr->exDLM.col(i);
    dgn->ex2DLM.col(i) += ctr->exDLM.col(i).array().square().matrix();
    sumDLM = ctr->exDLM.col(i).sum();
    dgn->cumDLM(i) += sumDLM;
    dgn->cum2DLM(i) += pow(sumDLM, 2);
  }
} // end dlmtreeRecDLM function

std::string modRuleStr(Node* n, modDat* Mod)
{
  std::string rule = "";
  if (n->depth == 0)
    return(rule);
  
  Node* parent = n->parent;
  int splitVar = parent->nodestruct->get(1);
  int splitVal = parent->nodestruct->get(2);
  std::vector<int> splitVec = parent->nodestruct->get2(1);
  
  rule += std::to_string(splitVar);
  if (Mod->varIsNum[splitVar]) {
    if (parent->c1 == n)
      rule += "<";
    else
      rule += ">=";
    rule += std::to_string(splitVal);
  } else {
    if (parent->c1 == n)
      rule += "[]";
    else
      rule += "][";
    for (int i : splitVec)
      rule += std::to_string(i) + ",";
    rule.pop_back();
  }
  if (parent->depth != 0)
    rule += "&" + modRuleStr(parent, Mod);
  
  return(rule);
} // end modRuleStr function
                 
                 
Eigen::VectorXd countMods(Node* tree, modDat* Mod)
{
  Eigen::VectorXd modCount(Mod->nMods); modCount.setZero();
  Eigen::VectorXd unavailProb(Mod->nMods); unavailProb.setZero();
  std::vector<int> unavail;  
  for (Node* tn : tree->listInternal()) {
    modCount(tn->nodestruct->get(1)) += 1.0;
    unavail.clear();
    unavailProb.setZero();
    for (int i = 0; i < Mod->nMods; ++i) {
      if (tn->nodestruct->get3(1)[i].size() == 0) {
        unavail.push_back(i);
        unavailProb(i) = Mod->modProb[i];
      }
    }
    if (unavail.size() > 0) {
      std::random_shuffle(unavail.begin(), unavail.end());
      double totProb = unavailProb.sum();
      int pseudoDraw = R::rgeom(std::max(0.00000001, 1 - totProb));
      int binomDraw = 0;
      if (pseudoDraw > 0) {
        for (int i : unavail) {
          binomDraw = R::rbinom(pseudoDraw, unavailProb(i) / totProb);
          if (binomDraw > 0)
            modCount(i) += binomDraw * 1.0;
          totProb -= unavailProb(i);
          pseudoDraw -= binomDraw;
          if (pseudoDraw < 1)
            break;
        } // end multinom
      } // end pseudoDraw
    } // end unavail
  } // end modCount
  return(modCount);
} // end countMods function

void drawTree(Node* tree, Node* n, double alpha, double beta)
{
  double logProb = log(alpha) - beta * log(1.0 + n->depth);
  if (log(R::runif(0, 1)) < logProb) {
    if (n->grow()) {
      if (n->depth > 0)
        n = n->proposed;
        
      tree->accept();
      drawTree(tree, n->c1, alpha, beta);
      drawTree(tree, n->c2, alpha, beta);
    }
  } // end grow tree
  return;
} // end drawTree function

void updateGPMats(Node* n, dlmtreeCtr* ctr)
{
  if (n->nodevals->updateXmat == 0)
    return;
  if (n->depth == 0) {
    n->nodevals->XtX = ctr->XtXall;
    n->nodevals->ZtXmat = ctr->ZtXall;
    n->nodevals->VgZtXmat = ctr->VgZtXall;
    n->nodevals->updateXmat = 0;
    return;
  }
  Node* par = n->parent;
  if (par->nodevals->updateXmat)
    updateGPMats(par, ctr);
  Node* sib = n->sib();
  std::vector<int> idx;
  if (n->nodevals->idx.size() <= sib->nodevals->idx.size()) {
    idx = n->nodevals->idx;
  } else {
    idx = sib->nodevals->idx;
  }
  
  Eigen::MatrixXd Xtemp(idx.size(), ctr->pX); Xtemp.setZero();
  Eigen::MatrixXd Ztemp(idx.size(), ctr->pZ); Ztemp.setZero();
  
  for (std::size_t i = 0; i < idx.size(); ++i) {
    Xtemp.row(i) = ctr->X.row(idx[i]);
    Ztemp.row(i) = ctr->Z.row(idx[i]);
  }
  
  if (n->nodevals->idx.size() <= sib->nodevals->idx.size()) {
    n->nodevals->XtX = Xtemp.transpose() * Xtemp;
    n->nodevals->ZtXmat = Ztemp.transpose() * Xtemp;
    n->nodevals->VgZtXmat = ctr->Vg * n->nodevals->ZtXmat;
    sib->nodevals->XtX = par->nodevals->XtX - n->nodevals->XtX;
    sib->nodevals->ZtXmat = par->nodevals->ZtXmat - n->nodevals->ZtXmat;
    sib->nodevals->VgZtXmat = par->nodevals->VgZtXmat - n->nodevals->VgZtXmat;
  } else {
    sib->nodevals->XtX = Xtemp.transpose() * Xtemp;
    sib->nodevals->ZtXmat = Ztemp.transpose() * Xtemp;
    sib->nodevals->VgZtXmat = ctr->Vg * sib->nodevals->ZtXmat;
    n->nodevals->XtX = par->nodevals->XtX - sib->nodevals->XtX;
    n->nodevals->ZtXmat = par->nodevals->ZtXmat - sib->nodevals->ZtXmat;
    n->nodevals->VgZtXmat = par->nodevals->VgZtXmat - sib->nodevals->VgZtXmat;
  }
  n->nodevals->updateXmat = 0;
  sib->nodevals->updateXmat = 0;
  
}