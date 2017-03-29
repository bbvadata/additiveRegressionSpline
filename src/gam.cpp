#include <RcppArmadillo.h>
using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::vec hatMatrix(arma::mat m, int length) {
  // Calculate the influence vector
  arma::vec diag = arma::diagvec(m * arma::inv(m.t()*m) * m.t());
  // Return first length values from the diagonal
  return(diag.subvec(0,length-1));
}

//[[Rcpp::export]]
arma::mat cubicSpline(NumericVector x, NumericVector knots) {
  // Create an empty matrix plus two columns:
  // one for intercept and another for x values
  arma::mat modelMatrix = arma::ones(x.length(), knots.length()+2);
  // Update second column with x values
  modelMatrix.col(1) = Rcpp::as<arma::colvec>(x);
  // Iteration over columns
  for(int i=0; i<knots.length(); i++) {
    // Apply on first knot
    double z = knots[i];
    // Vector to save one basis i.e.: one column of design matrix
    NumericVector basis(x.length());
    // Define iterators
    NumericVector::iterator it, out_it;
    for(it = x.begin(), out_it=basis.begin();it<x.end();++it,++out_it) {
      // Cubic spline
      *out_it = ((std::pow(z-0.5,2) - (1.0/12.0)) * (std::pow(*it-0.5,2) - (1.0/12.0))) /4 -
        (std::pow(std::abs(*it-z)-0.5,4) - (std::pow(std::abs(*it-z)-0.5,2)/2) + (7.0/240.0)) / 24.0;
    }
    // Update column from model matrix
    modelMatrix.col(i+2) = Rcpp::as<arma::colvec>(basis);
  }
  return(modelMatrix);
}

//[[Rcpp::export]]
List splineModel(NumericVector x, NumericVector knots, arma::mat y) {
  // Model matrix with cubic basis
  arma::mat modelMatrix = cubicSpline(x, knots);
  // Solve system of equations (fast mode)
  arma::mat betas = arma::solve(modelMatrix, y);
  // Return results
  List results;
  results["betas"] = betas;
  results["modelMatrix"] = modelMatrix;
  return(results);
}

//[[Rcpp::export]]
List splineModelPenalized(arma::vec y, NumericVector x, NumericVector knots, double lambda) {
  // Model matrix
  arma::mat modelMatrixKnots = arma::zeros(knots.length()+2, knots.length()+2);
  // Cubic splines on knots modelMatrixKnots
  arma::mat cubicS = cubicSpline(knots,knots);
  modelMatrixKnots.submat(2,2,modelMatrixKnots.n_rows-1,modelMatrixKnots.n_cols-1) = cubicS.submat(0,2,cubicS.n_rows-1,cubicS.n_cols-1);

  // Model matrix
  arma::mat  modelMatrix = cubicSpline(x, knots);
  // Eign decomposition
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec,modelMatrixKnots,"std");
  eigval = arma::pow(eigval, 0.5);
  // Replace NAN for 0
  eigval.replace(arma::datum::nan, 0);
  // Augmented model matrix
  arma::mat mS = eigvec * arma::diagmat(eigval) * eigvec.t();
  mS = mS * std::pow(lambda, 0.5);
  // Join two matrix
  arma::mat mSA = arma::join_cols(modelMatrix, mS);
  // Augmented y vector
  arma::mat yAug = arma::join_cols(y, arma::zeros(knots.length()+2));
  // Solve system of equations (fast mode)
  arma::mat betas;
  try {
    betas = arma::solve(mSA, yAug);
  } catch(...) {
    ::Rf_error("Solve can not find a solution for the system.");
  }
  // General Crodd Validation (GCV)) calculation
  arma::vec influence = hatMatrix(mSA, x.length());
  // Fitted values
  arma::vec fitted = mSA * betas;
  int n = x.length();
  fitted = fitted.subvec(0, n-1);
  double rss = sum(arma::pow(y - fitted, 2));
  double gcv = (n*rss) / std::pow((n-sum(influence)), 2);

  // Return results
  List results;
  results["betas"] = betas;
  results["modelatrix"] = mSA;
  results["rss"] = rss;
  results["gcv"] = gcv;
  return(results);
}



//[[Rcpp::export]]
List splineModelPMatrix(NumericVector x, NumericVector z,
                            NumericVector knots_x, NumericVector knots_z) {
  int q = (2*(knots_x.length()+2))-1;
  // Penalty matrix
  // Get penalties matrix 1
  arma::mat sCX = arma::zeros(knots_x.length()+2, knots_x.length()+2);
  arma::mat cubicSX = cubicSpline(knots_x,knots_x);
  sCX.submat(2,2,sCX.n_rows-1,sCX.n_cols-1) = cubicSX.submat(0,2,cubicSX.n_rows-1,cubicSX.n_cols-1);
  // Update matrix S with penalty matrix 1
  arma::mat sX = arma::zeros(q,q);
  sX.submat(1,1,knots_x.length()+1,knots_x.length()+1) = sCX.submat(1,1,sCX.n_rows-1,sCX.n_cols-1);

  // Get penalties matrix 2
  arma::mat sCZ = arma::zeros(knots_z.length()+2, knots_z.length()+2);
  arma::mat cubicSZ = cubicSpline(knots_z,knots_z);
  sCZ.submat(2,2,sCZ.n_rows-1,sCZ.n_cols-1) = cubicSZ.submat(0,2,cubicSZ.n_rows-1,cubicSZ.n_cols-1);
  // Update matrix S with penalty matrix 2
  arma::mat sZ = arma::zeros(q,q);
  sZ.submat(knots_x.length()+2,knots_x.length()+2,sZ.n_cols-1,sZ.n_rows-1) = sCZ.submat(1,1,sCZ.n_rows-1,sCZ.n_cols-1);

  // Model matrix
  int n = x.length();
  arma::mat modelMatrix = arma::ones(n,q);
  // Cubic splines smooths over knots
  arma::mat xSmooth = cubicSpline(x,knots_x);
  arma::mat zSmooth = cubicSpline(z,knots_z);
  // Update model matrix with smooths (we remove the intercept for xSmooth and zSmooth i.e.: dropping the first column)
  modelMatrix.submat(0,1,modelMatrix.n_rows-1,xSmooth.n_cols-1) = xSmooth.submat(0,1,xSmooth.n_rows-1,xSmooth.n_cols-1);
  modelMatrix.submat(0,zSmooth.n_cols,modelMatrix.n_rows-1,modelMatrix.n_cols-1) = zSmooth.submat(0,1,zSmooth.n_rows-1,zSmooth.n_cols-1);

  // Return results
  List results;
  results["modelMatrix"] = modelMatrix;
  results["sX"] = sX;
  results["sZ"] = sZ;
  results["sCX"] = sCX;
  results["sCZ"] = sCZ;
  results["cubicSX"] = cubicSX;
  results["cubicSZ"] = cubicSZ;
  results["xSmooth"] = xSmooth;
  results["zSmooth"] = zSmooth;
  return(results);
}

//[[Rcpp::export]]
List splineModelPMatrixFit(arma::vec y, NumericVector x,
                           arma::mat modelMatrix,
                           arma::mat sX,arma::mat sZ,
                           double lambda_x, double lambda_z) {
  int n = x.length();
  // Fit model
  arma::mat penaltyM = sX * lambda_x + sZ * lambda_z;
  // Sqrt of penalty matrix
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, penaltyM,"std");
  eigval = arma::pow(eigval, 0.5);
  // Replace NAN for 0
  eigval.replace(arma::datum::nan, 0);
  // Augmented model matrix
  arma::mat penaltyMSqrt = eigvec * arma::diagmat(eigval) * eigvec.t();
  // Augmented ModelMatrix: Row join between modelMatrix (with smoothings) and penalty matrix
  arma::mat augmentedModelMatrix = arma::join_cols(modelMatrix, penaltyMSqrt);
  // Augmented Y vector (with number of paramters)
  arma::vec yAugmented = arma::zeros(y.n_rows+modelMatrix.n_cols);
  yAugmented.subvec(0,y.n_rows-1) = y;
  // Solve system of equations (fast mode)
  arma::mat betas;
  try {
    betas = arma::solve(augmentedModelMatrix, yAugmented);
  } catch(...) {
    ::Rf_error("Solve can not find a solution for the system.");
  }
  // General Crodd Validation (GCV)) calculation
  arma::vec influence = hatMatrix(augmentedModelMatrix, x.length());
  // Fitted values
  arma::vec fitted = augmentedModelMatrix * betas;
  fitted = fitted.subvec(0, n-1);
  double rss = sum(arma::pow(y - fitted, 2));
  double gcv = (n*rss) / std::pow((n-sum(influence)), 2);

  // Return results
  List results;
  results["modelMatrix"] = modelMatrix;
  results["penaltyM"] = penaltyM;
  results["penaltyMSqrt"] = penaltyMSqrt;
  results["augmentedModelMatrix"] = augmentedModelMatrix;
  results["betas"] = betas;
  results["gcv"] = gcv;
  results["rss"] = rss;
  results["lambda_x"] = lambda_x;
  results["lambda_z"] = lambda_z;
  results["fitted"] = fitted;
  return(results);
}
