#include <RcppGSL.h>
#include <gsl/gsl_linalg.h>

using namespace Rcpp;
using namespace std;

double correlationCoefficient(NumericVector X, NumericVector Y, int n){
  
  double sum_X = 0, sum_Y = 0, sum_XY = 0;
  double squareSum_X = 0, squareSum_Y = 0;
  for (int i = 0; i < n; i++){
    sum_X = sum_X + X[i];
    sum_Y = sum_Y + Y[i];
    sum_XY = sum_XY + X[i] * Y[i];
    squareSum_X = squareSum_X + X[i] * X[i];
    squareSum_Y = squareSum_Y + Y[i] * Y[i];
  }
  double corr = (double)(n * sum_XY - sum_X * sum_Y) / sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
  return corr;
}

void makeSubSigma(NumericVector snp_50, NumericVector sigma_50, int subSize, int sampleSize){
  NumericVector X (sampleSize);
  NumericVector Y (sampleSize);
  
  for(int i=0;i<subSize; i++){
    for(int j=0; j<subSize; j++){
      for(int k=0; k<sampleSize; k++){
        X[k]=snp_50[i*sampleSize+k];
        Y[k]=snp_50[j*sampleSize+k];
      }
      sigma_50[i*subSize+j]=correlationCoefficient(X,Y,sampleSize);
    }
  }
}

void makeSigmaPositiveSemiDefinite(NumericVector sigma, int size) {
  int gsl_tmp = 0;
  double matDet  = 0;
  double addDiag = 0;
  bool positive = false;
  // int z = 0;
  gsl_matrix * tmpResultMatrix = gsl_matrix_calloc (size, size);
  gsl_permutation *p = gsl_permutation_alloc(size);
  do{
    for(int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        if(i==j)
          gsl_matrix_set(tmpResultMatrix,i,j,sigma[i*size+j]+addDiag);
        else
          gsl_matrix_set(tmpResultMatrix,i,j,sigma[i*size+j]);
      }
    }
    
    gsl_linalg_LU_decomp(tmpResultMatrix, p, &gsl_tmp);
    matDet = gsl_linalg_LU_det(tmpResultMatrix,gsl_tmp);
    // std::cout.precision(11);
    // std::cout << "matDet:" << matDet << "\n";
    
    if(matDet > 0 )
      positive = true;
    else {
      // z+=1;
      // std::cout << "number addDiag:" << z << "\n";
      addDiag+=0.1;
    }
  } while(!positive);
  
  // std::cout << "addDiag:" << addDiag << "\n";
  for(int i = 0; i < size*size; i++){
    if(i%(size+1) == 0)
      sigma[i] = sigma[i] + addDiag;
  }
}


// [[Rcpp::export]]
NumericMatrix makesigmasemiPDRcppGSL(NumericMatrix geno) {
  int N = geno.cols();
  int subsize = geno.rows();
  
  Function c("c");   
  
  NumericVector raw_sigmaMatrix (N * subsize);
  NumericVector snp_50 = c(transpose(geno));
  makeSubSigma(snp_50, raw_sigmaMatrix, subsize, N);
  makeSigmaPositiveSemiDefinite(raw_sigmaMatrix, subsize);
  
  NumericMatrix sigmaMatrix(subsize, subsize, raw_sigmaMatrix.begin() );
  
  return(sigmaMatrix);
}
