#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix minSSEC(NumericMatrix Q, NumericMatrix Z){
  int J0;
  int n = Z.nrow(), m = Z.ncol();
  int I = n-3;
  int nlambda = Q.ncol();
  nlambda = nlambda/n;
//  int nrowSSE = nlambda;
  int minI;
  
  NumericMatrix out(I,m);  
  NumericMatrix SSE(nlambda,m);
  NumericMatrix minSSE(2,m);  
  
  for(int k=0; k<nlambda; k++){
    for(int i=0; i<I; i++){
      for(int h=0; h<m; h++){
        for(int j=0; j<n; j++){
          J0 = n*k+j;
          out(i,h) += Q(i,J0)*Z(j,h);
        }
      }
    }
    
    for(int h=0; h<m; h++){
      for(int i=0; i<I; i++){
        SSE(k,h) += pow(out(i,h),2);
        out(i,h) = 0;
      }
    }
  }
  
  for(int j=0; j<m; j++){
    minI = 0;
    for(int i=1; i<nlambda; i++){
      if(SSE(minI,j)>SSE(i,j)){
        minI = i;
      }
    }
    minSSE(0,j) = SSE(minI,j);
    minSSE(1,j) = minI+1;
  }
  
  return(minSSE);
}
