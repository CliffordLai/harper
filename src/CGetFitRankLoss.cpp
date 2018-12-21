#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

# define M_PI           3.14159265358979323846  /* pi */

// /*Rank for ties "Average" */
// NumericVector Crank(NumericVector x){
//   int N = x.size();
//   NumericVector out(N);
//   NumericVector outrank(N);
//   
//   for(int k=0;k<N;k++){
//     outrank(k) = 1;
//     for(int h=0;h<N;h++){
//       if(h==k) 
//         continue;
//       if(x(h)<x(k))
//         outrank(k)++;
//       if(x(h)==x(k))
//         outrank(k) = outrank(k)+0.5;
//     }
//   }
//   return outrank;
// }

/*Rank for  ties "last" */
NumericVector Crank(NumericVector x){
  int N = x.size();
  NumericVector out(N);
  NumericVector outrank(N);
  
  for(int k=0;k<N;k++){
    outrank(k) = 1;
    for(int h=0;h<N;h++){
      if(h==k) 
        continue;
      if(x(h)<x(k))
        outrank(k)++;
      if(x(h)==x(k))
        if(k<h) 
          outrank(k) = outrank(k)+1;
    }
  }
  return outrank;
}

/*MSE*/
double MSE(NumericVector x, NumericVector y){
  int N = x.size();
  double out = 0;
  for(int i=0; i<N; i++){
    out += pow( (x(i)-y(i)), 2)/N;
  }
  return out;
}


double ParabolaMin(double a, double b, double c,
                                  double fa, double fb, double fc){
  double C1= (b-a)*(fb-fc);
  double C2= (b- c)*(fb- fa);
  double out = b - 0.5 * ((b- a)*C1-(b-c)*C2) / (C1-C2);
  return out;
}


// [[Rcpp::export]] 
NumericMatrix CMGetFitRankLoss(NumericMatrix Z, 
                               NumericVector t,
                               NumericVector lambda, 
                               NumericVector phi,
                               bool Parabola){
  int N= Z.nrow();
  int M = Z.ncol();
  int I = lambda.size();
  int J = phi.size();
  NumericMatrix RZ(N,M);
  NumericMatrix RZfit(I*J,N);
  NumericVector zs(N);
  NumericVector zhat(N);
  NumericMatrix mse(I,J);
  double an = pow((N+1.0),2)/4.0;
  double bn = (N+1.0)*sqrt(N-1.0)/12.0;
  NumericVector minmse(M,  pow((N+1.0),2)/2.0  ); //With a default value pow((N+1.0),2)/2.0
  NumericVector minI(M);
  NumericVector minJ(M);
  NumericMatrix out(3,M);
  
  for(int m=0; m<M; m++){
    zs = Z(_,m);
    RZ(_,m) = Crank(zs);
  }
  
  for(int i=0; i<I; i++){
    for(int j=0; j<J; j++){
      for(int k=0; k<N; k++){
        RZfit(i*J+j,k) = cos(2*M_PI*lambda(i)*t(k)+2*M_PI*phi(j)); 
        RZfit(i*J+j,k)  = round(RZfit(i*J+j,k) *100000.0)/100000.0; //The accurate digits of the sinusoid curve
      }
      zhat = RZfit(i*J+j,_);
      zhat = Crank(zhat);
      for(int k=0; k<N; k++){
        RZfit(i*J+j,k) = zhat(k);
      }
    }
  }
  
  for(int m=0; m<M; m++){
    zs = RZ(_,m);
    for(int i=0; i<I; i++){
      for(int j=0; j<J; j++){
        zhat = RZfit(i*J+j,_);
        mse(i,j) = MSE(zs,zhat);
        if(mse(i,j)<minmse(m) ){
          minmse(m) =  mse(i,j);
          minI(m) = i;
          minJ(m) = j;
        }
      }
    }
  }   
  
  if(!Parabola){
    for(int m=0; m<M; m++){
      out(0,m) = (((N+1.0)*(2*N+1.0)/3.0-minmse(m) )/2.0-an)/bn;
      out(1,m) = lambda(minI(m));
      out(2,m) = phi(minJ(m)); 
    }
  }else{
    double a,b,c,fa,fb,fc; 
    for(int m=0; m<M; m++){
      if(minI(m)>0 && minI(m)<(I-1) ){
        a = lambda(minI(m)-1);
        b = lambda(minI(m));
        c = lambda(minI(m)+1);
        fa = mse(minI(m)-1,minJ(m));
        fb = mse(minI(m),minJ(m));
        fc = mse(minI(m)+1,minJ(m));
      }else if(minI(m)==0){
        a = lambda(minI(m));
        b = lambda(minI(m)+1);
        c = lambda(minI(m)+2);
        fa = mse(minI(m),minJ(m));
        fb = mse(minI(m)+1,minJ(m));
        fc = mse(minI(m)+2,minJ(m));
      }else{
        a = lambda(minI(m)-2);
        b = lambda(minI(m)-1);
        c = lambda(minI(m));
        fa = mse(minI(m)-2,minJ(m));
        fb = mse(minI(m)-1,minJ(m));
        fc = mse(minI(m),minJ(m));
      }
      out(1,m) = ParabolaMin(a,b,c,fa,fb,fc);
      out(2,m) = phi(minJ(m));  
      
      for(int k=0; k<N;k++){
        zhat(k) = cos(2*M_PI*out(1,m)*t(k)+out(2,m)); 
        zhat(k) = round(zhat(k)*10000)/10000.0;
      }
      zs = RZ(_,m);
      zhat = Crank(zhat);
      
      //Check if it is the minimum.
      out(0,m) = MSE(zs,zhat);
      if(out(0,m)>minmse(m)){
        out(0,m) =  (((N+1.0)*(2*N+1.0)/3.0-minmse(m) )/2.0-an)/bn;
        out(1,m) = lambda(minI(m));
      }else{
        out(0,m) =  (((N+1.0)*(2*N+1.0)/3.0-out(0,m) )/2.0-an)/bn;
      }
  
    }
  }
  return(out);
}


