#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


// penalized likelihood ratio
double penLR(double x, int len, double pen, int n) {
  double th = double(len)/double(n);
  // double pen = sqrt(2*(1-log(th*(1-th)))); // for the sake of speedup, shift it to the input
  double aux = 2*n*(th*log(th/x) + (1-th)*log((1-th)/(1-x)));
  // return (aux>0?sqrt(aux):0) - pen;
  return sqrt(std::max(0.0,aux)) - pen;
}

// [[Rcpp::export(".msQuantile")]]
NumericVector msQuantile(IntegerVector left, IntegerVector right, int n, int nsim, bool isGen) {
  // output
  NumericVector stat(nsim);
  int nintv = left.size(); // number of intervals
  // // record the last simulation
  // static List last_simulation;  // a list of entries: n, nsim, stat
  
  
  // Rprintf("n = %d\n", n);
  // Rprintf("nsim = %d\n", nsim);
  // Rprintf("nintv = %d\n", nintv);
  
  // if (last_simulation.size() > 0 && n == (int)last_simulation[0] && nsim == (int)last_simulation[1] && nintv == (int)last_simulation[2]) {
  //     if (verbose) {
  //       Rprintf("Quantile is retrieved form the previous simulation!\n");
  //     }
  //     stat = last_simulation[3];
  // } else {
  //     if (verbose)
  //         Rprintf("Quantile simulation (only for the first time) might take a while ... ");
  // pre-store the penality for computational speedup, at expense of memory 
  double aux;
  NumericVector Pen(nintv);
  for (int i = 0; i < nintv; ++i) {
    aux = (right[i]-left[i])/double(n);
    Pen[i] = sqrt(2*(1-log(aux*(1-aux))));
  }
  // Monte-Carlo simulation
  NumericVector U(n);
  double pLR;
  GetRNGstate();
  for (int i = 0; i < nsim; ++i) {
    U = runif(n, 0, 1);
    std::sort(U.begin(), U.end());
    U.push_back(1);
    stat[i] = R_NegInf;
    for (int j = 0; j < nintv; ++j) {
      if (isGen) 
        pLR = std::max(penLR(U[right[j]+1]-U[left[j]],right[j]-left[j],Pen[j],n), penLR(U[right[j]]-U[left[j]+1],right[j]-left[j],Pen[j],n));
      else 
        pLR = penLR(U[right[j]]-U[left[j]], right[j]-left[j], Pen[j], n);
      if (!ISNA(pLR)) 
        stat[i] = stat[i]>pLR ? stat[i] : pLR;// stat[i] = std::fmax(stat[i], pLR);
    }
  }
  PutRNGstate();
  // // store simulated results
  // last_simulation = List::create(n, nsim, nintv, stat);
  // if (verbose)
  //     Rprintf("... end!\n");
  // }
  return stat;
}



