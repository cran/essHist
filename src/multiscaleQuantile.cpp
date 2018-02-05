#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


// penalized likelihood ratio
double penLR(double x, int len, double pen, int n) {
  double th = double(len)/double(n);
  // pen = sqrt(2*(1-log(th*(1-th)))); for the sake of speedup
  return sqrt(fmax(2*n*(th*log(th/x) + (1-th)*log((1-th)/(1-x))), 0)) - pen;
}

// [[Rcpp::export(".msQuantile")]]
NumericVector msQuantile(int n, int nsim) {
  // output
  NumericVector stat(nsim);

  // record the last simulation
  static List last_simulation;  // a list of entries: n, nsim, stat

  if (last_simulation.size() > 0 && n == (int)last_simulation[0] && nsim == (int)last_simulation[1]) {

    // Rprintf("Use the previous simulation!\n");

    stat = last_simulation[2];
  } else {
    Rprintf("Quantile simulation (only for the first time) might take a while ... ");
    
    // Rprintf("Start Monte-Carlo simulation ...\n");

    // compute Rivera & Walther (2013) system of intervals
    int   nintv = 0; // number of intervals
    int    lmin = 2;
    int    lmax = floor(log(n/log(n))/log(2));
    int      dl;
    double   ml;

    for (int l = lmin; l <= lmax; ++l) {
      ml = double(n) * powf(2.0, -l);
      dl = ceil(ml/(6*sqrt(double(l))));
      for (int len = ceil(ml/dl); len <= floor(2*ml/dl); ++len) {
        nintv += floor((n-1)/double(dl)) - len + 1;
      }
    }

    // Rprintf("nintv = %d\n", nintv);

    IntegerVector  St(nintv);
    IntegerVector  Ed(nintv);
    IntegerVector Len(nintv);
    NumericVector Pen(nintv);
    int cnt = 0;
    double aux;
    for (int l = lmin; l <= lmax; ++l) {
      ml = double(n) * powf(2.0, -l);
      dl = ceil(ml/(6*sqrt(double(l))));
      for (int len = ceil(ml/dl); len <= floor(2*ml/dl); ++len) {
        for (int st = 0; (st + len)*dl < n; ++st) {
           St[cnt] = st*dl;
           Ed[cnt] = (st+len)*dl;
          Len[cnt] = len*dl + 1;
          aux = Len[cnt]/double(n);
          Pen[cnt] = sqrt(2*(1-log(aux*(1-aux))));
          ++cnt;
        }
      }
    }

    // Rprintf("cnt = %d\n", cnt);

    // Monte-Carlo simulation
    NumericVector U(n);
    double pLR;
    GetRNGstate();
    for (int i = 0; i < nsim; ++i) {
      U = runif(n, 0, 1);
      sort(U.begin(), U.end());
      stat[i] = R_NegInf;
      for (int j = 0; j < nintv; ++j) {
        pLR = penLR(U[Ed[j]]-U[St[j]], Len[j], Pen[j], n);
        if (!ISNA(pLR)) {
          stat[i] = fmax(stat[i], pLR);
        }
      }
    }
    PutRNGstate();
    // store simulated results
    last_simulation = List::create(n, nsim, stat);
    Rprintf("... end!\n");
  }

  return stat;
}



