#include "stepHistogram.h"

/***************
 * class StepHistogram
 * Housen Li, 2016-2019
 ***************/

/*************
 * constructor for data points
 ****************/

//StepHistogram::StepHistogram(NumericVector xcs) : cs(xcs) {}

/*************
 * constructor for data points and bounds
 ****************/
StepHistogram::StepHistogram(NumericVector xlb, NumericVector xub, NumericVector xod, NumericVector xcc) : lb(xlb), ub(xub), od(xod), cc(xcc) {}

/*************
 * costBound
 * calculate cost of a block given bounds
 *
 * in:
 * startIndex : the first index in the block, 0 <= startIndex <= endIndex
 * endIndex : the last index in the block, startIndex <= endIndex <= N - 1
 * bound : the lower and upper bound for the mean in this interval
 *
 * out:
 * the cost functional of that block
 ****************/
double StepHistogram::costBound(int startIndex, int endIndex, const LUBound& bound) const { // default: (a,b] 
  double nDataPoints = cc[endIndex] - cc[startIndex];
  double height = nDataPoints / std::fabs(od[endIndex] - od[startIndex]) / cc[cc.size()-1];
  
  if (bound.lower > bound.upper) {
    return R_NaN;
  } else if (bound.lower > height || bound.upper < height) {
    return R_PosInf;
  } else {
    return -nDataPoints*std::log(height);
  }
}

/*************
 * estBound
 * calculate estimate on a block given bounds
 *
 * in:
 * startIndex : the first index in the block, 0 <= startIndex <= endIndex
 * endIndex : the last index in the block, startIndex <= endIndex <= N - 1
 * bound : the lower and upper bound for the mean in this interval
 *
 * out:
 * the estimate for that block
 ****************/
double StepHistogram::estBound(int startIndex, int endIndex, const LUBound& bound) const { 
  // default: (a,b] but for including the first observation, trick is to use startIndex = -1
  // include the first observation
  double nDataPoints  = cc[endIndex] - cc[startIndex];
  double height = nDataPoints / std::fabs(od[endIndex] - od[startIndex]) / cc[cc.size()-1];
  if (bound.lower > bound.upper) {
    return R_NaN;
  } else if (bound.lower > height || bound.upper < height) {
    return R_NaN;
  } else {
    return height;
  }
}


/*************
 * bounded
 * compute optimal solution with minimal jumps fulfilling bounds
 *
 * out:
 * DataFrame : a list giving the solution ( left-continuous, to be compatible with intervals (a,b] )
 ****************/
DataFrame StepHistogram::bounded(Bounds& B) const {
  // allocate storage
  unsigned int s = 0;         // number of jumps
  unsigned int N = od.size(); // number of data points
  // unsigned int Nstar = 0;     // number of different values in data (could be different from N in case of discrete data)
  unsigned int locR;          // 1st location we cannot reach in the current rightward search
  unsigned int locL;          // 1st location we cannot reach in the last rightward search
  unsigned int auxlocR;       // auxilary storage for 'locR'
  int isstop;                 // boolen variable for stoping criteria
  //    unsigned int auxKL;            // auxilary variable for computing KL (needed for confidence interval of locations)
  //    unsigned int* const stId = (unsigned int*) R_alloc(N, sizeof(unsigned int)); // stId[i] is the start index of i-th different values of data
  // // edId[i] is the end index of i-th different values of data (i.e. {0,...N-1} in case of continuous data)
  // unsigned int* const edId = (unsigned int*) R_alloc(N, sizeof(unsigned int)); 
  unsigned int* const minJ = (unsigned int*) R_alloc(N, sizeof(unsigned int)); // the minimal number of jumps to reach k, indexed by k
  double* const J = (double*) R_alloc(N, sizeof(double)); // cost for optimal solution with s jumps over ( 0, ..., edId[k] ], indexed by k
  // ( edId[L[k]], ..., edId[k] ] is the last segment of the optimal solution over ( 0, ..., edId[k] ], indexed by k
  int* const L = (int*) R_alloc(N, sizeof(int));          
  double* const V = (double*) R_alloc(N, sizeof(double)); // estimate on last constant interval ( edId[l], ..., edId[k] ], indexed by k
  double curJ = R_PosInf;                                 // the cost for current solution with last jump at l
  double curD = R_PosInf;                                 // cost for constant estimate on interval ( edId[l], ..., edId[k] ]
  // needed for confidence interval
  //    unsigned int* const K = (unsigned int*) R_alloc(N + 2, sizeof(unsigned int)); // maximal k s.t. s jumps are sufficient for feasible solution over (0, ... ,k], index s runs from -2 to n-1, requires offset of 2:
  //    unsigned int const Koffset = 2;
  //    int* const KL = (int*) R_alloc(N - 1, sizeof(int));     // minimal k s.t. s jumps are unnecessary for feasible solution over (0, ... ,k-1], index s runs from 1 to n-1, requires offset of -1:
  //    int const KLoffset = -1;
  
  unsigned int k; // index for interval ( 0, ..., edId[k] ]
  int l;          // edId[l] index of last jump, i.e. right index of second but last block
  int noUpdate;   // whether there is an update in each search or not, serving as a criterion for the existence problem
  
  // inititalize
  //    K[-2 + Koffset] = 0; K[-1 + Koffset] = 0; // for "negative number of jumps"
  for (k = 0; k < N; ++k) {
    J[k]    = R_PosInf;
    minJ[k] = N;
  }
  // first search
  J[0]     = 0;
  minJ[0]  = -1;
  locR     = N;
  locL     = 1;
  noUpdate = 1;
  for (k = locL; k < N; ++k) { // find constant solution over ( 0, ..., k ] and k >= locL = 1
    for (l = k - 1; l > 0; --l)  // precompute bounds on ( l, ..., k ] for l > 0
      B.current(l, k);
    curD = costBound(0, k, B.current(0, k));
    if (R_FINITE(curD)) { // constant solution on (0, ..., k] exists
      noUpdate = 0;
      curJ = curD;
      if (curJ < J[k]) { // improvement
        J[k] = curJ;
        L[k] = 0;
        V[k] = estBound(0, k, B.current(0, k)); 
      }
      //            K[0 + Koffset] = k;
      minJ[k] = 0;
    } else {
      if (k < locR)  // 1st location we cannot reach
        locR = k;
      if (ISNAN(curD))  // no constant solution on (0, ..., k] possible
        break;
    }
  }
  if (1 == noUpdate)
    Rcpp::stop("Search No. 1: No solution exists!");
  
  // later searches if needed
  if (N == minJ[N-1]) { // found no feasible solution on (0, ..., N-1] note that edId[Nstar-1] = N-1
    // calculations with at least one jump
    for (s = 1; s < N; ++s) { // try to find solution with s jumps
      noUpdate = 1;
      auxlocR  = N;
      //            auxKL      = N;
      for (k = locR; k < N; ++k) {
        if (J[k] == R_PosInf) {
          B.reset();
          isstop = 1;
          for(l = k - 1; l >= (int) locL; --l) {
            for (int ll = l + 1; ll <= (int)k; ++ll)
              B.current(l, ll);
            if (minJ[l] == s-1) {
              curD = costBound(l, k, B.current(l, k));
              if (ISNAN(curD)) { // no constant solution on (l, ..., k] possible
                break;
              } else {
                isstop = 0;
                if (R_FINITE(curD)) { // constant solution on (l, ..., k] exists
                  noUpdate = 0;
                  curJ = J[l] + curD;
                  if (curJ < J[k]) { // improvement
                    J[k] = curJ;
                    L[k] = l;
                    V[k] = estBound(l, k, B.current(l, k));
                  }
                  //                                    if (l < (int)auxKL) auxKL = l;
                }
              }
            }
          }
          
          if (R_FINITE(J[k])) {
            //                        K[s + Koffset] = k;
            minJ[k] = s;
          } else {
            if (k < auxlocR)  // 1st location we cannot reach
              auxlocR = k;
            if (isstop == 1)  // no feasible solution with s jumps on (0, ..., k]
              break;
          }
        }
      }
      //            KL[s + KLoffset] = auxKL;
      locL = locR;
      locR = auxlocR;
      if (minJ[N - 1] < N) 
        break; // found a feasible solution on (0, ..., N-1]
      if (1 == noUpdate) 
        Rcpp::stop("Search No. %d: No solution exists!", s+1);
    }
  }
  
  // return result
  IntegerVector rightIndex(s + 1);
  NumericVector value(s + 1);
  //    IntegerVector jumpLeftBound(s + 1);
  //    IntegerVector jumpRightBound(s + 1);
  
  k = N - 1;
  for (int i = s; i >= 0; --i) {
    rightIndex[i] = k + 1;  // turn from C-style index into R-style index
    value[i]      = V[k];
    //        if (i == (int) s) {
    //            jumpLeftBound[i]  = N;
    //            jumpRightBound[i] = N;
    //        } else {
    //            jumpLeftBound[i]  = KL[i] + 1;  // turn from C-style index into R-style index
    //            jumpRightBound[i] = K[i + Koffset] + 1;  // turn from C-style index into R-style index
    //        }
    k = L[k];
  }
  
  return DataFrame::create(Named("rightIndex") = rightIndex, Named("value") = value);
  //Named("cost") = J[N-1], Named("rightIndexLeftBound") = jumpLeftBound, Named("rightIndexRightBound") = jumpRightBound
}



/*************
 * boundedHistogram
 * function to be called from R, wrapper for StepHistogram::bounded
 * computes bounded solution
 *
 * in:
 * orderedData : a numeric vector, the ordered data samples
 * start : for every possible left index where intervals with this left index start in the list of intervals (increasing in left indices), NA if none
 * rightIndex : right indices of the intervals in the list, increasing for each left index
 * lower : the lower bounds for the estimator at the respective interval
 * upper : the upper bounds for the estimator at the respective interval
 *
 * out:
 * a list conprimising the right indices and values of the solution, and the associated cost
 ****************/

// [[Rcpp::export(".boundedHistogram")]]
DataFrame boundedHistogram(NumericVector orderedData, NumericVector cumCount, IntegerVector start, IntegerVector rightIndex, NumericVector lower, NumericVector upper) {
  // initialise object
  StepHistogram data = StepHistogram(lower, upper, orderedData, cumCount);
  
  // check lengths
  if(orderedData.size() <= 1) 
    Rcpp::stop("there must be more than one block");
  if (cumCount.size() != orderedData.size())
    Rcpp::stop("lengths of 'cumCount' and 'orderedData' must match!");
  if(start.size() != (int) orderedData.size()) 
    Rcpp::stop("length of start must match orderedData's");
  if(lower.size() != upper.size()) 
    Rcpp::stop("lower must have same length as upper");
  if(upper.size() != rightIndex.size()) 
    Rcpp::stop("upper must have same length as rightIndex");
  
  Bounds B = Bounds(start, rightIndex, lower, upper);
  
  // run algorithm
  return data.bounded(B); // the optimal feasible solution using minimal number of jumps
}









