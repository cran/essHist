#include "stepHistogram.h"

/***************

 * class StepHistogram

 * Housen Li & Hannes Sieling, 2016

 ***************/

/*************

 * constructor for data points

 ****************/

//StepHistogram::StepHistogram(NumericVector xcs) : cs(xcs) {}


/*************

 * constructor for data points and bounds

 ****************/

StepHistogram::StepHistogram(NumericVector xcs, NumericVector xlb, NumericVector xub) : lb(xlb), ub(xub), cs(xcs) {}

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

double StepHistogram::costBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const {

  double nDataPoints  = endIndex - startIndex + 1;
  double startLoc; // = (startIndex == 0) ? cs[startIndex] : cs[startIndex] - cs[startIndex-1];
  double endLoc;   // = (endIndex == 0)   ? cs[endIndex]   : cs[endIndex] - cs[endIndex-1];

  unsigned int N = cs.size();

  // middle point as separation location

  if (startIndex == 0) {
    startLoc = (4*cs[0] - cs[1])/2;
  } else if (startIndex == 1) {
    startLoc = cs[1]/2;
  } else {
    startLoc = (cs[startIndex] - cs[startIndex-2])/2;
  }

  if (endIndex == 0) {
    endLoc = cs[1]/2;
  } else if (endIndex == N-1) {
    endLoc = (3*cs[N-1] - 4*cs[N-2] + cs[N-3])/2;
  } else {
    endLoc = (cs[endIndex+1] - cs[endIndex-1])/2;
  }

  double height = nDataPoints / fabs(endLoc - startLoc) / N;

  if (bound.lower > bound.upper) {
    return R_NaN;
  } else if (bound.lower > height || bound.upper < height) {
    return R_PosInf;
  } else {
    return -nDataPoints*log(height);
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

double StepHistogram::estBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const {

  double nDataPoints  = endIndex - startIndex + 1;
  double startLoc; // = (startIndex == 0) ? cs[startIndex] : cs[startIndex] - cs[startIndex-1];
  double endLoc;   // = (endIndex == 0)   ? cs[endIndex]   : cs[endIndex] - cs[endIndex-1];

  unsigned int N = cs.size();

  // middle point as separation location

  if (startIndex == 0) {
    startLoc = (4*cs[0] - cs[1])/2;
  } else if (startIndex == 1) {
    startLoc = cs[1]/2;
  } else {
    startLoc = (cs[startIndex] - cs[startIndex-2])/2;
  }

  if (endIndex == 0) {
    endLoc = cs[1]/2;
  } else if (endIndex == N-1) {
    endLoc = (3*cs[N-1] - 4*cs[N-2] + cs[N-3])/2;
  } else {
    endLoc = (cs[endIndex+1] - cs[endIndex-1])/2;
  }

  double height = nDataPoints / fabs(endLoc - startLoc) / N;

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

 * DataFrame : a list giving the solution

 ****************/

DataFrame StepHistogram::bounded(Bounds& B) const {

  // allocate storage

  unsigned int s = 0;                                        // number of jumps
  unsigned int N = cs.size();                                // number of data points
  unsigned int locR = N;                                     // 1st location we cannot reach in the current rightward search
  unsigned int locL = 0;                                     // 1st location we cannot reach in the last rightward search

  unsigned int auxlocR;                                      // auxilary storage for 'locR'

  int isstop;                                                // boolen variable for stoping criteria

  unsigned int auxKL;                                        // auxilary variable for computing KL

  unsigned int* const minJ = (unsigned int*) R_alloc(N, sizeof(unsigned int)); // the minimal number of jumps to reach k, indexed by k


  double* const J = (double*) R_alloc(N, sizeof(double)); // cost for optimal solution with s jumps over [0, ..., k], indexed by k
  int* const L = (int*) R_alloc(N, sizeof(int));          // last jump for optimal solution with s jumps over [0, ..., k], indexed by k
  double* const V = (double*) R_alloc(N, sizeof(double)); // estimate on last constant interval [l+1, ..., k], indexed by k
  double curJ = R_PosInf;                                 // the cost for current solution with last jump at l
  double curD = R_PosInf;                                 // cost for constant estimate on interval [l+1, ..., k]

  unsigned int* const K = (unsigned int*) R_alloc(N + 2, sizeof(unsigned int)); // maximal k s.t. s jumps are sufficient for feasible solution over [0, ... ,k], index s runs from -2 to n-1, requires offset of 2:
  unsigned int const Koffset = 2;
  int* const KL = (int*) R_alloc(N - 1, sizeof(int));     // minimal k s.t. s jumps are unnecessary for feasible solution over [0, ... ,k-1], index s runs from 1 to n-1, requires offset of -1:
  int const KLoffset = -1;
  unsigned int k;                                         // index for interval [0, ..., k]
  int l;                                                  // index of last jump, i.e. right index of second but last block


  // inititalize

  K[-2 + Koffset] = 0; K[-1 + Koffset] = 0; // for "negative number of jumps"

  for (k = 0; k < N; ++k) {
    J[k]    = R_PosInf;
    minJ[k] = N;
  }

  for(k = 0; k < N; k++) { // find constant solution over [0, ..., k]
    for(l = k - 1; l > 0; l--)  // precompute bounds on [l, k] for l > 0
      B.current(l, k);
    
    curD = costBound(0, k, B.current(0, k));

    if (R_FINITE(curD)) { // constant solution on []0, ..., k] exists
      curJ = curD;
      if(curJ < J[k]) { // improvement
        J[k] = curJ;
        L[k] = -1;
        V[k] = estBound(0, k, B.current(0, k));
      }
      K[0 + Koffset] = k;
      minJ[k] = 0;
    } else {
      if (k < locR)  // 1st location we cannot reach
        locR = k;

      if (ISNAN(curD))  // no constant solution on [0, ..., k] possible
        break;
    }
  }

  if (minJ[N - 1] == N) { // found no feasible solution on [0, ..., N-1]
    // calculations with at least one jump
    for (s = 1; s < N; s++) { // try to find solution with s jumps
      auxlocR = N;
      auxKL   = N;
      for (k = locR; k < N; ++k) {
        if (J[k] == R_PosInf) {
          B.reset();
          isstop = 1;
          for(l = k - 1; l >= (int) locL; --l) {
            for (int ll = l + 1; ll <= (int)k; ++ll)
              B.current(l + 1, ll);

            if (minJ[l] == s-1) {
              curD = costBound(l + 1, k, B.current(l + 1, k));
              if (ISNAN(curD)) { // no constant solution on [l+1, ..., k] possible
                break;
              } else {
                isstop = 0;
                if (R_FINITE(curD)) { // constant solution on [l+1, ..., k] exists
                  curJ = J[l] + curD;
                  if (curJ < J[k]) { // improvement
                    J[k] = curJ;
                    L[k] = l;
                    V[k] = estBound(l + 1, k, B.current(l + 1, k));
                  }
                  if (l < (int)auxKL) auxKL = l;
                }
              }
            }
          }

          if (R_FINITE(J[k])) {
            K[s + Koffset] = k;
            minJ[k] = s;
          } else {
            if (k < auxlocR)  // 1st location we cannot reach
              auxlocR = k;
            
            if (isstop == 1)  // no feasible solution with s jumps on [0, ..., k]
              break;
          }
        }
      }

      KL[s + KLoffset] = auxKL;
      locL = locR;
      locR = auxlocR;

      if (minJ[N - 1] < N) break; // found a feasible solution on [0, ..., N-1]
    }
  }

  // return result
  IntegerVector rightIndex(s + 1);
  NumericVector value(s + 1);
  IntegerVector jumpLeftBound(s + 1);
  IntegerVector jumpRightBound(s + 1);

  k = N - 1;
  for (int i = s; i >= 0; i--) {
    rightIndex[i] = k + 1;  // turn from C-style index into R-style index
    value[i] = V[k];
    if (i == (int) s) {
      jumpLeftBound[i] = N;
      jumpRightBound[i] = N;
    } else {
      jumpLeftBound[i] = KL[i] + 1;  // turn from C-style index into R-style index
      jumpRightBound[i] = K[i + Koffset] + 1;  // turn from C-style index into R-style index
    }
    k = L[k];
  }

  return DataFrame::create(Named("rightIndex") = rightIndex, Named("value") = value, Named("rightIndexLeftBound") = jumpLeftBound, Named("rightIndexRightBound") = jumpRightBound); //Named("cost") = J[N-1]
}



/*************

 * boundedHistogram

 * function to be called from R, wrapper for StepHistogram::bounded

 * computes bounded solution

 *

 * in:

 * cumSum : a numeric vector, the cumulative sums of the signal

 * start : for every possible left index where intervals with this left index start in the list of intervals (increasing in left indices), NA if none

 * rightIndex : right indices of the intervals in the list, increasing for each left index

 * lower : the lower bounds for the estimator at the respective interval

 * upper : the upper bounds for the estimator at the respective interval

 *

 * out:

 * a list conprimising the right indices and values of the solution, and the associated cost

 ****************/

// [[Rcpp::export(".boundedHistogram")]]

DataFrame boundedHistogram(NumericVector cumSum, IntegerVector start, IntegerVector rightIndex, NumericVector lower, NumericVector upper) {

  // initialise object

  StepHistogram data = StepHistogram(cumSum, lower, upper);



  // check lengths

  if(cumSum.size() <= 1) perror("there must be more than one block");

  if(start.size() != (int) cumSum.size()) perror("length of start must match cumSum's");

  if(lower.size() != upper.size()) perror("lower must have same length as upper");

  if(upper.size() != rightIndex.size()) perror("upper must have same length as rightIndex");



  Bounds B = Bounds(start, rightIndex, lower, upper);



  // run algorithm

  return data.bounded(B); // the optimal feasible solution using minimal number of jumps

}









