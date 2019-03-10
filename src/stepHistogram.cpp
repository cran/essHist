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
StepHistogram::StepHistogram(NumericVector xod, NumericVector xlb, NumericVector xub) : lb(xlb), ub(xub), od(xod) {}

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
double StepHistogram::costBound(int startIndex, int endIndex, const LUBound& bound) const { // default: (a,b] but for including the first observation, trick is to use startIndex = -1
    double nDataPoints  = endIndex - startIndex; // endIndex - startIndex + 1;
    double height = nDataPoints / std::fabs(od[endIndex] - od[std::max(startIndex,0)]) / od.size();
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
double StepHistogram::estBound(int startIndex, int endIndex, const LUBound& bound) const { // default: (a,b] but for including the first observation, trick is to use startIndex = -1
    double nDataPoints  = endIndex - startIndex; // endIndex - startIndex + 1;
    double height = nDataPoints / std::fabs(od[endIndex] - od[std::max(startIndex,0)]) / od.size();
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
    unsigned int s = 0;            // number of jumps
    unsigned int N = od.size();    // number of data points
    unsigned int Nstar = 0;        // number of different values in data (could be different from N in case of discrete data)
    unsigned int locR;             // 1st location we cannot reach in the current rightward search
    unsigned int locL;             // 1st location we cannot reach in the last rightward search
    unsigned int auxlocR;          // auxilary storage for 'locR'
    int isstop;                    // boolen variable for stoping criteria
//    unsigned int auxKL;            // auxilary variable for computing KL (needed for confidence interval of locations)
//    unsigned int* const stId = (unsigned int*) R_alloc(N, sizeof(unsigned int)); // stId[i] is the start index of i-th different values of data
    unsigned int* const edId = (unsigned int*) R_alloc(N, sizeof(unsigned int)); // edId[i] is the end index of i-th different values of data (i.e. {0,...N-1} in case of continuous data)
    unsigned int* const minJ = (unsigned int*) R_alloc(N, sizeof(unsigned int)); // the minimal number of jumps to reach k, indexed by k
    double* const J = (double*) R_alloc(N, sizeof(double)); // cost for optimal solution with s jumps over ( 0, ..., edId[k] ], indexed by k
    int* const L = (int*) R_alloc(N, sizeof(int));          // ( edId[L[k]], ..., edId[k] ] is the last segment of the optimal solution over ( 0, ..., edId[k] ], indexed by k
    double* const V = (double*) R_alloc(N, sizeof(double)); // estimate on last constant interval ( edId[l], ..., edId[k] ], indexed by k
    double curJ = R_PosInf;                                 // the cost for current solution with last jump at l
    double curD = R_PosInf;                                 // cost for constant estimate on interval ( edId[l], ..., edId[k] ]
    // needed for confidence interval
//    unsigned int* const K = (unsigned int*) R_alloc(N + 2, sizeof(unsigned int)); // maximal k s.t. s jumps are sufficient for feasible solution over (0, ... ,k], index s runs from -2 to n-1, requires offset of 2:
//    unsigned int const Koffset = 2;
//    int* const KL = (int*) R_alloc(N - 1, sizeof(int));     // minimal k s.t. s jumps are unnecessary for feasible solution over (0, ... ,k-1], index s runs from 1 to n-1, requires offset of -1:
//    int const KLoffset = -1;
    
    unsigned int k;                                         // index for interval ( 0, ..., edId[k] ]
    int l;                                                  // edId[l] index of last jump, i.e. right index of second but last block
    int numNewStep;                                         // number of new steps made in each search, serving as a criterion for the existence problem
    
    // inititalize
    //      for discrete data
    for (k = 0; k < N-1; ++k)
        if (od[k] != od[k+1])
            edId[Nstar++] = k;
    edId[Nstar++] = N-1;
    //    K[-2 + Koffset] = 0; K[-1 + Koffset] = 0; // for "negative number of jumps"
    for (k = 0; k < Nstar; ++k) {
        J[k]    = R_PosInf;
        minJ[k] = Nstar;
    }
    J[0]       = 0;
    minJ[0]    = -1;
    locR       = Nstar;
    locL       = 1;
    numNewStep = 0;
    
    // first search
    for (k = locL; k < Nstar; ++k) { // find constant solution over ( 0, ..., edId[k] ] and k >= locL = 1
        for (l = k - 1; l > 0; --l)  // precompute bounds on ( edId[l], ..., edId[k] ] for l > 0
            B.current(edId[l], edId[k]);
        curD = costBound(-1, edId[k], B.current(0, edId[k])); // '-1' instead of '0' is to include the first observation
        if (R_FINITE(curD)) { // constant solution on (0, ..., k] exists
            ++numNewStep;
            curJ = curD;
            if (curJ < J[k]) { // improvement
                J[k] = curJ;
                L[k] = 0;
                V[k] = estBound(-1, edId[k], B.current(0, edId[k])); // '-1' instead of '0' is to include the first observation
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
    if (0 == numNewStep)
        Rcpp::stop("Search No. 1: No solution exists!");
    
    // later searches if needed
    if (Nstar == minJ[Nstar-1]) { // found no feasible solution on (0, ..., N-1] note that edId[Nstar-1] = N-1
        // calculations with at least one jump
        for (s = 1; s < Nstar; ++s) { // try to find solution with s jumps
            numNewStep = 0;
            auxlocR    = Nstar;
//            auxKL      = N;
            for (k = locR; k < Nstar; ++k) {
                if (J[k] == R_PosInf) {
                    B.reset();
                    isstop = 1;
                    for(l = k - 1; l >= (int) locL; --l) {
                        for (int ll = l + 1; ll <= (int)k; ++ll)
                            B.current(edId[l], edId[ll]);
                        
                        if (minJ[l] == s-1) {
                            curD = costBound(edId[l], edId[k], B.current(edId[l], edId[k]));
                            if (ISNAN(curD)) { // no constant solution on (l, ..., k] possible
                                break;
                            } else {
                                isstop = 0;
                                if (R_FINITE(curD)) { // constant solution on (l, ..., k] exists
                                    ++numNewStep;
                                    curJ = J[l] + curD;
                                    if (curJ < J[k]) { // improvement
                                        J[k] = curJ;
                                        L[k] = l;
                                        V[k] = estBound(edId[l], edId[k], B.current(edId[l], edId[k]));
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
            
            if (minJ[Nstar - 1] < Nstar) break; // found a feasible solution on (0, ..., N-1]
            
            if (0 == numNewStep) {
                Rcpp::stop("Search No. %d: No solution exists!", s+1);
            }
        }
    }
    
    // return result
    IntegerVector rightIndex(s + 1);
    NumericVector value(s + 1);
//    IntegerVector jumpLeftBound(s + 1);
//    IntegerVector jumpRightBound(s + 1);
    
    k = Nstar - 1;
    for (int i = s; i >= 0; --i) {
        rightIndex[i] = edId[k] + 1;  // turn from C-style index into R-style index
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
DataFrame boundedHistogram(NumericVector orderedData, IntegerVector start, IntegerVector rightIndex, NumericVector lower, NumericVector upper) {
    // initialise object
    StepHistogram data = StepHistogram(orderedData, lower, upper);
    
    // check lengths
    if(orderedData.size() <= 1) Rcpp::stop("there must be more than one block");
    if(start.size() != (int) orderedData.size()) Rcpp::stop("length of start must match orderedData's");
    if(lower.size() != upper.size()) Rcpp::stop("lower must have same length as upper");
    if(upper.size() != rightIndex.size()) Rcpp::stop("upper must have same length as rightIndex");
    
    Bounds B = Bounds(start, rightIndex, lower, upper);
    
    // run algorithm
    return data.bounded(B); // the optimal feasible solution using minimal number of jumps
}









