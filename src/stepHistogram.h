#include "bounds.h"

/***************

* class StepHistogram

* Housen Li & Hannes Sieling, 2016

***************/

class StepHistogram {

public:

  NumericVector lb; // lower bound

  NumericVector ub; // upper bound

  NumericVector cs; // cumulative sums

  // initialization
    
//  StepHistogram(NumericVector xcs); // constructor for n data points, with pointers to arrays of that length

  StepHistogram(NumericVector xcs, NumericVector xlb, NumericVector xub); // constructor for n data points and bounds

  // computations

  double cost(unsigned int startIndex, unsigned int endIndex) const; // calculate cost of a block

  double costBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const; // calculate cost of a block given bounds

  double estBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const; // corresponding estimate

  DataFrame bounded(Bounds& B) const; // compute optimal solution with minimal jumps fulfilling bounds
};

