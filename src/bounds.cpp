#include <cmath>
#include "bounds.h"

/***************
 * class LUBound
 * specifies lower and upper bound
 * Housen Li & Hannes Sieling, 2016 (based on stepR package on CRAN)
 ***************/

/*************
 * constructor given lower and upper bound
 ****************/
LUBound::LUBound(double lb, double ub) : lower(lb), upper(ub) {}

/*************
 * constructor with default bounds set to -+Inf, i.e. always fulfilled
 ****************/
LUBound::LUBound() : lower(R_NegInf), upper(R_PosInf) {}

/*************
 * add
 * add more constraints
 ****************/
void LUBound::add(double lb, double ub) {
  lower = lower>lb ? lower : lb;
  upper = upper<ub ? upper : ub;
}
void LUBound::add(LUBound b) {
  add(b.lower, b.upper);
}

/*************
 * feasible
 * check feasibility of constraints
 ****************/
bool LUBound::feasible() {
  return lower <= upper;
}

/***************
 * class Bounds
 * stores lower and upper bounds over intervals
 * Housen Li & Hannes Sieling, 2016 (on the basis of stepR package)
 ***************/

/*************
 * constructor
 ****************/
Bounds::Bounds(IntegerVector xstart, IntegerVector xri, NumericVector xlb, NumericVector xub): N(xstart.size()), li(xstart.begin()), Ni(xri.size()), ri(xri.begin()), lb(xlb.begin()), ub(xub.begin()) {
  // ensure there are bounds
  if (Ni < 1) Rf_error("no bounds specified!");
  
  // allocate arrays
  nexti = (int*) R_alloc(N, sizeof(int));
  cri   = (int*) R_alloc(N, sizeof(int));
  cb    = (LUBound*) R_alloc(N, sizeof(LUBound));
  
  // initialize ci, cri and cb, and check whether bounds can be fulfilled; corresponds to k = 0
  for(unsigned int i = 0; i < N; i++) {
    // initialize bound to be infinite
    cb[i] = LUBound();
    // current right index equals left index
    cri[i] = (int) i;
    // go through all intervals [i,i]
    for(nexti[i] = li[i]; nexti[i] != R_NaInt && nexti[i] < (int) Ni && ri[nexti[i]] == (int) i; nexti[i]++) {
      // check whether we've gone too far
      if(i < N - 1 && li[i+1] != R_NaInt && nexti[i] >= li[i+1]) {
        nexti[i] = R_NaInt; // indicate that there are no more constraints for this left index
        break;
      }
      cb[i].add(lb[nexti[i]], ub[nexti[i]]); // add constraint
    }
    // check whether we've gone too far
    if(nexti[i] >= (int) Ni) {
      nexti[i] = R_NaInt; // indicate that there are no more constraints for this left index
    }
    //     // check well-behaviour of indices
    //     if(nexti[i] != R_NaInt && nexti[i] >= Ni) error("index %d of interval with left index %d is >= number of intervals %d!", nexti[i], i, Ni);
    // check feasibility
    if(!cb[i].feasible()) Rf_error("Bounds not feasible at index %d!", i);
  }
  //************ Housen
  
  // store initial state
  //      allocate arrays
  nexti0 = (int*) R_alloc(N, sizeof(int));
  cri0   = (int*) R_alloc(N, sizeof(int));
  cb0    = (LUBound*) R_alloc(N, sizeof(LUBound));
  //      copy values
  for (unsigned int i = 0; i < N; ++i) {
    nexti0[i] = nexti[i];
    cri0[i]   = cri[i];
    cb0[i]    = LUBound(cb[i].lower, cb[i].upper);
  }
  
  //*******************
}

// some comments below is to ensure that it also works for search leftwards
/*************
 * current
 * return current bound for interval (l,r], assuming it has been computed for (l,r-1] and (l+1,r] before
 ****************/
LUBound Bounds::current(unsigned int l, unsigned int r) {
  // check whether these indices may be asked for
  if (r >= N  || r <= l) 
    Rprintf("indices must fulfill l %d < r %d < N %d", l, r, N);
  //    if ((int) r < cri[l]) error("for l %d we are already at cri %d, i.e. beyond r %d", l, cri[l], r);
  //    if ((int) r > cri[l] + 1) error("for l %d we are at cri %d, i.e. r %d is too far", l, cri[l], r);
  // check whether we have already computed this one
  if (cri[l] == (int) r) {
    return cb[l];
  }
  // add already computed bound on [li+1,ri]
  //    if (l < N - 1 && cri[l + 1] != (int) r) {
  //        error("bound for l + 1 = %d and r = %d needs to be available, but is at cri %d!", l + 1, r, cri[l + 1]);
  //    } else {
  cb[l].add(cb[l + 1]); // add constraint
  //    }
  // add all intervals (l,r] in the list
  for(; nexti[l] != R_NaInt && nexti[l] < (int) Ni && ri[nexti[l]] == (int) r; nexti[l]++) {
    // check whether we've gone too far
    if(l < N - 1 && li[l+1] != R_NaInt && nexti[l] >= li[l+1]) {
      nexti[l] = R_NaInt; // indicate that there are no more constraints for this left index
      break;
    }
    cb[l].add(lb[nexti[l]], ub[nexti[l]]); // add constraint
  }
  // update current right index
  cri[l] = r;
  // return bound on this interval
  return cb[l];
}


/*************
 * reset
 * set back to the initial state
 ****************/
void Bounds::reset() {
  for (unsigned int i = 0; i < N; ++i) {
    nexti[i]    = nexti0[i];
    cri[i]      = cri0[i];
    cb[i].lower = cb0[i].lower;
    cb[i].upper = cb0[i].upper;
  }
}
