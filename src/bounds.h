#include <Rcpp.h>

using namespace Rcpp;



/***************

 * class LUBound

 * specifies lower and upper bound

 * Housen Li & Hannes Sieling, 2016 (based on stepR package)

 ***************/



class LUBound {

public:

  double lower; // lower bound

  double upper; // upper bound



  // initialization

  LUBound();                     // default constructor setting bounds to -+Inf, i.e. always fulfilled

  LUBound(double lb, double ub); // constructor given lower and upper bound



  // add more constraints

  void add(double lb, double ub);

  void add(LUBound b);



  // check feasibility

  bool feasible();

};





/***************

 * class Bounds

 * stores lower and upper bounds over intervals

 * Housen Li & Hannes Sieling, 2016 (based on stepR package)

 ***************/



class Bounds {

private:

  // given bounds

  unsigned int N;   // number of data points

  int* const li;    // for each left index, where does it start, NA if no interval starting with it; has length N

  unsigned int Ni;  // number of intervals

  int* const ri;    // right indices; has length Ni

  double* const lb; // lower bound; has lenght Ni

  double* const ub; // upper bound; has length Ni

  // currently computed bounds, always just for a certain right index

  int* nexti;   // next index to look at for each left index; has length N

  int* cri;     // current right index for each left index; has length N

  LUBound* cb;  // current bound for each left index; has length N

  // store the initial state

  int* nexti0;  // next index to look at for each left index; has length N

  int* cri0;    // current right index for each left index; has length N

  LUBound* cb0; // current bound for each left index; has length N



public:

  // initialization

  // Bounds(unsigned int n, int* xstart, unsigned int ni, int* xri, double* xlb, double * xub); // constructor

  Bounds(IntegerVector xstart, IntegerVector xri, NumericVector xlb, NumericVector xub); // constructor



  LUBound current(unsigned int li, unsigned int ri); // return current bound for interval [li,ri], assuming it has been computed for [li,ri-1] and [li+1,ri] before

  void reset(); // set back to the initial state

};

