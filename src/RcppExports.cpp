// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// msQuantile
NumericVector msQuantile(IntegerVector left, IntegerVector right, int n, int nsim, bool isGen);
RcppExport SEXP _essHist_msQuantile(SEXP leftSEXP, SEXP rightSEXP, SEXP nSEXP, SEXP nsimSEXP, SEXP isGenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type left(leftSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type right(rightSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< bool >::type isGen(isGenSEXP);
    rcpp_result_gen = Rcpp::wrap(msQuantile(left, right, n, nsim, isGen));
    return rcpp_result_gen;
END_RCPP
}
// boundedHistogram
DataFrame boundedHistogram(NumericVector orderedData, NumericVector cumCount, IntegerVector start, IntegerVector rightIndex, NumericVector lower, NumericVector upper);
RcppExport SEXP _essHist_boundedHistogram(SEXP orderedDataSEXP, SEXP cumCountSEXP, SEXP startSEXP, SEXP rightIndexSEXP, SEXP lowerSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type orderedData(orderedDataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cumCount(cumCountSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type rightIndex(rightIndexSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(boundedHistogram(orderedData, cumCount, start, rightIndex, lower, upper));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_essHist_msQuantile", (DL_FUNC) &_essHist_msQuantile, 5},
    {"_essHist_boundedHistogram", (DL_FUNC) &_essHist_boundedHistogram, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_essHist(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
