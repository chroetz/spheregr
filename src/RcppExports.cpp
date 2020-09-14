// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// get_v_rcpp
NumericVector get_v_rcpp(NumericVector x, NumericVector a);
RcppExport SEXP _spheregr_get_v_rcpp(SEXP xSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(get_v_rcpp(x, a));
    return rcpp_result_gen;
END_RCPP
}
// angle2R3_1
NumericVector angle2R3_1(NumericVector x);
RcppExport SEXP _spheregr_angle2R3_1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(angle2R3_1(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spheregr_get_v_rcpp", (DL_FUNC) &_spheregr_get_v_rcpp, 2},
    {"_spheregr_angle2R3_1", (DL_FUNC) &_spheregr_angle2R3_1, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_spheregr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
