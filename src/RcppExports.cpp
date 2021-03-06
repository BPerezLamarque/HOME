// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// vmult
NumericVector vmult(NumericVector m1, NumericMatrix m2);
RcppExport SEXP _HOME_vmult(SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(vmult(m1, m2));
    return rcpp_result_gen;
END_RCPP
}
// vvmult
NumericVector vvmult(NumericVector m1, NumericVector m2);
RcppExport SEXP _HOME_vvmult(SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(vvmult(m1, m2));
    return rcpp_result_gen;
END_RCPP
}
// vvvmult
double vvvmult(NumericVector m1, NumericVector m2);
RcppExport SEXP _HOME_vvvmult(SEXP m1SEXP, SEXP m2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type m1(m1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m2(m2SEXP);
    rcpp_result_gen = Rcpp::wrap(vvvmult(m1, m2));
    return rcpp_result_gen;
END_RCPP
}
// llpruning
double llpruning(NumericVector nodes, NumericVector edges, NumericVector el, NumericMatrix L, NumericVector eig_val, NumericMatrix eig_vect, NumericMatrix ivp, NumericVector propinv);
RcppExport SEXP _HOME_llpruning(SEXP nodesSEXP, SEXP edgesSEXP, SEXP elSEXP, SEXP LSEXP, SEXP eig_valSEXP, SEXP eig_vectSEXP, SEXP ivpSEXP, SEXP propinvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type nodes(nodesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type edges(edgesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type el(elSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eig_val(eig_valSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type eig_vect(eig_vectSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ivp(ivpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type propinv(propinvSEXP);
    rcpp_result_gen = Rcpp::wrap(llpruning(nodes, edges, el, L, eig_val, eig_vect, ivp, propinv));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_HOME_vmult", (DL_FUNC) &_HOME_vmult, 2},
    {"_HOME_vvmult", (DL_FUNC) &_HOME_vvmult, 2},
    {"_HOME_vvvmult", (DL_FUNC) &_HOME_vvvmult, 2},
    {"_HOME_llpruning", (DL_FUNC) &_HOME_llpruning, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_HOME(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
