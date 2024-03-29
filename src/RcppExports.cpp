// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ncs
Rcpp::List ncs(IntegerMatrix dat, IntegerMatrix nminusy, IntegerVector z, int nspp, int nloc, int ngroup);
RcppExport SEXP _EcoCluster_ncs(SEXP datSEXP, SEXP nminusySEXP, SEXP zSEXP, SEXP nsppSEXP, SEXP nlocSEXP, SEXP ngroupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type nminusy(nminusySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type nspp(nsppSEXP);
    Rcpp::traits::input_parameter< int >::type nloc(nlocSEXP);
    Rcpp::traits::input_parameter< int >::type ngroup(ngroupSEXP);
    rcpp_result_gen = Rcpp::wrap(ncs(dat, nminusy, z, nspp, nloc, ngroup));
    return rcpp_result_gen;
END_RCPP
}
// rmultinom1
IntegerVector rmultinom1(NumericMatrix prob, NumericVector randu);
RcppExport SEXP _EcoCluster_rmultinom1(SEXP probSEXP, SEXP randuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type randu(randuSEXP);
    rcpp_result_gen = Rcpp::wrap(rmultinom1(prob, randu));
    return rcpp_result_gen;
END_RCPP
}
// getql
Rcpp::List getql(IntegerVector z, IntegerVector w, IntegerMatrix dat, int ngrloc, int ngrspp);
RcppExport SEXP _EcoCluster_getql(SEXP zSEXP, SEXP wSEXP, SEXP datSEXP, SEXP ngrlocSEXP, SEXP ngrsppSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type w(wSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type ngrloc(ngrlocSEXP);
    Rcpp::traits::input_parameter< int >::type ngrspp(ngrsppSEXP);
    rcpp_result_gen = Rcpp::wrap(getql(z, w, dat, ngrloc, ngrspp));
    return rcpp_result_gen;
END_RCPP
}
// convertSBtoNormal
NumericVector convertSBtoNormal(NumericVector v);
RcppExport SEXP _EcoCluster_convertSBtoNormal(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(convertSBtoNormal(v));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EcoCluster_ncs", (DL_FUNC) &_EcoCluster_ncs, 6},
    {"_EcoCluster_rmultinom1", (DL_FUNC) &_EcoCluster_rmultinom1, 2},
    {"_EcoCluster_getql", (DL_FUNC) &_EcoCluster_getql, 5},
    {"_EcoCluster_convertSBtoNormal", (DL_FUNC) &_EcoCluster_convertSBtoNormal, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_EcoCluster(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
