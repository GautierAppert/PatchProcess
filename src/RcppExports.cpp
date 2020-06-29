// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// computeLoop
void computeLoop(NumericVector& m, NumericVector& v, NumericVector& w, NumericVector& C, IntegerVector& A, IntegerVector& B, IntegerVector& Wsx, List& paramList);
RcppExport SEXP _PatchProcess34_computeLoop(SEXP mSEXP, SEXP vSEXP, SEXP wSEXP, SEXP CSEXP, SEXP ASEXP, SEXP BSEXP, SEXP WsxSEXP, SEXP paramListSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type C(CSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type A(ASEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type B(BSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Wsx(WsxSEXP);
    Rcpp::traits::input_parameter< List& >::type paramList(paramListSEXP);
    computeLoop(m, v, w, C, A, B, Wsx, paramList);
    return R_NilValue;
END_RCPP
}
// combineLabels
void combineLabels(IntegerVector& Swx, IntegerVector& Qwx, IntegerVector& Rwx, int S_height, int Q_height, int width);
RcppExport SEXP _PatchProcess34_combineLabels(SEXP SwxSEXP, SEXP QwxSEXP, SEXP RwxSEXP, SEXP S_heightSEXP, SEXP Q_heightSEXP, SEXP widthSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type Swx(SwxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Qwx(QwxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Rwx(RwxSEXP);
    Rcpp::traits::input_parameter< int >::type S_height(S_heightSEXP);
    Rcpp::traits::input_parameter< int >::type Q_height(Q_heightSEXP);
    Rcpp::traits::input_parameter< int >::type width(widthSEXP);
    combineLabels(Swx, Qwx, Rwx, S_height, Q_height, width);
    return R_NilValue;
END_RCPP
}
// computeMerge
void computeMerge(IntegerVector& A, IntegerVector& Asizes, IntegerVector& B, IntegerVector& Bsizes, IntegerVector& Twx, IntegerVector& Jt, List& paramList);
RcppExport SEXP _PatchProcess34_computeMerge(SEXP ASEXP, SEXP AsizesSEXP, SEXP BSEXP, SEXP BsizesSEXP, SEXP TwxSEXP, SEXP JtSEXP, SEXP paramListSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type A(ASEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Asizes(AsizesSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type B(BSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Bsizes(BsizesSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Twx(TwxSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Jt(JtSEXP);
    Rcpp::traits::input_parameter< List& >::type paramList(paramListSEXP);
    computeMerge(A, Asizes, B, Bsizes, Twx, Jt, paramList);
    return R_NilValue;
END_RCPP
}
// computeContext
void computeContext(IntegerVector& Jt_v, IntegerVector& Cont_v, IntegerVector& ContOrder_v, int number_of_labels, int number_of_merged_labels);
RcppExport SEXP _PatchProcess34_computeContext(SEXP Jt_vSEXP, SEXP Cont_vSEXP, SEXP ContOrder_vSEXP, SEXP number_of_labelsSEXP, SEXP number_of_merged_labelsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type Jt_v(Jt_vSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Cont_v(Cont_vSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type ContOrder_v(ContOrder_vSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_labels(number_of_labelsSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_merged_labels(number_of_merged_labelsSEXP);
    computeContext(Jt_v, Cont_v, ContOrder_v, number_of_labels, number_of_merged_labels);
    return R_NilValue;
END_RCPP
}
// applySyntax
void applySyntax(IntegerVector& Twx_v, IntegerVector& Gtt_v, IntegerVector& Jt_v, IntegerVector& Gt_v, IntegerVector& Uwx_v, IntegerVector& Axw_v, IntegerVector& Axu_v, int number_of_syntax_labels, int number_of_labels, int number_of_patches, int number_of_images);
RcppExport SEXP _PatchProcess34_applySyntax(SEXP Twx_vSEXP, SEXP Gtt_vSEXP, SEXP Jt_vSEXP, SEXP Gt_vSEXP, SEXP Uwx_vSEXP, SEXP Axw_vSEXP, SEXP Axu_vSEXP, SEXP number_of_syntax_labelsSEXP, SEXP number_of_labelsSEXP, SEXP number_of_patchesSEXP, SEXP number_of_imagesSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type Twx_v(Twx_vSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Gtt_v(Gtt_vSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Jt_v(Jt_vSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Gt_v(Gt_vSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Uwx_v(Uwx_vSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Axw_v(Axw_vSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type Axu_v(Axu_vSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_syntax_labels(number_of_syntax_labelsSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_labels(number_of_labelsSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_patches(number_of_patchesSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_images(number_of_imagesSEXP);
    applySyntax(Twx_v, Gtt_v, Jt_v, Gt_v, Uwx_v, Axw_v, Axu_v, number_of_syntax_labels, number_of_labels, number_of_patches, number_of_images);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PatchProcess34_computeLoop", (DL_FUNC) &_PatchProcess34_computeLoop, 8},
    {"_PatchProcess34_combineLabels", (DL_FUNC) &_PatchProcess34_combineLabels, 6},
    {"_PatchProcess34_computeMerge", (DL_FUNC) &_PatchProcess34_computeMerge, 7},
    {"_PatchProcess34_computeContext", (DL_FUNC) &_PatchProcess34_computeContext, 5},
    {"_PatchProcess34_applySyntax", (DL_FUNC) &_PatchProcess34_applySyntax, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_PatchProcess34(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
