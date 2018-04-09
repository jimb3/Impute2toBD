// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// FindImpute2SNPData
Rcpp::List FindImpute2SNPData(std::string& i2file, Rcpp::StringVector& snpColNames, const int dataStart);
RcppExport SEXP _Impute2toBD_FindImpute2SNPData(SEXP i2fileSEXP, SEXP snpColNamesSEXP, SEXP dataStartSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type i2file(i2fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector& >::type snpColNames(snpColNamesSEXP);
    Rcpp::traits::input_parameter< const int >::type dataStart(dataStartSEXP);
    rcpp_result_gen = Rcpp::wrap(FindImpute2SNPData(i2file, snpColNames, dataStart));
    return rcpp_result_gen;
END_RCPP
}
// ReadImpute2Subjects
Rcpp::List ReadImpute2Subjects(std::string& i2file, const int dataStart, const bool usesFID);
RcppExport SEXP _Impute2toBD_ReadImpute2Subjects(SEXP i2fileSEXP, SEXP dataStartSEXP, SEXP usesFIDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type i2file(i2fileSEXP);
    Rcpp::traits::input_parameter< const int >::type dataStart(dataStartSEXP);
    Rcpp::traits::input_parameter< const bool >::type usesFID(usesFIDSEXP);
    rcpp_result_gen = Rcpp::wrap(ReadImpute2Subjects(i2file, dataStart, usesFID));
    return rcpp_result_gen;
END_RCPP
}
// CountImpute2Subjects
Rcpp::List CountImpute2Subjects(std::string& i2file, const int dataStart, const int format);
RcppExport SEXP _Impute2toBD_CountImpute2Subjects(SEXP i2fileSEXP, SEXP dataStartSEXP, SEXP formatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type i2file(i2fileSEXP);
    Rcpp::traits::input_parameter< const int >::type dataStart(dataStartSEXP);
    Rcpp::traits::input_parameter< const int >::type format(formatSEXP);
    rcpp_result_gen = Rcpp::wrap(CountImpute2Subjects(i2file, dataStart, format));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP _Impute2toBD_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}
// Impute2toBD_C
std::string Impute2toBD_C(std::string& impute2File, bool header, Rcpp::IntegerVector& snpInfoCols, int dataStartCol, std::string& binaryDosageFile);
RcppExport SEXP _Impute2toBD_Impute2toBD_C(SEXP impute2FileSEXP, SEXP headerSEXP, SEXP snpInfoColsSEXP, SEXP dataStartColSEXP, SEXP binaryDosageFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type impute2File(impute2FileSEXP);
    Rcpp::traits::input_parameter< bool >::type header(headerSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector& >::type snpInfoCols(snpInfoColsSEXP);
    Rcpp::traits::input_parameter< int >::type dataStartCol(dataStartColSEXP);
    Rcpp::traits::input_parameter< std::string& >::type binaryDosageFile(binaryDosageFileSEXP);
    rcpp_result_gen = Rcpp::wrap(Impute2toBD_C(impute2File, header, snpInfoCols, dataStartCol, binaryDosageFile));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Impute2toBD_FindImpute2SNPData", (DL_FUNC) &_Impute2toBD_FindImpute2SNPData, 3},
    {"_Impute2toBD_ReadImpute2Subjects", (DL_FUNC) &_Impute2toBD_ReadImpute2Subjects, 3},
    {"_Impute2toBD_CountImpute2Subjects", (DL_FUNC) &_Impute2toBD_CountImpute2Subjects, 3},
    {"_Impute2toBD_rcpp_hello", (DL_FUNC) &_Impute2toBD_rcpp_hello, 0},
    {"_Impute2toBD_Impute2toBD_C", (DL_FUNC) &_Impute2toBD_Impute2toBD_C, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_Impute2toBD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
