#include <Rcpp.h>
using namespace Rcpp;

//' Function to read an Impute 2 file
//'
//' Function to read an Impute 2 file
//'
//' @return
//' 0 success
//' 1 failure
//' @export
// [[Rcpp::export]]
List rcpp_hello() {
  CharacterVector x = CharacterVector::create("foo", "bar");
  NumericVector y   = NumericVector::create(0.0, 1.0);
  List z            = List::create(x, y);
  return z;
}

//' Function to convert and Impute 2 file to binary dosage format
//'
//' Function to convert and Impute 2 file to binary dosage format
//'
//' @param impute2File
//' Name of Impute 2 file
//' @param header
//' Does the impute 2 file have a header
//' @param snpInfoCols
//' Column numbers for chromosome, SNP name, location in BP, allele 1 and allele 2, respecitively
//' @param dataStartCol
//' Column where dosage/genetic probabilities start
//' @param binaryDosageFile
//' Name of output file
//' @return
//' String reporting status
//' @export
// [[Rcpp::export]]
std::string Impute2toBD_C(std::string &impute2File, bool header, Rcpp::IntegerVector &snpInfoCols, int dataStartCol, std::string &binaryDosageFile) {
  std::string status = "Successfully converted";

  if (header)
    Rcpp::Rcout << "Impute 2 file has a header" << std::endl;
  else
    Rcpp::Rcout << "Impute 2 file does not have a header" << std::endl;
  Rcpp::Rcout << snpInfoCols.size() << "\tcolumn numbers provided" << std::endl;
  Rcpp::Rcout << "Genetic data starts in column:\t" << dataStartCol << std::endl;
  return status;
}
