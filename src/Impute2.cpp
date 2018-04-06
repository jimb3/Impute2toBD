#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "Impute2.h"

//' Function to find the SNP data columns from an Impute 2 file
//'
//' Function to find the SNP data columns from an Impute 2 file.
//' Checks that the SNP name column exists and that there are no
//' duplicate columns.
//'
//' @param i2file
//' Name of impute 2 file
//' @param snpColNames
//' Vector of strings that are the names of the columns that contain
//' SNP name, chromosome, location in base pairs, reference allele,
//' and alternate allele, respectively. Unused values must be an empty
//' string.
//' @param dataStart
//' Column where the genetic data starts. Only the columns before
//' the genetic data are checked for SNP data
//' @return
//' List containing status and an
//' Integer vector of length 5 with column numbers for
//' SNP name, chromosome, location in BP, reference allele
//' and alternate allele. 0 is used when value passed in snpCol is empty.
//' @export
// [[Rcpp::export]]
Rcpp::List FindImpute2SNPData(std::string &i2file, Rcpp::StringVector &snpColNames, const int dataStart) {
  Rcpp::IntegerVector z(5);
  std::string status = "Good";
  Rcpp::List retval = Rcpp::List::create(Rcpp::Named("snpCol") = z, Rcpp::Named("status") = status);
  std::vector<std::string> cSNPColNames = Rcpp::as<std::vector<std::string> >(snpColNames);
  std::string firstLine;
  std::ifstream infile;
  std::istringstream instring;
  std::string colName;
  int i, j;

// Open the file
  infile.open(i2file.c_str());
  // Did file open successfully?
  if (!infile.good()) {
    retval["status"] = "Unable to open Impute 2 file";
    return retval;
  }
  // Read the first line
  std::getline(infile, firstLine);
  // Was there an error reading the first line?
  if (!infile.good()) {
    retval["status"] = "Error reading header";
    infile.close();
    return retval;
  }
  infile.close();

  instring.str(firstLine);
  // Find the columns with the SNP data
  // Only need to check the columns before where the data starts
  for (i = 0; i < dataStart; ++i) {
    instring >> colName;
    // Was there an error reading from the string? Most likely error is reached end of string.
    if (!instring.good()) {
      retval["status"] = "Error reading header";
      infile.close();
      return retval;
    }
    for (j = 0; j < 5; ++j) {
      if (cSNPColNames[j] == "")
        continue;
      if (colName == cSNPColNames[j]) {
        // Has the column already been found
        if (z[j] != 0) {
          retval["status"] = "More than one column found with column name " + snpColNames[j];
          return retval;
        }
        z[j] = i + 1;
        break;
      }
    }
  }
  for (i = 0; i < 5; ++i) {
    if (cSNPColNames[i] != "" && z[i] == 0) {
      retval["status"] = "Not all columns found";
      return retval;
    }
  }
  retval["snpCol"] = z;
  // Has a column been found named 'SNP'?
  if (z[0] == 0)
    retval["status"] = "Column for SNP name not found";

  return retval;
}

//' Function to read the subject IDs from an Impute 2 file
//'
//' Function to read the subject IDs from an Impute 2 file. This
//' requires that the file have a header
//'
//' @param i2file
//' Name of impute 2 file
//' @param dataStart
//' Column where the genetic data starts. This is also the column
//' where the subject IDs start.
//' @param usesFID
//' Indicator if family IDs are included in the header
//' @return
//' List containing status and a
//' data frame with the subject and family IDs. If no family ID is
//' in the file these values are set to blank, "", not missing.
//' @export
// [[Rcpp::export]]
Rcpp::List ReadImpute2Subjects(std::string &i2file, const int dataStart, const bool usesFID) {
  int numSub = 0;
  Rcpp::DataFrame subjects;
  Rcpp::List retval = Rcpp::List::create(Rcpp::Named("status") = "Good",
                                         Rcpp::Named("numSub") = numSub,
                                         Rcpp::Named("subjects") = subjects);
  std::ifstream infile;
  std::istringstream instring;
  std::string firstLine;
  std::string sub, fam;
  int i;

  // Open the file
  infile.open(i2file.c_str());
  // Did it open successfully?
  if (!infile.good()) {
    retval["status"] = "Unable to open Impute 2 file";
    return retval;
  }
  // Read first line
  std::getline(infile, firstLine);
  // Was it read successfully?
  if (!infile.good()) {
    retval["status"] = "Error reading header";
    infile.close();
    return retval;
  }
  infile.close();

  // Read the columns names before the subject IDs
  instring.str(firstLine);
  for (i = 1; i < dataStart; ++i) {
    instring >> sub;
  }
  if (!instring.good()) {
    retval["status"] = "Error reading column names before IDs";
    return retval;
  }

  instring >> sub;
  // Read until eof or an error occurs - an error can occur is there are blank spaces at the end of the line
  while(instring.good() && !instring.eof()) {
    ++numSub;
    instring >> sub;
  }
  // If the end of line was reached and no error occurred add the last item read
  // The error is probably blank spaces at the end of the line.
  if (instring.eof() && !instring.bad()) {
    ++numSub;
  }
  // If there are family IDs, the number of values read in must be even, FID and IID pairs
  if (usesFID) {
    if (numSub & 1) {
      retval["status"] = "Odd number of columns read, should be even if both subject and family ID are used";
      return retval;
    }
    numSub /= 2;
  }
  retval["numSub"] = numSub;
  // Were any subjects read in?
  if (numSub == 0) {
    retval["status"] = "No subjects found";
    return retval;
  }

  // Allocate the memory
  Rcpp::StringVector iid(numSub);
  Rcpp::StringVector fid(numSub);
  // Reset the string stream
  instring.clear();
  instring.str(firstLine);
  // Read the columns names before the subject IDs
  instring.str(firstLine);
  for (i = 1; i < dataStart; ++i) {
    instring >> sub;
  }
  for (i = 0; i < numSub; ++i) {
    if (usesFID) {
      instring >> sub;
      fid[i] = sub;
    }
    instring >> sub;
    iid[i] = sub;
  }
  subjects = Rcpp::DataFrame::create(Rcpp::Named("FID") = fid, Rcpp::Named("IID") = iid, Rcpp::Named("stringsAsFactors") = false);
  retval["subjects"] = subjects;

  return retval;
}

//' Function to count the subject in an Impute 2 file
//'
//' Function to count the subject in an Impute 2 file. This
//' requires that the file not have a header. No error checking
//' of the genetic values are done.
//'
//' @param i2file
//' Name of impute 2 file
//' @param dataStart
//' Column where the genetic data starts. This is also the column
//' where the subject IDs start.
//' @param format
//' Indicator if family IDs are included in the header
//' @return
//' List containing status and the number of subjects found
//' @export
// [[Rcpp::export]]
Rcpp::List CountImpute2Subjects(std::string &i2file, const int dataStart, const int format) {
  int numSub = 0;
  Rcpp::List retval = Rcpp::List::create(Rcpp::Named("status") = "Good",
                                         Rcpp::Named("numSub") = numSub);
  std::ifstream infile;
  std::istringstream instring;
  std::string firstLine;
  std::string junk;
  int i;

  // Open the file
  infile.open(i2file.c_str());
  // Did it open successfully?
  if (!infile.good()) {
    retval["status"] = "Unable to open Impute 2 file";
    return retval;
  }
  // Read first line
  std::getline(infile, firstLine);
  // Was it read successfully?
  if (!infile.good()) {
    retval["status"] = "Error reading first line";
    infile.close();
    return retval;
  }
  infile.close();

  // Read the columns names before the subject IDs
  instring.str(firstLine);
  for (i = 1; i < dataStart; ++i) {
    instring >> junk;
  }
  if (!instring.good()) {
    retval["status"] = "Error reading columns before genetic values";
    return retval;
  }

  instring >> junk;
  // Read until eof or an error occurs - an error can occur is there are blank spaces at the end of the line
  while(instring.good() && !instring.eof()) {
    ++numSub;
    instring >> junk;
  }
  // If the end of line was reached and no error occurred add the last item read
  // The error is probably blank spaces at the end of the line.
  if (instring.eof() && !instring.bad()) {
    ++numSub;
  }
  // Was the number of values read a multiple of the format?
  if ((numSub % format) != 0) {
    retval["status"] = "Number of values read in not a multiple of format";
    return retval;
  }
  retval["numSub"] = (int)(numSub / format);

  // Were any subjects read in?
  if (numSub == 0) {
    retval["status"] = "No subjects found";
    return retval;
  }

  return retval;
}
