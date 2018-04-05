#ifndef IMPUTE2_H

// Find the columns with the SNP data
Rcpp::List FindImpute2SNPData(std::string &i2file, const int dataStart);

// Read the subject IDs from the header of an Impute 2 file
Rcpp::List ReadImpute2Subjects(std::string &i2file, const int dataStart, const bool usesFID);

// Count the number of subjects in an Impute 2 file
Rcpp::List CountImpute2Subjects(std::string &i2file, const int dataStart, const int format);

#endif
