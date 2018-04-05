#' Function to get info about Impute 2 file
#'
#' Function to get info about Impute 2 file
#'
#' @param i2file
#' Impute 2 file name
#' @param header
#' Does the Impute 2 file have a header?
#' @param snpCol
#' Must be provided if there is no header.
#' Vector of column numbers for SNP name, chromosome, location in BP,
#' reference allele, and alternate allele respectively. All must be provided
#' a value of 0 indicates the value isn't in the file. SNP name column is required.
#' If provided and there is a header, column names in header are ignored.
#' If not provided, column names in header are used to find appropriate columns
#' In this case the column names must be 'SNP', 'CHR', 'BP', 'A1' and 'A2'.
#' @param startcol
#' Column number where dosages or genetic probabilities begin. Must be provided
#' @param format
#' Number of genetic probabilities in file, must be 1, 2 or 3
#' 1 indicates that only the dosage is in the file
#' @param usesFID
#' Does the impute file contain both the family and subject ID? Not used if
#' there is not header
#' @param sep
#' Separator used in file. Currently only space and tabs are allowed
#' @export
GetI2Info <- function(i2file, header = TRUE, snpCol, startCol, format = 2L, usesFID = TRUE, sep = '\t') {
  fileinfo <- list(status = "Good",
                   filename = "",
                   snpCol = rep(0L,5),
                   numSubjects = 0,
                   subjects = data.frame(FID = character(), IID = character(), stringsAsFactors = FALSE),
                   snps = data.frame(SNP = character(), CHR = character(), BP = integer(), A1 = character(), A2 = character(), stringsAsFactors = FALSE))
  # Is a file name provided?
  if (missing(i2file) == TRUE) {
    fileinfo$status <- "Missing file name"
#    stop("No file name specified")
    return (fileinfo)
  } else {
    # Is file name a string
    if (is.character(i2file) == FALSE) {
      fileinfo$status <- "File name is not a string"
      return (fileinfo)
    }
    fileinfo$filename <- i2file
  }
  # Is the header a logical value?
  if (is.logical(header) == FALSE) {
    fileinfo$status <- "Value for header is not a logical value"
    return (fileinfo)
  }
  # If there is no header, is there column numbers proviced?
  if (header == FALSE) {
    if (missing(snpCol) == TRUE) {
      fileinfo$status <- "No snp data columns specified"
      return (fileinfo)
    }
  }
  # Is a genetic data starting column number provided?
  if (missing(startCol) == TRUE) {
    fileinfo$status <- "No data start column specified"
    return (fileinfo)
  }
  # The starting column value must be an integer
  if ((startCol == as.integer(startCol)) == FALSE) {
    fileinfo$status <- "Starting data column must be an integer"
    return (fileinfo)
  }
  startCol <- as.integer(startCol)
  # Starting column number must be greater than 1 - 1 column is needed to identify the SNP
  if (startCol < 2) {
    fileinfo$status <- "Starting data column must be greater than 1"
    return (fileinfo)
  }
  # If column numbers are provided, is it an integer vector of length 5?
  if (missing(snpCol) == FALSE) {
    if (is.vector(snpCol) == FALSE) {
      fileinfo$status <- "SNP column number not a vector"
      return (fileinfo)
    }
    if (length(snpCol) != 5) {
      fileinfo$status <- "Length of SNP columns must be 5"
      return (fileinfo)
    }
    if (all(snpCol == as.integer(snpCol)) == FALSE) {
      fileinfo$status <- "SNP column number not an integer vector"
      return (fileinfo)
    }
    snpCol <- as.integer(snpCol)
    # Column numbers must be nonnegative
    if (min(snpCol) < 0) {
      fileinfo$status <- "SNP column numbers cannot be negative"
      return (fileinfo)
    }
    # SNP name column must exist, id est, > 0
    if (min(snpCol[1]) < 1) {
      fileinfo$status <- "SNP column for SNP name must be positive"
      return (fileinfo)
    }
    # SNP column numbers must be unique
    snpColNonZeros <- snpCol[snpCol > 0]
    snpUnique <- unique(snpColNonZeros)
    if (length(snpColNonZeros) != length(snpUnique)) {
      fileinfo$status <- "Nonzero SNP column numbers must be unique"
      return (fileinfo)
    }
    # Is the starting data column number greater than the largest
    # column number in  the snp column array?
    if (startCol <= max(snpCol)) {
      fileinfo$status <- "Data start column must be greater than largest SNP column value"
      return (fileinfo)
    }
  }
  # Is format an integer from 1 to 3?
  if ((as.integer(format) == format) == FALSE) {
    fileinfo$status <- "Value for format must be an integer from 1 to 3"
    return (fileinfo)
  }
  format <- as.integer(format)
  if (format < 1 | format > 3) {
    fileinfo$status <- "Value for format must be an integer from 1 to 3"
    return (fileinfo)
  }
  # Is usesFID a logical value?
  if (is.logical(usesFID) == FALSE) {
    fileinfo$status <- "Value for usesFID is not a logical value"
    return (fileinfo)
  }

  # Find the column numbers for the SNP data
  # This is only done if there is a header and snpCol was not provided
  if (missing(snpCol) == TRUE) {
    z <- FindImpute2SNPData(i2file, startCol)
    # Were the values successfully found, if not return the error
    if (z$status != "Good") {
      fileinfo$status <- z$status;
      return (fileinfo)
    }
    snpCol <- z$snpCol
  }
  fileinfo$snpCol = snpCol

  # Read subject IDs if there is a header otherwise count number of subjects
  if (header == TRUE) {
    z <- ReadImpute2Subjects(i2file, startCol, usesFID)
    if (z$status != "Good") {
      fileinfo$status <- z$status;
      return (fileinfo)
    }
    fileinfo$subjects <- z$subjects
    fileinfo$numSubjects <- z$numSub
  } else {
    z <- CountImpute2Subjects(i2file, startCol, format)
    if (z$status != "Good") {
      fileinfo$status <- z$status;
      return (fileinfo)
    }
    fileinfo$numSubjects <- z$numSub
  }

  # Read in the information about the SNPs
  numCol <- fileinfo$numSubjects * format + startCol - 1
  colClasses <- c(rep("NULL", numCol))
  snpColNonZeros <- snpCol[snpCol > 0]
  colClasses[snpColNonZeros] <- "character"
  if (snpCol[3] > 0)
    colClasses[snpCol[3]] <- "integer"
  snpColNames = c("SNP", "CHR", "BP", "A1", "A2")
  if (header == TRUE)
    snps <- read.table(i2file, sep = sep, skip = 1, colClasses = colClasses, quote = "")
  else
    snps <- read.table(i2file, sep = sep, header = FALSE, colClasses = colClasses, quote = "")
  colnames(snps) <- snpColNames[snpColNonZeros]
  fileinfo$snps <- snps

  return (fileinfo)
}
