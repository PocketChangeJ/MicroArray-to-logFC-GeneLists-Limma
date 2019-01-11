### Limma and curation of resulting Gene Lists ###
#####
#
# Input Files
# a. Gene Expression Data, tab-separated (colnames, rownames and data)
# b. Classlabels, tab separated (0 or 1, last line empty)
# c. Platform, tab separated (2 columns; probe ids and gene symbols)
#
# Pre checking - do manually:
# a. Set respective working directory (setwd())
# b. boxplot and decide on setting toNormalize = TRUE or FALSE
#
### ~ ###

#####
### FUNCTIONS ###
#####

# reading input files into matrices
parseFile <- function(filename, h){
  matrix <- as.matrix(read.delim(filename, header = h))
  return(matrix)
}

checkLog <- function(intensities){
  toLog <- FALSE
  if (max(intensities) > 20) toLog <- TRUE
  return(toLog)
}

# handle data based on toLog, toNormalize or both, or un-log, normalize and re-log
handleData <- function(intensities, norm, log){
  if ( norm && log ){
    cat(sprintf("Normalizing and doing a log2 on the data.\n"))
    intensities <- normalizeQuantiles(intensities)
    intensities <- log2(intensities)
  } else if ( !norm && log ){
    cat(sprintf("Doing a log2 only on the data.\n"))
    intensities <- log2(intensities)
  } else if ( norm && !log ){ # need to un-log, normalize and then re-log
    cat(sprintf("Doing an Un-log, then normalizing and finally re-doing a log2 on the data.\n"))
    intensities <- 2 ^ intensities
    intensities <- normalizeQuantiles(intensities)
    intensities <- log2(intensities)
  }
  return(intensities)
}

runLimma <- function(intensities, classes){
  # transposing classes, as need for limma's model design
  classes <- t(classes)
  colnames(classes) <- "Group"
  design <- cbind(Intercept = 1, classes)
  fit <- lmFit(intensities, design)
  # Moderated t-statistic
  fit <- eBayes(fit)
  topTable <- topTable(fit, coef = 2, n = nrow(intensities))
  limma <- cbind(rownames(topTable), topTable$P.Value, topTable$logFC)
  colnames(limma) <- cbind("probeIDs","p-value","logFC")
  cat(sprintf("Limma finished.\n"))
  return(limma)
}

pvalueClear <- function(l, cut){
  # data type change
  l <- data.frame(l)
  l[, 1] <- as.matrix(as.character(l[, 1]))
  l[, 2] <- as.matrix(as.double(as.character(l[, 2])))
  l[, 3] <- as.matrix(as.double(as.character(l[, 3])))
  # actual cutoff
  l <- l[l[, 2] < cut, ] # delete where pvalue > 0.05
  cat(sprintf("pvalue entries > %f removed.\n", cut))
  return(l)
}

# if remaining Gene Symbols found both up- and down- regulated remove from list
removeDualSignProbes <- function(l){
  signMatrix <- matrix("" , nrow = 0, ncol = 2) # gene, FC value
  positions <- vector() # initialize empty vector of row positions which will be finally removed
  for (i in 1:nrow(l)){
    position <- match(l[i, 4], signMatrix[, 1])
    if (is.na(position)){ # if not checked yet add in signMatrix
      signMatrix <- rbind(signMatrix, c(as.character(l[i, 4]), l[i, 3]))
    } else if (sign(as.double(signMatrix[position, 2])) + sign(l[i, 3]) == 0){
      # remove all entries of this gene from l
      pos <- which(l[, 4] == l[i, 4])
      positions <- c(positions, pos)
      cat(sprintf("Dual sign Gene found and will be removed: %s\n", l[i, 4])) # debug
    }
  }
  # if is necessary else it deletes everything
  if (length(positions) > 0) l <- l[-as.integer(positions), ]
  return(l)
}

removeDuplicates <- function(l){
  # sort by pvalue
  l <- l[order(l[, 2], decreasing = FALSE), ]
  # remove duplicates on gene symbol
  l <- l[!duplicated(l[, 4]),]
  cat(sprintf("Duplicate Gene Symbol entries removed.\n"))
  return(l)
}

sortByAbsLogFC <- function(l){
  l <- l[order(abs(l[, 3]), decreasing = TRUE), ]
  cat(sprintf("Sorted by descending absolute logFC values.\n"))
  return(l)
}

#####
### MAIN ###
#####

# library
library(limma)

# directory
setwd("...")

# input variables
folder <- "experiment"
series_matrix_file <- paste(folder, "series_matrix.tsv", sep = "/")
classlabels_file <- paste(folder, "classlabels.tsv", sep = "/")
platform_file <- paste(folder, "platform.tsv", sep = "/")
# output filename
outfile <- paste(folder, "limma.tsv", sep = "/")

# parse files, filename + boolean header
series_matrix <- parseFile(series_matrix_file, TRUE) # character or double
classlabels <- parseFile(classlabels_file, FALSE) # integer
platform <- parseFile(platform_file, TRUE) # character (watch out for dates!)

# breaking series_matrix components to probe IDS and gene intensities
probeIDs <- matrix(as.character(series_matrix[, 1]))
colnames(probeIDs) <- "probeIDs"
temp <- as.matrix(series_matrix[,-1])
gene_intensities <- matrix(as.double(temp),
                       nrow = nrow(temp), ncol = ncol(temp))
colnames(gene_intensities) <- colnames(temp)
rownames(gene_intensities) <- probeIDs
remove(temp) # clearing memory

# boxplot and decide toNormalize = TRUE or FALSE
boxplot(gene_intensities)
toNormalize <- FALSE # TRUE

# check if intensities are log-ed
toLog <- checkLog(gene_intensities) # TRUE or FALSE

# handle input Data based on their values
gene_intensities <- handleData(gene_intensities, toNormalize, toLog)

# run Limma with gene intensities and classlabels
limma <- runLimma(gene_intensities, classlabels)

# clear p-value > 0.05
cutoff <- 0.05
limma <- pvalueClear(limma, cutoff)

# assigning gene symbols
colnames(platform)[1] <- "probeIDs"
limma <- merge(limma, platform)
# removing values without Gene Symbols first
limma <- limma[limma[, 4] != "", ] # --- or " " or ""

# experimental function
limma <- removeDualSignProbes(limma)

# removing lower pvalue duplicates
limma <- removeDuplicates(limma)

# sorting by descending absolute logFC
limma <- sortByAbsLogFC(limma)

# printing outfile
write.table(limma, outfile, quote = FALSE, row.names = FALSE, sep = "\t")

#####
### SCRIPT END ###
