#' Check that the table contains the mandatory columns
#' and keep only rows (cohorts) with the specified ancestry \code{pop}.
#' If \code{pop == "ALL"}, all rows are kept.
#'
#' @param df is the loaded data frame with individual cohorts information
#' @param ancest is the ancestry analyzed
#' @param pheno is the studied phenotype
#' @param expo is the studied exposure
#'
#' @return The data frame with only rows corresponding to the studied population \code{ancest}
#' and colums corresponding to the studied phenotype \code{pheno} and exposure \code{expo}
#'
#' @examples
#' data("COHORT")
#' df <- preparePhenoTable(COHORT, "EUR", "pheno1", "expo1")
#'
preparePhenoTable <- function(df, ancest, pheno, expo) {
  colsToSearch <- c("Cohort", "ANCESTRY", paste0(pheno, c("_N", "_Mean", "_SD")))
  missingCols <- colsToSearch[!(colsToSearch %in% colnames(df))]
  if (length(missingCols) > 0)
    stop(paste0("Missing column '", missingCols[i], "' in the cohort dataframe"), 
         call. = FALSE)
  
  if (ancest != "ALL") {
    colsToGrep <- c(colsToSearch, paste(pheno, expo, sep = "_"), 
                    paste0(expo, c("_Mean", "_SD")))
    df <- df[df$ANCESTRY == ancest, 
             grepl(paste(colsToGrep, collapse = "|"), colnames(df))]
    if (nrow(df) == 0) {
      stop(paste0("Cannot find ancestry ", ancest, " in the cohort dataframe"), 
           call. = FALSE)
    }
  } else {
    colsToGrep <- c(colsToSearch, paste(pheno, expo, sep = "_"))
    df <- df[, grepl(paste(colsToGrep,  collapse = "|"), colnames(df))]
  }
  
  stats::na.omit(df)
}

#' Check that the table contains the mandatory columns
#' and keep only rows (cohorts) with the specified ancestry \code{pop},
#' the specified phenoype \code{pheno}
#' and the specified exposure \code{expo}
#'
#' @param df is the loaded data frame with individual cohorts information
#' @param ancest is the ancestry analyzed
#' @param pheno is the studied phenotype
#' @param expo is the studied exposure
#'
#' @return The data frame with only rows corresponding to the studied population 
#' \code{ancest} and to the studied phenotype \code{pheno} and exposure \code{expo}
#'
#' @examples
#' data("GWAS")
#' df <- prepareResultTable(GWAS, "EUR", "pheno1", "expo1")
#'
prepareResultTable <- function(df, ancest, pheno, expo) {
  df <- subset(df, POP == ancest)
  if (nrow(df) == 0) {
    stop(paste0("Cannot find ancestry ", ancest, " in in the results dataframe"), 
         call. = FALSE)
  }
  df <- subset(df, PHENO == pheno & EXPO == expo)
  if (nrow(df) == 0) {
    stop(paste("Cannot find phenotype - exposure", paste(pheno, expo), 
               "in the results dataframe"), call. = FALSE)
  }
  df
}

#' Calculate the mean and variance of a quantitative variable in a pooled sample 
#' of several cohorts
#'
#' @param N is a vector of sample size in each cohorts
#' @param m is a vector of the mean of the variable in each individual cohort
#' @param v is a vector of the standard deviation of the variable in each 
#'   individual cohorts
#'
#' @return A vector of length 2 which first element is the mean and 
#'   second element is the variance
#'
#' @examples
#' sample_sizes <- c(250, 1000, 10000, 7500)
#' means <- c(2.5, 2.28, 2.32, 2.42)
#' sds <- c(1.05, 1.1, 0.98, 0.94)
#' parameters <- calculateContParams(sample_sizes, means, sds)
#'
#' @export
#'
calculateContParams <- function(N, m, v) {
  gmean <- sum(N * m) / sum(N)
  aa <- (N - 1) * v^2
  bb <- N * (m - gmean)^2
  c(gmean, sum(aa + bb) / (sum(N) - 1))
}

#' Calculate the mean and the variance of the exposure
#'
#' @param df is the dataframe with the cohort information
#' @param pheno is the studied outcome
#' @param expo is the studied exposure
#'
#' @return A vector of length 2 which first element is the mean and 
#'   second element is the variance
#'
#' @examples
#' # Case where E is quantitative
#' datafr <- data.frame(floor(rnorm(5, 5000, 2000)), runif(5, 2.5, 3), runif(5, 1.2, 1.4))
#' colnames(datafr) <- c("pheno_N", "pheno_expo_Mean", "pheno_expo_SD")
#' params <- calculateExpoParams(df = datafr, pheno = "pheno", expo = "expo")
#' # Case where E is binary
#' datafr <- data.frame(floor(rnorm(5, 5000, 2000)), floor(runif(5, 1000, 3000)))
#' colnames(datafr) <- c("pheno_N", "pheno_expo_P")
#' params <- calculateExpoParams(df = datafr, pheno = "pheno", expo = "expo")
#'
#' @export
#' 
calculateExpoParams <- function(df, pheno, expo) {
  if (any(grepl(paste0(pheno, "_", expo, "_P"), colnames(df)))) {
    n <- df[, grepl("_N", colnames(df))]
    nexp <- df[, grepl(paste0(expo, "_P"), colnames(df))]
    meanval <- sum(nexp) / sum(n)
    return(c(meanval, meanval * (1 - meanval)))
  }
  else if (any(grepl(paste0(pheno, "_", expo, "_Mean"), colnames(df))) &
           any(grepl(paste0(pheno, "_", expo, "_SD"), colnames(df)))) {
    n <- df[, grepl("_N", colnames(df))]
    m <- df[, grepl(paste0(expo, "_Mean"), colnames(df))]
    v <- df[, grepl(paste0(expo, "_SD"), colnames(df))]
    return(calculateContParams(n, m, v))
  }
  else {
    stop("Cannot calcuate exposure parameters.\nCheck columns names.")
  }
}

#' Derive betas in the standardized model from betas in the general model
#'
#' @param betaG is the vector of main genetic effects in the general model
#' @param betaINT is the vector of interaction effects in the general model
#' @param maf is the vector of the variants' frequency
#' @param meanE is the mean of the exposure
#' @param varE is the variance of the exposure
#' @param type designates the coefficients to standardize:
#'  "G" for the main genetic effect and "I" for the interaction effects
#'
#' @return The vector of standardized effects
#'
#' @examples
#' betaGs <- rnorm(10, 0, 0.1)
#' betaIs <- rnorm(10, 0, 0.05)
#' mafs <- runif(10, 0.05, 0.95)
#' meanE <- runif(1, -2, 2)
#' varE <- runif(1, 0.5, 1.5)
#' std_betaG <- standardizeBeta(betaGs, betaIs, mafs, meanE, varE, "G")
#' std_betaI <- standardizeBeta(betaGs, betaIs, mafs, meanE, varE, "I")
#'
#' tryCatch(
#'   std_betaI <- standardizeBeta(betaGs, betaIs, mafs, meanE, varE, 
#'                                "anything different from G or I"),
#'   error = function(e) print(e))
#'
#' @export
#' 
standardizeBeta <- function(betaG,betaINT,maf,meanE,varE,type) {
  if (type == "I") {
    return(betaINT * sqrt(2 * maf * (1 - maf)) * sqrt(varE))
  }
  else if (type == "G") {
    return((betaG + betaINT * meanE) * sqrt(2 * maf * (1 - maf)))
  }
  else {
    stop("type must be either \"I\" or \"G\"", call. = FALSE)
  }
}

#' Calculate the fraction of phenotypic variance explained by as set of 
#' significant genetic effect and/or interaction effects.
#'
#' This version does not take into account potential biases in neither the 
#' estimation of effect sizes nor the estimation of the correlation matrix
#'
#' @param std_bG is the vector of standardized genetic effect sizes
#' @param std_bI is the vector of standardized interaction effect sizes
#' @param matcor is the genotype correlation matrix
#' @param varY is the phenotypic variance in the pooled sample
#' @param type indicates wether only genetic or interactions effects should be 
#'   considered or if both should be considered jointly.
#'
#' @return The fraction of phenotypic variance explained by user-specified effects.
#'
#' @examples
#' std_bG <- rnorm(5, 0, 0.01)
#' std_bI <- rnorm(5, 0, 0.001)
#' matcor <- cor(matrix(runif(5*5, 1, 5), nrow = 5))
#' varY <- 2.25
#' calculateVarFrac(std_bG, std_bI, matcor, varY, "G")
#' calculateVarFrac(std_bG, std_bI, matcor, varY, "I")
#' calculateVarFrac(std_bG, std_bI, matcor, varY, "J")
#'
#' @importFrom MASS ginv
#' 
#' @export
#' 
calculateVarFrac <- function(std_bG, std_bI, matcor, varY, type) {
  if (type == "G") {
    std_bI <- rep(0, length(std_bG))
  }
  else if (type == "I") {
    std_bG <- rep(0, length(std_bI))
  }
  else if (type != "J") {
    stop("type must be in c(\"G\", \"I\", \"J\"", call. = FALSE)
  }
  return((crossprod(t(crossprod(std_bG, ginv(matcor))), std_bG) +
            crossprod(t(crossprod(std_bI, ginv(matcor))), std_bI)) / varY)
}

#' Calculate the fraction of phenotypic variance explained by as set of 
#' significant genetic effect and/or interaction effects.
#'
#' This version takes into account potential rank deficiencies
#' in the correlation matrix and noise in the  effect size estimation.
#'
#' For more details, see Shi et al., Am. J. Hum. Genet., 2016
#'
#' @param std_bG is the vector of standardized genetic effect sizes
#' @param std_bI is the vector of standardized interaction effect sizes
#' @param matcor is the genotype correlation matrix
#' @param varY is the phenotypic variance in the pooled sample
#' @param N is the total sample size
#' @param type indicates wether only genetic or interactions effects should be 
#'   considered or if both should be considered jointly.
#'
#' @return The fraction of phenotypic variance explained by user-specified effects.
#'
#' @examples
#' std_bG <- rnorm(5, 0, 0.01)
#' std_bI <- rnorm(5, 0, 0.001)
#' matcor <- cor(matrix(runif(5*5, 1, 5), nrow = 5))
#' varY <- 2.25
#' N <- 100000
#' calculateVarFrac_v2(std_bG, std_bI, matcor, varY, N, "G")
#' calculateVarFrac_v2(std_bG, std_bI, matcor, varY, N, "I")
#' calculateVarFrac_v2(std_bG, std_bI, matcor, varY, N, "J")
#'
#' @importFrom MASS ginv
#' 
#' @export
#' 
calculateVarFrac_v2 <- function(std_bG, std_bI, matcor, varY, N, type) {
  if (type == "G") {
    std_bI <- rep(0, length(std_bG))
  }
  else if (type == "I") {
    std_bG <- rep(0, length(std_bI))
  }
  else if (type != "J") {
    stop("type must be in c(\"G\", \"I\", \"J\")", call. = FALSE)
  }
  q = qr(matcor)$rank
  return((N * (crossprod(t(crossprod(std_bG, ginv(matcor))), std_bG) + 
                 crossprod(t(crossprod(std_bI, ginv(matcor))), std_bI)) -
            q) / ((N - q) * varY))
}
