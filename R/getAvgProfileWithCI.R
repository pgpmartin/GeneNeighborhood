#' Compute the average profile and its confidence interval for a selection of features
#'
#' @param FeatureProfiles Matrix (features in row, position in columns) or RleList with coverage profiles for a set of features
#' @param selFeatures Character or integer vector defining a selection of features from \code{FeatureProfiles}
#' @param pos (Optional) Integer vector with values defining the genomic positions represented in \code{FeatureProfiles}
#' @param conf Numeric in ]0,1] representing the confidence level used to compute the confidence interval
#' @param type Character string indicating the method used to aggregate values from all the \code{selFeatures}. Methods are:
#'  \itemize{
#'    \item \code{"mean"} (Default). Compute the average value at each position.
#'    \item \code{"freq"} Compute the proportion of signal at each position. Sum of "freq" signal accross positions is 1.
#'    \item \code{"avgnorm"} Same as \code{"freq"} but mutiplied by the number of positions. Mean of "avgnorm" signal accross positions is 1.
#'    \item \code{"sum"} Compute the sum of signal at each position.
#'    \item \code{"n"} Compute the number of non missing observations at each genomic positions.
#'  }
#'
#' @return a data frame with 4 columns:
#'   \itemize{
#'   \item \code{Position}: pos values
#'   \item \code{Profile}: values obtained by aggregating as defined by \code{type} on all \code{selFeatures}
#'   \item \code{Upper}: Upper bound of the confidence interval of \code{Profile}
#'   \item \code{Lower}: Lower bound of the confidence interval of \code{Profile}
#'   }
#'
#' @export
#'
#' @importFrom S4Vectors elementNROWS
#' @importFrom stats var qnorm qt
#' @importFrom matrixStats colSds
#' @importFrom binom binom.confint
#'
#' @examples
#' ## Extract the profiles around (+/-50bp) the first 200 genes. We bin the genes in 3 bins only.
#'   top200Prof <- annotationCoverageAroundFeatures(Genegr,
#'                                                  features=1:200,
#'                                                  sidedist = 50,
#'                                                  usePercent = TRUE,
#'                                                  nbins=3)
#' ## Assemble the metagene profiles
#'   metaProf <- assembleMetagene(top200Prof)
#' ## Select a set of 50 random genes within the top200
#'   set.seed(123)
#'   randGenes <- sample(names(Genegr)[1:200], 50)
#' ##Get the average profile for these genes
#'   avgProf_sense <- getAvgProfileWithCI(metaProf$Metaprofiles_Sense,
#'                                        selFeatures = randGenes,
#'                                        pos = c(-50:0, 1:3, 0:50))
#'   head(avgProf_sense)
#' ## A quick look at the profile:
#'   plot(1:105, avgProf_sense$Profile,
#'        type="l", axes = FALSE,
#'        xlab="Position", ylab="Percent")
#'   lines(1:105, avgProf_sense$Upper, lty=2, col="grey")
#'   lines(1:105, avgProf_sense$Lower, lty=2, col="grey")
#'   axis(side=2)
#'   axis(side=1, at=c(1,51,55,105),
#'        labels=c("-50bp", "TSS", "TES", "+50bp"),
#'        las=2)

#'
getAvgProfileWithCI <- function(FeatureProfiles,
                                selFeatures = NULL,
                                pos = NULL,
                                conf = 0.95,
                                type = c("mean", "freq", "avgnorm", "sum", "n")) {


# Check arguments
  ## Check FeatureProfiles
  isRleList <- is(FeatureProfiles,"RleList")

  if (!isRleList && !is.matrix(FeatureProfiles))
    stop("FeatureProfiles should be an RleList or a matrix")

  if (isRleList && length(FeatureProfiles)==0) {
    stop("length of FeatureProfiles is 0")
  }

  if (isRleList &&
        stats::var(S4Vectors::elementNROWS(FeatureProfiles)) != 0) {
      stop("All elements of FeatureProfiles should have the same length")
    }

  ## Check type
  type <- match.arg(type)

  ## Check conf
  if (!is.numeric(conf) || (conf<=0 || conf > 1)) {
    stop("conf should be a numeric in ]0,1]")
  }

  ## Check selFeatures
  if (!missing(selFeatures) && length(selFeatures) == 1 && is.na(selFeatures)) {
    selFeatures <- NULL
  }

  if (!missing(selFeatures) && !is.null(selFeatures) &&
      !is.numeric(selFeatures) && !is.character(selFeatures)) {
    stop("selFeatures should be a numeric or a character vector")
  }

  if (!missing(selFeatures) && is.character(selFeatures) && any(is.na(selFeatures))) {
    stop("selFeatures should not contain missing values")
  }

  nfeatures <- if (isRleList) length(FeatureProfiles) else nrow(FeatureProfiles) # Number of features positions
  if (is.numeric(selFeatures)) {
    if (any(is.na(selFeatures))) {
      stop("selFeatures should not contain missing values")
    } else {
      if (min(selFeatures)<1 || max(selFeatures) > nfeatures) {
       stop("selFeatures is out of range")
      }
    }
  }

  ## Check pos and define it if NULL
  if (length(pos)>1 && any(is.null(pos))) {
    stop("pos contains NULL values")
  }

  if (length(pos)>1 && any(!is.finite(pos))) {
    stop("pos contains non finite values")
  }

  npos = if (isRleList) length(FeatureProfiles[[1]]) else ncol(FeatureProfiles) # Number of genomic positions
  if (is.null(pos))
    pos <- 1:npos
  if (length(pos)!=npos) {
    stop("pos length is not equal to the number of positions in FeatureProfiles")
  }

# Convert FeatureProfiles to a matrix if necessary
  if (isRleList) {
    FeatureProfiles <- RleList2matrix(FeatureProfiles)
  }

#Determine if FeaturesProfiles is binary (see https://stackoverflow.com/a/23276062)
  isBinary <- identical(as.vector(FeatureProfiles),
                        as.numeric(as.logical(FeatureProfiles)))

  ##Add rownames if NULL:
  featNames <- rownames(FeatureProfiles)
  if (is.null(featNames)) {
    rownames(FeatureProfiles) <-
      featNames <-
        paste0("Feature", 1:nrow(FeatureProfiles))
  }

# Data selection
  ## If selFeatures is NULL (or a single NA), select all features
  if (is.null(selFeatures)) {
    selFeatures <- 1:nrow(FeatureProfiles)
  }

  ## Remove from selFeatures the features with no profile
  if (is.character(selFeatures)) {
    isfeatName <- selFeatures %in% featNames

    if (all(!isfeatName)) {
    stop("none of the selfeatures are found in FeatureProfiles names")
    }

    if (any(!isfeatName)) {
      message(sum(!isfeatName), " selfeatures without profiles are excluded")
      selFeatures <- selFeatures[isfeatName]
    }
  }

  ## Select the Features
  FeatureProfiles <- FeatureProfiles[selFeatures,]


# Define the function applied to the columns
  ## mean by columns (default, type=="mean")
  if (type == "mean") {
    avgFUN <- function(data) {colMeans(data, na.rm = TRUE)} #faster than colSums(data)/nrow(data)
  }

  ## sum by columns (type=="sum")
  if (type == "sum") {
    avgFUN <- function(data) {colSums(data, na.rm = TRUE)}
  }

  ## number of (non NA) observations by columns (type=="n")
  if (type=="n") {
    avgFUN <- function(data) {colSums(!is.na(data))}
  }

  ## Frequence (type=="freq", proportion of signal in each column)
    ### Note that the values obtained depend on the grand sum, i.e. on the number of columns / window width.
    ### Do not use this to compare profiles of different width
  if (type == "freq"){
    avgFUN <- function(data) {
      ss <- colSums(data, na.rm = TRUE)
      ss <- ss / sum(ss, na.rm = TRUE)
    }
  }

  ## Normalized by the average (type=="avgnorm")
    ### Average signal by column divided by grand mean of the signal
    ### This is an analog to type=="freq" but more independent from the size of the window
    ### Use this (or the default type=="mean"...) rather than type="freq" to compare profiles of different width
  if (type=="avgnorm"){
    avgFUN <- function(data){
      ss <- colSums(data, na.rm = TRUE)
      ss <- ss * ncol(data) / sum(ss, na.rm = TRUE)
    }
  }

# Apply the function
  avgprof <- avgFUN(FeatureProfiles)

#Get confidence intervals.
  ## No confidence intervals for type=="sum" or "n"
  if (type %in% c("sum", "n")) {
    CI <- cbind(Upper = rep(NA, nrow(FeatureProfiles)),
                Lower = rep(NA, nrow(FeatureProfiles)))
  } else {

  ## Number of observations in each column (if n>30 we use a normal approximation)
    nobs <- colSums(!is.na(FeatureProfiles))
    NormalApprox <- min(nobs)>30

  ## When the matrix is NOT binary:
    if (!isBinary) {
    ### Get the SEM
      SEM <- matrixStats::colSds(FeatureProfiles, na.rm = TRUE) / sqrt(nobs)
      if (type == "freq") {
        SEM <- SEM * colSums(!is.na(FeatureProfiles)) / sum(FeatureProfiles, na.rm = TRUE)
      }
      if (type == "avgnorm") {
        SEM <- SEM * colSums(!is.na(FeatureProfiles)) * ncol(FeatureProfiles) / sum(FeatureProfiles, na.rm = TRUE)
      }
    ### Get the confidence intervals
      alpha <- 1 - conf
      if (NormalApprox) {
        CI <- cbind(Upper = avgprof + stats::qnorm(1-alpha/2) * SEM,
                    Lower = avgprof + stats::qnorm(alpha/2) * SEM)
      } else {
        CI <- cbind(Upper = avgprof + stats::qt(1-alpha/2, df = pmax(1, nobs-1)) * SEM,
                    Lower = avgprof + stats::qt(alpha/2, df = pmax(1, nobs-1)) * SEM)
      }
    } else {
    ## When the matrix is binary:
    ## For details, see https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
      if (NormalApprox) {
        SD <- sqrt(avgprof * (1-avgprof) / nobs)
        if (type == "freq") { ## unlikely choice if isBinary
          SD <- SD * colSums(!is.na(FeatureProfiles)) / sum(FeatureProfiles, na.rm = TRUE)
        }
        if (type == "avgnorm") { ## unlikely choice if isBinary
          SD <- SD * colSums(!is.na(FeatureProfiles)) * ncol(FeatureProfiles) / sum(FeatureProfiles, na.rm = TRUE)
        }
      ### Get the confidence intervals
        alpha <- 1 - conf
        CI <- cbind(Upper = avgprof + stats::qnorm(1-alpha/2) * SD,
                    Lower = avgprof + stats::qnorm(alpha/2) * SD)
      } else {
      ### When n<30 we use the conservative Clopper–Pearson interval (method="exact")
        isOne <- colSums(FeatureProfiles, na.rm = TRUE)
        CI <- binom::binom.confint(isOne,
                                   nobs,
                                   conf.level = conf,
                                   methods = "exact")[, c("upper", "lower")]
        colnames(CI) <- c("Upper", "Lower")
      }
    }
  }

return(data.frame(Position = pos,
                  Profile = avgprof,
                  CI))
}
