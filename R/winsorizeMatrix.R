#' @title Winsorizing of an entire matrix
#'
#' @description Winsorizing of an entire matrix. Here, winsorization consists in replacing extreme values in a matrix by the corresponding upper and lower quantiles.
#'
#' @param mat A matrix to winsorize
#' @param low Numeric. A value in [0,100]. All the values below \code{quantile(mat,low/100)} are replaced by \code{quantile(mat,low/100)}
#' @param high Numeric. A value in [0,100]. All the values above \code{quantile(mat,high/100)} are replaced by \code{quantile(mat,high/100)}
#' @param verbose Logical. Default to FALSE. If \code{verbose=TRUE}, the function will print statistics on the number of values replaced.
#'
#' @return A matrix with extreme values replaced by the corresponding quantiles.
#'
#' @importFrom stats quantile
#'
#' @export
#'
#' @examples
#' ## Create a matrix with random number from N(0,1):
#'   set.seed(123)
#'   x <- matrix(rnorm(1e4), ncol=100)
#' ## Replace the top and bottom 1% extreme values by the corresponding quantiles:
#'   y <- winsorizeMatrix(x, low=1, high=99, verbose = TRUE)
#' ## Add some NA in x:
#'   x[x>2] <- NA
#' ## winsorizeMatrix still works
#'   y <- winsorizeMatrix(x, low=1, high=99, verbose = TRUE)

winsorizeMatrix <- function(mat,
                            low = 1,
                            high = 99,
                            verbose = FALSE)
{
# Some controls
  if (length(low)!=1 || length(high)!=1) {
    stop("high and low must be single numeric values")
  }

  if (is.numeric(low) && is.numeric(high) && low>high) {
    stop("low must be less than high")
  }

  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("mat should be a matrix")
  }

# Convert mat if it is a data.frame
  if (is.data.frame(mat)) {
    mat <- data.matrix(mat)
    message("mat is a data frame and was converted to a matrix with data.matrix")
  }

# Initialize
  res <- mat

# Get lower quantile
  if (is.numeric(low)) {
    if (low<0 || low>100) {
      stop("lower threshold must be in [0:100]")
    } else {
      quantlow <- stats::quantile(mat, probs = low/100, na.rm = TRUE)
    }
  } else {
    stop("low should be a numeric")
  }

# Get upper quantile
  if (is.numeric(high)){
    if (high<0 | high>100){
      stop("higher threshold must be in [0:100]")
    } else {
      quanthigh <- stats::quantile(mat, probs = high/100, na.rm = TRUE)
    }
  } else {
    stop("high should be a numeric")
  }

# Replace extreme values by corresponding lower/upper quantile
  res[res<quantlow & !is.na(res)] <- quantlow
  res[res>quanthigh & !is.na(res)] <- quanthigh

# Verbose
  if (verbose) {
    ## Number of values affected by the changes:
      cat("Number of values replaced by the ", high, "th percentile: ",sum(mat>quanthigh & !is.na(mat)), "\n", sep="")
      cat("Number of values replaced by the ", low, "th percentile: ",sum(mat<quantlow & !is.na(mat)), "\n", sep="")
    ## Percentage of values affected by the changes in each column:
    diffmats <- abs(mat-res) > 1e-15
    percValuesChanged <- 100 * colMeans(diffmats, na.rm = TRUE)
    cat("Summary for the percentage of replaced values by columns:\n")
    print(summary(percValuesChanged))
  }

# Return modified matrix
  return(res)
}
