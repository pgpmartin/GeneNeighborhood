% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{matrix2RleList}
\alias{matrix2RleList}
\title{Convert an \code{matrix} to an \code{RleList}}
\usage{
matrix2RleList(mat = NULL)
}
\arguments{
\item{mat}{A \code{\link{matrix}} to convert as an
\code{\link[IRanges]{RleList}}}
}
\value{
an \code{\link{RleList}}
}
\description{
Convert an \code{matrix} to an \code{RleList}.
             if mat is not a matrix, a conversion will be attempted.
}
\examples{
matrix2RleList(matrix(c(rep(1:3, 3), rep(4:6, 3)), nrow=3))

}
\author{
Pascal GP Martin
}
