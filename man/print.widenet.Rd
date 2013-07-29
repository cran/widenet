\name{print.widenet}
\alias{print.widenet}

\title{
Print Method for widenet Objects
}

\description{
This function just calls \code{print(summary(x))}. See \code{\link{summary.widenet}}.
}

\usage{
\method{print}{widenet}(x, digits, ...)
}

\arguments{
  \item{x}{
The \code{"widenet"} object to be printed.
}
  \item{digits}{
Passed to \code{print.summary.widenet}.
}
  \item{\dots}{
Passed to \code{print.summary.widenet}.
}
}
\value{
Returns \code{x} invisibly.
}

\author{
Stephan Ritter, with design contributions from Alan Hubbard.

Much of the code (and some help file content) is adapted from the \pkg{glmnet} package, whose authors are Jerome Friedman, Trevor Hastie and Rob Tibshirani.
}

\seealso{
\code{\link{widenet}}, \code{\link{summary.widenet}}
}
