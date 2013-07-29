\name{summary.widenet}
\alias{summary.widenet}
\alias{print.summary.widenet}

\title{
Generate and print summaries of class \code{"widenet"} objects.
}

\description{
Print summaries of widenet objects.
}

\usage{
\method{summary}{widenet}(object, ...)
\method{print}{summary.widenet}(x, digits = 3, ...)
}

\arguments{
  \item{object}{
The \code{"widenet"} object to summarize.
}
  \item{x}{
An object of class \code{"summary.widenet"}.
}
  \item{digits}{
Passed to \code{print.default}.
}
  \item{\dots}{
For the summary method: ignored.
For the print method: passed to \code{print.default}.
}
}

\value{
For the summary method: an object of type \code{"summary.widenet"} containing a subset of the elements of \code{object}.

For the print method: returns \code{x} invisibly.
}

\author{
Stephan Ritter, with design contributions from Alan Hubbard.

Much of the code (and some help file content) is adapted from the \pkg{glmnet} package, whose authors are Jerome Friedman, Trevor Hastie and Rob Tibshirani.
}

\seealso{
\code{\link{widenet}}, \code{\link{print.widenet}}
}
