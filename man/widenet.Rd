\name{widenet}
\alias{widenet}

\title{
Extends the relaxnet Package with Polynomial Basis Expansions
}

\description{
Expands the basis according to the \code{order} argument, then runs relaxnet in order to select a subset of the basis functions. Multiple values of \code{order} and \code{alpha} (the elastic net tuning parameter) may be specified, leading to selection of a specific value by cross-validation.
}

\usage{
widenet(x, y, family = c("gaussian", "binomial"),
        order = 1:3,
        alpha = 1,
        nfolds = 10,
        foldid,
        screen.method = c("none", "cor", "ttest"),
        screen.num.vars = 50,
        multicore = FALSE,
        mc.cores,
        mc.seed = 123,
        ...)
}

\arguments{
  \item{x}{
Input matrix, each row is an observation vector. Sparse matrices are not yet supported for the \code{widenet} function. Must have unique colnames.
}
  \item{y}{
Response variable. Quantitative for \code{family="gaussian"}. For \code{family="binomial"} should be either a factor with two levels, or a two-column matrix of counts or proportions.
}
  \item{family}{
Response type (see above).
}
  \item{order}{
The order of basis expansion. Elements must be in the set \code{c(1, 2, 3)}. If there is more than one element, cross-validation is used to chose the order with best cross-validated performance.
}
  \item{alpha}{
The elastic net mixing parameter, see \code{\link{glmnet}}. If there is more than one element, cross-validation is used to chose the value with best cross-validated performance.
}
  \item{nfolds}{
Number of folds - default is 10. Although \code{nfolds} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is \code{nfolds=3}.
}
  \item{foldid}{
An optional vector of values between 1 and \code{nfold} identifying what fold each observation is in. If supplied, \code{nfolds} can be missing.
}
  \item{screen.method}{
The method to use to screen variables before basis expansion is applied. Default is no screening. \code{"cor"} = correlation, i.e. bivariate correlation with the outcome. ttest is meant for binary outcomes (\code{family = "binomial"}). The screening methods are adapted from the \pkg{SuperLearner} package, the author of which is Eric Polley.
}
  \item{screen.num.vars}{
The number of variables (columns of \code{x} to screen in when using screening.
}
  \item{multicore}{
Should execution be parallelized over cv folds (for \code{cv.relaxnet}) or over alpha values (for \code{cv.alpha.relaxnet}) using multicore functionality from R's parallel package?
}
  \item{mc.cores}{
Number of cores/cpus to be used for multicore processing. Parallelization is over cross-validation folds.
}
  \item{mc.seed}{
Integer value with which to seed the RNG when using parallel processing (internally, \code{\link[base]{RNGkind}} will be called to set the RNG to \code{"L'Ecuyer-CMRG"}). Will be ignored if \code{multicore} is \code{FALSE}. If \code{mulicore} is \code{FALSE}, one should be able to get reprodicible results by setting the seed normally (with \code{\link[base]{set.seed}}) prior to running.
}
  \item{\dots}{
Further arguments passed to \code{\link{relaxnet}} or \code{\link{cv.relaxnet}}, which should also be passed on to \code{\link{glmnet}}. Use with caution as this has not been tested.
}
}


\details{
The \code{type.measure} argument has not yet been implemented. For type = gaussian models, mean squared error is used, and for type = binomial, binomial deviance is used.}


\value{
Returns and object of class \code{"widenet"} with the following elements:

  \item{call}{
A copy of the call which generated this object
}
  \item{order}{
The value of the \code{order} argument
}
  \item{alpha}{
The value of the \code{alpha} argument
}
  \item{screen.method}{
The value of the \code{screen.method} argument
}
  \item{screened.in.index}{
A vector which indexes the columns of \code{x}, indicating those variables which were screened in for the run on the full data
}
  \item{colsBinary}{
A vector of length \code{ncol(x)} representing which of the columns of \code{x} contained binary data. These columns will be represented by a 2. The other columns will have a 3.
}
  \item{cv.relaxnet.results}{
A list of lists containing \code{"cv.relaxnet"} objects, one for each combination of values of alpha and order.
}
  \item{min.cvm.mat}{
A matrix containing the minimum cross-validated risk for each combination of values of alpha and order
}
  \item{which.order.min}{
The order which "won" the cross-validation, i.e. resulted in minimum cross-validated risk.
}
  \item{which.alpha.min}{
The alpha value which "won" the cross-validation.
}
  \item{total.time}{
Total time in seconds to produce this result.
}
}

\references{
Stephan Ritter and Alan Hubbard, Tech report (forthcoming).
}


\author{
Stephan Ritter, with design contributions from Alan Hubbard.

Much of the code (and some help file content) is adapted from the \pkg{glmnet} package, whose authors are Jerome Friedman, Trevor Hastie and Rob Tibshirani.
}

\seealso{
\code{\link{predict.widenet}}, \code{\link{relaxnet}}, \code{\link{cv.relaxnet}}
}

\note{
This is a preliminary release and several additional features are planned for later versions.
}

\examples{

n <- 300
p <- 5

set.seed(23)
x <- matrix(rnorm(n*p), n, p)

colnames(x) <- paste("x", 1:ncol(x), sep = "")

y <- x[, 1] + x[, 2] + x[, 3] * x[, 4] + x[, 5]^2 + rnorm(n)

widenet.result <- widenet(x, y, family = "gaussian",
                          order = 2, alpha = 0.5)

summary(widenet.result)
coefs <- drop(predict(widenet.result, type = "coef"))
coefs[coefs != 0]
}