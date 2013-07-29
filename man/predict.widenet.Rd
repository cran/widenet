\name{predict.widenet}
\alias{predict.widenet}

\title{
Predict Method for \code{"widenet"} Objects
}

\description{
Similar to other predict methods, this functions predicts fitted values, logits,
coefficients and more from a fitted \code{"widenet"} object.
}

\usage{
\method{predict}{widenet}(object,
        newx,
        order = object$which.order.min,
        alpha.val = object$which.alpha.min,
        type = c("link", "response", "coefficients", "nonzero", "class"),
	...)
}

\arguments{
  \item{object}{
The \code{"widenet"} object from which to get predictions.
}
  \item{newx}{
Matrix of new values for \code{x} at which predictions are to be made. Must be a matrix; sparse matrices are not yet supported. This argument is not used for \code{type=c("coefficients","nonzero")}
}
  \item{order}{
Specifies the order for which predictions are required. Must be equal to one of the elements of the \code{order} argument for the \code{object}. The default is to use that order which "won" the cross-validation.
}
  \item{alpha.val}{
Value of \code{alpha} at which predictions are required. Default is to use that value which "won" the cross-validation.
}
  \item{type}{
See \code{link[glmnet]{predict.glmnet}}.
}
  \item{\dots}{
Further arguments passed to \code{predict.cv.relaxnet}, which should be passed through to \code{predict.glmnet}. However, in the current version, these are not guaranteed to work correctly (for example, \code{offset} has not yet been implemented for widenet).
}
}

\value{
The object returned depends on type.
}

\author{
Stephan Ritter, with design contributions from Alan Hubbard.

Much of the code (and some help file content) is adapted from the \pkg{glmnet} package, whose authors are Jerome Friedman, Trevor Hastie and Rob Tibshirani.
}

\seealso{
\code{\link{widenet}}, \code{\link{predict.cv.relaxnet}}, \code{\link{predict.glmnet}}
}
