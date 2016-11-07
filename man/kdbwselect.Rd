\name{kdbwselect}
\alias{kdbwselect}
\alias{print.kdbwselect}
\alias{summary.kdbwselect}
\title{Bandwidth Selection Procedures for Kernel Density Estimation}

\description{
\code{\link{kdbwselect}} implements bandwidth selectors for kernel density estimation and inference procedures developed in  Calonico, Cattaneo and Farrell (2016a).

Companion commands are: \code{\link{kdrobust}} for point estimation and inference procedures.

For more details, and related Stata and R packages useful for empirical analysis, visit \url{https://sites.google.com/site/rdpackages/}
}
\usage{
kdbwselect(x, c, deriv = 0, p = 2, kernel = "epan", 
bwselect = "mse", all = FALSE, subset = NULL)
}
\arguments{
  \item{x}{the data from which the estimate is to be computed.}
  \item{c}{specifies the evaluation point in \code{x}.}
  \item{p}{specifies the order of the kernel used to construct the point-estimator; default is \code{p = 2}.}
  \item{deriv}{specifies the order of the derivative of the density function to be estimated. Default is \code{deriv=0}.}
  \item{kernel}{is the kernel function used to construct the kernel density estimator(s). Options are \code{triangular} (default option), \code{epanechnikov} and \code{uniform}.}
\item{bwselect}{specifies the bandwidth selection procedure to be used. Options are:

\code{mse} MSE-optimal bandwidth selector for the kernel density estimator.

\code{cer} CER-optimal bandwidth selector for the kernel density estimator.

Note: MSE = Mean Square Error; CER = Coverage Error Rate.
Default is \code{bwselect=mse}. For details on implementation see Calonico, Cattaneo and Farrell (2016a), and the companion software articles.}

  \item{all}{if specified, \code{\link{kdbwselect}} reports all available bandwidth selection procedures.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
}

\value{
  \item{bws}{matrix containing the estimated bandwidths for each selected procedure.}
  \item{bwselect}{bandwidth selection procedure employed.}
  \item{kernel}{kernel function used to construct the kernel density estimator(s).}
  \item{p}{order of the kernel used to construct the point-estimator.}
  
}

\references{
Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2016. On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference. Working Paper.
\url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Farrell_2016_JASA.pdf}.

Calonico, S., M. D. Cattaneo, M. H. Farrell, and R. Titiunik. 2016b. rdrobust: Software for Regression Discontinuity Designs. Working Paper.
\url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Farrell-Titiunik_2016_Stata.pdf}.

Calonico, S., M. D. Cattaneo, and R. Titiunik. 2014a. Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs. Econometrica 82(6): 2295-2326.
\url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf}.

Calonico, S., M. D. Cattaneo, and R. Titiunik. 2014b. Robust Data-Driven Inference in the Regression-Discontinuity Design. Stata Journal 14(4): 909-946.
\url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_Stata.pdf}.

Calonico, S., M. D. Cattaneo, and R. Titiunik. 2015b. rdrobust: An R Package for Robust Nonparametric Inference in Regression-Discontinuity Designs. R Journal 7(1): 38-51.
\url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2015_R.pdf}.
}
\author{
Sebastian Calonico, University of Miami, Coral Gables, FL.  \email{scalonico@bus.miami.edu}.

Matias D. Cattaneo, University of Michigan, Ann Arbor, MI.  \email{cattaneo@umich.edu}.

Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.
}

\examples{
x<-runif(1000,-1,1)
kdbwselect(x, c=0)
}
\seealso{
  \code{\link{kdrobust}}
}
