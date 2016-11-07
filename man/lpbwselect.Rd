\name{lpbwselect}
\alias{lpbwselect}
\alias{print.lpbwselect}
\alias{summary.lpbwselect}
\title{Bandwidth Selection Procedures for Local Polynomial Regression Estimators}

\description{
\code{\link{lpbwselect}} implements bandwidth selectors for local polynomial regression point estimators and inference procedures developed in  Calonico, Cattaneo and Farrell (2016a).

Companion commands are: \code{\link{lprobust}} for point estimation and inference procedures.

For more details, and related Stata and R packages useful for empirical analysis, visit \url{https://sites.google.com/site/rdpackages/}
}
\usage{
lpbwselect(y, x, c, p = 1, q = 2, 
deriv = 0, kernel = "epa", bwselect = "mse", scaleregul = 1, 
vce = "nn", nnmatch = 3, all = FALSE, subset = NULL)
}
\arguments{
  \item{y}{is the dependent variable.}
  \item{x}{is the independent variable).}
    \item{c}{specifies the evaluation point in \code{x}.}
  \item{p}{specifies the order of the local-polynomial used to construct the point-estimator; default is \code{p = 1} (local linear regression).}
  \item{q}{specifies the order of the local-polynomial used to construct the bias-correction; default is \code{q = 2} (local quadratic regression).}
  \item{deriv}{specifies the order of the derivative of the regression functions to be estimated. Default is \code{deriv=0}.}
  \item{kernel}{is the kernel function used to construct the local-polynomial estimator(s). Options are \code{triangular} (default option), \code{epanechnikov} and \code{uniform}.}
\item{bwselect}{specifies the bandwidth selection procedure to be used. Options are:

\code{mse} MSE-optimal bandwidth selector for the local polynomial regression estimator.

\code{cer} CER-optimal bandwidth selector for the local polynomial regression estimator.

Note: MSE = Mean Square Error; CER = Coverage Error Rate.
Default is \code{bwselect=mse}. For details on implementation see Calonico, Cattaneo and Farrell (2016a), and the companion software articles.}
  \item{scaleregul}{specifies scaling factor for the regularization term added to the denominator of the bandwidth selectors. Setting \code{scaleregul = 0} removes the regularization term from the bandwidth selectors; default is \code{scaleregul = 1}.}
  \item{vce}{specifies the procedure used to compute the variance-covariance matrix estimator.
Options are:

\code{nn} for heteroskedasticity-robust nearest neighbor variance estimator with \code{nnmatch} the (minimum) number of neighbors to be used.

\code{hc0} for heteroskedasticity-robust plug-in residuals variance estimator without weights.

\code{hc1} for heteroskedasticity-robust plug-in residuals variance estimator with \code{hc1} weights.

\code{hc2} for heteroskedasticity-robust plug-in residuals variance estimator with \code{hc2} weights.

\code{hc3} for heteroskedasticity-robust plug-in residuals variance estimator with \code{hc3} weights.

Default is \code{vce=nn}.}
  
    \item{nnmatch}{to be combined with for \code{vce=nn} for heteroskedasticity-robust nearest neighbor variance estimator with \code{nnmatch} indicating the minimum number of neighbors to be used. Default is \code{nnmatch=3}}
  \item{all}{if specified, \code{\link{lpbwselect}} reports all available bandwidth selection procedures.}
  \item{subset}{an optional vector specifying a subset of observations to be used.}
}

\value{
  \item{bws}{matrix containing the estimated bandwidths for each selected procedure.}
  \item{bwselect}{bandwidth selection procedure employed.}
  \item{kernel}{kernel function used to construct the local-polynomial estimator(s).}
  \item{p}{order of the local-polynomial used to construct the point-estimator.}
  \item{q}{order of the local-polynomial used to construct the bias-correction estimator.}
}

\references{
Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2016. On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference. Working Paper.
\url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Farrell_2016_JASA.pdf}.

Calonico, S., M. D. Cattaneo, M. H. Farrell, and R. Titiunik. 2016a. Regression Discontinuity Designs using Covariates. Working Paper.
\url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Farrell-Titiunik_2016_wp.pdf}.

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
y<-5+3*x+rnorm(1000)
lpbwselect(y,x, c=0)
}
\seealso{
  \code{\link{lprobust}}
}
