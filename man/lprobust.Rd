\name{lprobust}
\alias{lprobust}
\alias{print.lprobust}
\alias{summary.lprobust}
\title{Local-Polynomial Estimation with Robust Confidence Intervals}
\description{
\code{\link{lprobust}} implements local polynomial point estimators with robust bias-corrected confidence intervals and inference procedures developed in Calonico, Cattaneo and Farrell (2016a). It also computes alternative estimation and inference procedures available in the literature.

For more details, and related Stata and R packages useful for empirical analysis, visit \url{https://sites.google.com/site/rdpackages/}
}
\usage{
lprobust(y, x, c, p = 1, q = 2, deriv = 0, 
h = NULL, b = NULL, rho = NULL, kernel = "epa", bwselect = "mse", scaleregul = 1, 
vce = "nn", nnmatch = 3, level = 95, all = FALSE, subset = NULL)
}
\arguments{
  \item{y}{is the dependent variable.}
  \item{x}{is the independent variable.}
    \item{c}{specifies the evalution point in \code{x}.}
  \item{p}{specifies the order of the local-polynomial used to construct the point-estimator; default is \code{p = 1} (local linear regression).}
  \item{q}{specifies the order of the local-polynomial used to construct the bias-correction; default is \code{q = 2} (local quadratic regression).}
  \item{deriv}{specifies the order of the derivative of the regression functions to be estimated. Default is \code{deriv=0}}
  \item{h}{specifies the main bandwidth used to construct the LPR point estimator. If not specified, bandwidth \code{h} is computed by the companion command \code{\link{lpbwselect}}.}
  \item{b}{specifies the bias bandwidth used to construct the bias-correction estimator. If not specified, bandwidth \code{b} is computed by the companion command \code{\link{lpbwselect}}.}
  \item{rho}{specifies the value of \code{rho}, so that the bias bandwidth \code{b} equals \code{h/rho}. Default is \code{rho = 1} if \code{h} is specified but \code{b} is not.}
    \item{kernel}{is the kernel function used to construct the local-polynomial estimator(s). Options are \code{triangular} (default option), \code{epanechnikov} and \code{uniform}.}
  
\item{bwselect}{specifies the bandwidth selection procedure to be used. By default it computes both \code{h} and \code{b}, unless \code{rho} is specified, in which case it only computes \code{h} and sets \code{b=h/rho}.}
Options are:

\code{mse} MSE-optimal bandwidth selector for the local polynomial regression estimator.

\code{cer} CER-optimal bandwidth selector for the local polynomial regression estimator.

Note: MSE = Mean Square Error; CER = Coverage Error Rate.
Default is \code{bwselect=mse}. For details on implementation see Calonico, Cattaneo and Titiunik (2014a), Calonico, Cattaneo and Farrell (2016a),
and Calonico, Cattaneo, Farrell and Titiunik (2016), and the companion software articles.

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
  \item{level}{sets the confidence level for confidence intervals; default is \code{level = 95}.}
  \item{all}{if specified, \code{\link{lprobust}} reports three different procedures:
  
  (i) conventional LPR estimates with conventional standard errors.
  
  (ii) bias-corrected estimates with conventional standard errors.
  
  (iii) bias-corrected estimates with robust standard errors.}
  
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  
}

\value{
  \item{c}{cutoff value.}
  \item{N}{sample size used.}
  \item{N}{overall sample size.}
  \item{p}{order of the polynomial used for estimation of the regression function.}
  \item{q}{order of the polynomial used for estimation of the bias of the regression function.}
  \item{h}{bandwidth used for estimation of the regression function.}
  \item{b}{bandwidth used for estimation of the bias of the regression function estimator.}
  \item{coef}{vector containing conventional and bias-corrected local-polynomial estimates.}
  \item{se}{vector containing conventional and robust standard errors of the local-polynomial estimates.}
  \item{pv}{vector containing the p-values associated with conventional, bias-corrected and robust local-polynomial estimates.}
  \item{ci}{matrix containing the confidence intervals associated with conventional, bias-corrected and robust local-polynomial estimates.}
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
lprobust(y,x, c=0)
}
\keyword{ LPR }
\keyword{ Robust Estimation }
\seealso{
  \code{\link{lpbwselect}}
}