\name{kdrobust}
\alias{kdrobust}
\alias{print.kdrobust}
\alias{summary.kdrobust}
\title{Kernel Density Estimation with Robust Confidence Intervals}
\description{
\code{\link{kdrobust}} implements kernel density estimators with robust bias-corrected confidence intervals and inference procedures developed in Calonico, Cattaneo and Farrell (2016a). It also computes alternative estimation and inference procedures available in the literature.

For more details, and related Stata and R packages useful for empirical analysis, visit \url{https://sites.google.com/site/rdpackages/}
}
\usage{
kdrobust(x, c , deriv = 0, p = 2, h = NULL, b = NULL, rho = NULL, 
kernel = "epan", bwselect = "mse", level = 95, all = FALSE, subset = NULL)
}
\arguments{
  \item{x}{the data from which the estimate is to be computed.}
  \item{c}{specifies the evalution point in \code{x}.}
  \item{p}{specifies the order of the kernel used to construct the point-estimator; default is \code{p = 2}.}
  \item{deriv}{specifies the order of the derivative of the density function to be estimated. Default is \code{deriv=0}.}
  \item{h}{specifies the main bandwidth used to construct the KD point estimator. If not specified, bandwidth \code{h} is computed by the companion command \code{\link{kdbwselect}}.}
  \item{b}{specifies the bias bandwidth used to construct the bias-correction estimator. If not specified, bandwidth \code{b} is computed by the companion command \code{\link{kdbwselect}}.}
  \item{rho}{specifies the value of \code{rho}, so that the bias bandwidth \code{b} equals \code{h/rho}. Default is \code{rho = 1} if \code{h} is specified but \code{b} is not.}
  \item{kernel}{is the kernel function used to construct the kernel density estimator(s). Options are \code{triangular} (default option), \code{epanechnikov} and \code{uniform}.}
  
\item{bwselect}{specifies the bandwidth selection procedure to be used. By default it computes both \code{h} and \code{b}, unless \code{rho} is specified, in which case it only computes \code{h} and sets \code{b=h/rho}.}
Options are:

\code{mse} MSE-optimal bandwidth selector for the kernel density estimator.

\code{cer} CER-optimal bandwidth selector for the kernel density estimator.

Note: MSE = Mean Square Error; CER = Coverage Error Rate.
Default is \code{bwselect=mse}. For details on implementation see Calonico, Cattaneo and Titiunik (2014a), Calonico, Cattaneo and Farrell (2016a),
and Calonico, Cattaneo, Farrell and Titiunik (2016), and the companion software articles.

  \item{level}{sets the confidence level for confidence intervals; default is \code{level = 95}.}
  \item{all}{if specified, \code{\link{kdrobust}} reports three different procedures:
  
  (i) conventional LPR estimates with conventional standard errors.
  
  (ii) bias-corrected estimates with conventional standard errors.
  
  (iii) bias-corrected estimates with robust standard errors.}
  
  \item{subset}{an optional vector specifying a subset of observations to be used.}
  
}

\value{
  \item{c}{cutoff value.}
  \item{N}{sample size used.}
  \item{p}{order of the kernel used for estimation of the density function.}
  \item{h}{bandwidth used for estimation of the KDE.}
  \item{b}{bandwidth used for estimation of the bias of the KDE.}
  \item{coef}{vector containing conventional and bias-corrected KDE.}
  \item{se}{vector containing conventional and robust standard errors of the KDE.}
  \item{pv}{vector containing the p-values associated with conventional, bias-corrected and robust KDE.}
  \item{ci}{matrix containing the confidence intervals associated with conventional, bias-corrected and robust KDE.}
}
\references{
Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2016. On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference. Working Paper.
\url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Farrell_2016_JASA.pdf}.

Calonico, S., M. D. Cattaneo, M. H. Farrell, and R. Titiunik. 2016b. rdrobust: Software for Regression Discontinuity Designs. Working Paper.
\url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Farrell-Titiunik_2016_Stata.pdf}.

Calonico, S., M. D. Cattaneo, and R. Titiunik. 2014a. Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs. Econometrica 82(6): 2295-2326.
\url{http://www-personal.umich.edu/~cattaneo/papers/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf}.

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
kdrobust(x, c=0)
}
\keyword{ LPR }
\keyword{ Robust Estimation }
\seealso{
  \code{\link{lpbwselect}}
}