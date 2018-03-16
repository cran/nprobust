\name{lpbwselect}
\alias{lpbwselect}
\title{Bandwidth Selection Procedures for Local Polynomial Regression Estimation and Inference}

\description{
\code{\link{lpbwselect}} implements bandwidth selectors for local polynomial regression point estimators and inference procedures developed in Calonico, Cattaneo and Farrell (2018). It also implements other bandwidth selectors available in the literature. See Wand and Jones (1995) and Fan and Gijbels (1996) for background references.

Companion commands: \code{\link{lprobust}} for local polynomial point estimation and inference procedures.

A detailed introduction to this command is given in Calonico, Cattaneo and Farrell (2017). For more details, and related Stata and R packages useful for empirical analysis, visit \url{https://sites.google.com/site/nppackages/}.
}
\usage{
lpbwselect(y, x, eval = NULL, neval = NULL, p = NULL, deriv = NULL,
rho = NULL, kernel = "epa", bwselect = "mse-dpi", 
bwcheck = NULL, bwregul = 1, imsegrid = 30, vce = "nn", nnmatch = 3, 
interior = FALSE, subset = NULL)
}
\arguments{
  \item{y}{dependent variable.}
  
  \item{x}{independent variable.}
  
  \item{eval}{vector of evaluation point(s). By default it uses 30 equally spaced points over to support of \code{x}.}
  
  \item{neval}{number of quantile-spaced evaluation points on support of \code{x}. Default is \code{neval=30}.}
  
  \item{p}{polynomial order used to construct point estimator; default is \code{p = 1} (local linear regression).}
  
  \item{deriv}{derivative order of the regression function to be estimated. Default is \code{deriv=0} (regression function).}
  
  \item{rho}{\code{rho}=\code{h}/\code{b}. Default is \code{rho = 1} if \code{h} is specified but \code{b} is not.}
  
  \item{kernel}{kernel function used to construct local polynomial estimators. Options are \code{epa} for the epanechnikov kernel, \code{tri} for the triangular kernel and \code{uni} for the uniform kernel. Default is \code{kernel = epa}.}
  
  \item{bwselect}{bandwidth selection procedure to be used via \code{\link{lpbwselect}}. By default it computes both \code{h} and \code{b}, unless \code{rho} is specified, in which case it only computes \code{h} and sets \code{b=h/rho}. Options are:
  
    \code{mse-dpi} second-generation DPI implementation of MSE-optimal bandwidth. Default option.
	
	\code{mse-rot} ROT implementation of MSE-optimal bandwidth.
	
	\code{imse-dpi} second-generation DPI implementation of IMSE-optimal bandwidth (computed using grid of evaluation points selected).
	
	\code{imse-rot} ROT implementation of IMSE-optimal bandwidth (computed using grid of evaluation points selected).
	
	\code{ce-dpi} second generation DPI implementation of CE-optimal bandwidth.
	
	\code{ce-rot} ROT implementation of CE-optimal bandwidth.
	
	\code{all} reports all available bandwidth selection procedures.
	
	Note: MSE = Mean Square Error; IMSE = Integrated Mean Squared Error; CE = Coverage Error; DPI = Direct Plug-in; ROT = Rule-of-Thumb. For details on implementation see Calonico, Cattaneo and Farrell (2017).}
	
  \item{bwcheck}{if a positive integer is provided, then the selected bandwidth is enlarged so that at least \code{bwcheck} effective observations are available at each evaluation point. Default is \code{bwcheck = NULL}.}
  
  \item{bwregul}{specifies scaling factor for the regularization term added to the denominator of bandwidth selectors. Setting \code{bwregul = 0} removes the regularization term from the bandwidth selectors. Default is \code{bwregul = 1}.}
  \item{imsegrid}{number of evaluations points used to compute the IMSE bandwidth selector. Default is \code{imsegrid = 30}.}
  
  \item{vce}{procedure used to compute the variance-covariance matrix estimator. Options are:
  
    \code{nn} heteroskedasticity-robust nearest neighbor variance estimator with \code{nnmatch} the (minimum) number of neighbors to be used. Default choice.
	
	\code{hc0} heteroskedasticity-robust plug-in residuals variance estimator without weights.
	
	\code{hc1} heteroskedasticity-robust plug-in residuals variance estimator with \code{hc1} weights.
	
	\code{hc2} heteroskedasticity-robust plug-in residuals variance estimator with \code{hc2} weights.
	
	\code{hc3} heteroskedasticity-robust plug-in residuals variance estimator with \code{hc3} weights.}

  \item{nnmatch}{to be combined with for \code{vce=nn} for heteroskedasticity-robust nearest neighbor variance estimator with \code{nnmatch} indicating the minimum number of neighbors to be used. Default is \code{nnmatch=3}}.
  
  \item{interior}{if TRUE, all evaluation points are assumed to be interior points. This option affects only data-driven bandwidth selection via \code{\link{lpbwselect}}. Default is \code{interior = FALSE}.}
  
  \item{subset}{optional rule specifying a subset of observations to be used.}
}
\value{
\item{Estimate}{A matrix containing \code{grid} (grid points), \code{h} and \code{b} (bandwidths), \code{N}
  (sample size)}
\item{opt}{A list containing options passed to the function.}
}
\references{
Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018. \href{https://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2018_JASA.pdf}{On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference.} Journal of the American Statistical Association, forthcoming.

Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2017. \href{http://sites.google.com/site/nppackages/nprobust/Calonico-Cattaneo-Farrell_2017_nprobust.pdf}{nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference.} Working paper.

Fan, J., and Gijbels, I. 1996. Local polynomial modelling and its applications, London: Chapman and Hall.

Wand, M., and Jones, M. 1995. Kernel Smoothing, Florida: Chapman & Hall/CRC.
}
\author{
Sebastian Calonico, University of Miami, Coral Gables, FL.  \email{scalonico@bus.miami.edu}.

Matias D. Cattaneo, University of Michigan, Ann Arbor, MI.  \email{cattaneo@umich.edu}.

Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}. 
}
\examples{
x   <- runif(500)
y   <- sin(4*x) + rnorm(500)
est <- lpbwselect(y,x)
summary(est)
}
\seealso{
  \code{\link{lprobust}}
}