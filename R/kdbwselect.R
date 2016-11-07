### version 0.0.1  06Nov2016

kdbwselect = function(x, c, deriv=0, p=2, kernel="epan", bwselect="mse", all=FALSE, subset = NULL){
  
  if (!is.null(subset)) x <- x[subset]
  na.ok <- complete.cases(x) 
  N=length(x)

  if (kernel=="epanechnikov" | kernel=="epa") {
    kernel_type = "Epanechnikov"
    C=2.34
  }  else if (kernel=="uniform" | kernel=="uni") {
    kernel_type = "Uniform"
    C=1.843
  }   else  {
    kernel_type = "Triangular"
    C=2.576
  }
  
  h.mse.rot = sd(x)*C*N^(-1/(2*p+1))
  h.cer.rot = h.mse.rot*N^(-(p-2)/((1+2*p)*(1+p+2)))
  b.mse.rot=h.mse.rot
  if (bwselect=="dpi") h.cer.dpi = h.cer.den(x,c,h.mse.rot)
                                                                              
  if (all=="FALSE"){
    results = matrix(NA,1,2)
    colnames(results)=c("h","b")
    rownames(results)=bwselect
    if  (bwselect=="mse" | bwselect=="") results[1,] = c(h.mse.rot, b.mse.rot)
    if  (bwselect=="cer") results[1,] = c(h.cer.dpi,  h.cer.dpi)
  }

  if (all=="TRUE"){
    bwselect="All"
    results = matrix(NA,2,2)
    colnames(results)=c("h","b")
    rownames(results)=c("mse","cer") 
    results[1,] =c(h.mse.rot, b.mse.rot)
    results[2,] =c(h.cer.dpi, h.cer.dpi)
  }
  
  tabl1.str=matrix(NA,3,1)
  dimnames(tabl1.str) <-list(c("BW Selector", "Number of Obs", "Kernel Type"), rep("", dim(tabl1.str)[2]))
  tabl1.str[1,1]=bwselect
  tabl1.str[2,1]=N
  tabl1.str[3,1]=kernel_type
  
  tabl2.str=matrix(NA,3,1)
  dimnames(tabl2.str) <-list(c("BW Selector", "Number of Obs", "Kernel Type"), rep("", dim(tabl2.str)[2]))
  tabl2.str[1,1]=bwselect
  tabl2.str[2,1]=N
  tabl2.str[3,1]=kernel_type
  
  
  bws=results
  out = list(tabl1.str=tabl1.str,tabl2.str=tabl2.str,bws=bws,bws,bwselect=bwselect,kernel=kernel_type)
  out$call <- match.call()
  class(out) <- "kdbwselect"
  return(out)
}


#kdbwselect <- function(y,x, ...) UseMethod("kdbwselect")

#kdbwselect.default <- function(y,x, ...){
#  est <- kdbwselectEst(y,x, ...)
#  est$call <- match.call()
#  class(est) <- "kdbwselect"
#  est
#}

print.kdbwselect <- function(x,...){
  cat("Call:\n")
  print(x$call)
  print(x$tabl1.str,quote=F)  
  cat("\n")
  print(x$tabl2.str,quote=F) 
  cat("\n")
  print(x$bws)  
}

summary.kdbwselect <- function(object,...) {
  TAB <- object$bws
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.kdbwselect"
  res
}

#print.summary.kdbwselect <- function(x, ...){
#  cat("Call:\n")
#  print(x$call)
#  cat("\n")
#  printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
#}