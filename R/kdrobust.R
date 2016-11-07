### version 0.0.1  06Nov2016

kdrobust = function(x, c, deriv=0, p=2, h=NULL, b=NULL, rho=NULL, 
                    kernel="epan", bwselect="mse", level=95, all=FALSE, subset = NULL) {
  
  if (!is.null(subset)) x <- x[subset]
  na.ok = complete.cases(x) 
  x = x[na.ok]
  
  kernel = tolower(kernel)
  #bwselect = toupper(bwselect)
  
  x_min = min(x);  x_max = max(x)
  N = length(x)
  quant = -qnorm(abs((1-(level/100))/2))
  
  #####################################################   CHECK ERRORS
  exit=0
    #if (kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
    #  print("kernel incorrectly specified")
    #  exit = 1
    #}
    
  if  (bwselect!="mse" & bwselect!="cer" & bwselect!=""){
    print("bwselect incorrectly specified")  
    exit = 1
  }
  
    if (c<=x_min | c>=x_max){
      print("c should be set within the range of x")
      exit = 1
    }
    
    if (deriv<0){
      print("p,q,deriv and matches should be positive integers")
      exit = 1
    }
    
    
    if (level>100 | level<=0){
      print("level should be set between 0 and 100")
      exit = 1
    }
    
    if (!is.null(rho)){  
       if (rho<0){
          print("rho should be greater than 0")
          exit = 1
        }
    }
  
    if (exit>0) stop()
    if (!is.null(h)) bwselect = "Manual"

    
  if (kernel=="epanechnikov" | kernel=="epa") {
    kernel_type = "Epanechnikov"
  }   else if (kernel=="uniform" | kernel=="uni") {
    kernel_type = "Uniform"
  }   else  {
    kernel_type = "Triangular"
  }

  ############################################################################################
    #print("Preparing data.") 
  if (!is.null(h) & is.null(rho) & is.null(b)) {
    rho = 1
    b = h
  }
  if (!is.null(h) & !is.null(rho) ) b = h/rho
  
  if (is.null(h)) {
    kdbws=kdbwselect(x=x, c=c, deriv=deriv, p=p, bwselect=bwselect, kernel=kernel)
    h = kdbws$bws[1]
    b = kdbws$bws[2]
    rho = h/b
  }
  

  ####################################################
  
  u = (x-c)/h
  N = length(x)
  K.d = kern.fun(u     , r=p,   d=deriv, kernel=kernel)
  L.r = kern.fun(rho*u , r=p+2, d=p,     kernel=kernel)
  
  K = K.d$k.x
  M = K - rho^(1+p)*L.r$k.x*L.r$m.k
  
  f.hat.us = mean(K)/h
  f.hat.bc = mean(M)/h
  
  se.us = sqrt((mean((K^2)) - mean(K)^2)/(N*h^2))
  se.bc = sqrt((mean((M^2)) - mean(M)^2)/(N*h^2))

  f.hat = c(f.hat.us, f.hat.bc, f.hat.bc)
  se  = c(se.us, se.us, se.bc)
  t  =  f.hat/se
  pv = 2*pnorm(-abs(t))
  ci = matrix(NA,nrow=3,ncol=2)
  rownames(ci)=c("Conventional","Bias-Corrected","Robust")
  colnames(ci)=c("Lower","Upper")
  ci[1,] = c(f.hat[1] - quant*se[1], f.hat[1] + quant*se[1])
  ci[2,] = c(f.hat[2] - quant*se[2], f.hat[2] + quant*se[2])
  ci[3,] = c(f.hat[3] - quant*se[3], f.hat[3] + quant*se[3])
    
  #print("Estimation Completed.") 
    coef=matrix(f.hat,3,1)
    se  =matrix(se,3,1)
    z   =matrix(t,3,1)
    pv  =matrix(pv,3,1)
    ci=ci

  bws=matrix(c(h,b),1,2)
  rownames(coef)=rownames(se)=rownames(se)=rownames(z)=rownames(pv)=c("Conventional","Bias-Corrected","Robust")
  colnames(coef)="Coeff"
  colnames(se)="Std. Err."
  colnames(z)="z"
  colnames(pv)="P>|z|"
  colnames(bws)=c("left","right")
  rownames(ci)=c("Conventional","Bias-Corrected","Robust")
  colnames(ci)=c("CI Lower","CI Upper")
    
  tabl1.str=matrix(NA,3,1)
  rownames(tabl1.str)=c("Number of Obs", "BW Type", "Kernel Type")
  dimnames(tabl1.str) <-list(c("Number of Obs", "BW Type", "Kernel Type"), rep("", dim(tabl1.str)[2]))
  tabl1.str[1,1]=N
  tabl1.str[2,1]=bwselect
  tabl1.str[3,1]=kernel_type
  
  tabl2.str=matrix(NA,4,1)
  colnames(tabl2.str)=c("Left")
  rownames(tabl2.str)=c("Number of Obs", "BW Loc Poly (h)","BW Bias (b)","rho (h/b)")
  tabl2.str[1,]=formatC(c(N),digits=0, format="f")
  tabl2.str[2,]=formatC(c(h),digits=4, format="f")
  tabl2.str[3,]=formatC(c(b),digits=4, format="f")
  tabl2.str[4,]=formatC(c(h/b),digits=4, format="f")
  
  tabl3.str=matrix("",2,6)
  colnames(tabl3.str)=c("Coef","Std. Err.","z","P>|z|","CI Lower","CI Upper")
  rownames(tabl3.str)=c("Conventional", "Robust")
  tabl3.str[1,1]  =formatC(coef[1],digits=4, format="f")
  tabl3.str[1,2]  =formatC(se[1],  digits=4, format="f")
  tabl3.str[1,3]  =formatC(z[1],   digits=4, format="f")
  tabl3.str[1,4]  =formatC(pv[1],  digits=4, format="f")
  tabl3.str[1,5:6]=formatC(ci[1,], digits=4, format="f")
  tabl3.str[2,4]  =formatC(pv[3],  digits=4, format="f")
  tabl3.str[2,5:6]=formatC(ci[3,] ,digits=4, format="f")
  
  if (all==TRUE){
    tabl3.str=formatC(cbind(coef,se,z,pv,ci),digits=4, format="f")                   
    colnames(tabl3.str)=c("Coef","Std. Err.","z","P>|z|","CI Lower","CI Upper")
  }

  out=list(tabl1.str=tabl1.str,tabl2.str=tabl2.str,tabl3.str=tabl3.str,coef=coef,bws=bws,se=se,z=z,pv=pv,ci=ci,h=h,b=b,rho=rho,N=N)
  out$call <- match.call()
  class(out) <- "kdrobust"
  return(out)
}

#kdrobust <- function(y,x, ...) UseMethod("kdrobust")

#kdrobust.default <- function(y,x,  ...){
#  est <- kdrobustEst(y,x, ...)
#  est$call <- match.call()
#  class(est) <- "kdrobust"
#  est
#}

print.kdrobust <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nSummary:\n")
  print(x$tabl1.str,quote=F)  
  cat("\n")
  print(x$tabl2.str,quote=F)
  cat("\nEstimates:\n")
  print(x$tabl3.str,quote=F)
}

summary.kdrobust <- function(object,...) {
  TAB <- cbind(Estimate    =object$coef,
               "Std. Error"=object$se,
               "z"         =object$z,
               "Pr(>|z|)"  =object$pv,
               "95% CI"    =object$ci)
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.kdrobust"
  res
}

#print.summary.kdrobust <- function(x, ...){
#  cat("Call:\n")
#  print(x$call)
#  cat("\n")
#  printCoefmat(x$coef)
#  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
#}