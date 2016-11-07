### version 0.0.1  06Nov2016

lprobust = function(y, x, c, p=1, q=2, deriv=0, h=NULL, b=NULL, rho=NULL, 
                    kernel="epa", bwselect="mse", scaleregul=1, vce="nn",  nnmatch=3, level=95, all=FALSE, subset = NULL) {
  
  covs = NULL
  cluster=NULL
  
  if (!is.null(subset)) {
    x <- x[subset]
    y <- y[subset]
  }
  
  na.ok <- complete.cases(x) & complete.cases(y)
  
  if (!is.null(cluster)){
    if (!is.null(subset))  cluster <- cluster[subset]
    na.ok <- na.ok & complete.cases(cluster)
  } 
  
  if (!is.null(covs)){
    if (!is.null(subset))  covs <- covs[subset]
    na.ok <- na.ok & complete.cases(covs)
  } 
  
  x = matrix(x[na.ok])
  y = matrix(y[na.ok])
  
  if (!is.null(covs))    covs    = matrix(   covs[na.ok])
  if (!is.null(cluster)) cluster = matrix(cluster[na.ok])
  
  if (vce=="nn") {
    order_x = order(x)
    x = x[order_x,,drop=FALSE]
    y = y[order_x,,drop=FALSE]
    
    if (!is.null(covs))    covs    =    covs[order_x,,drop=FALSE]
    if (!is.null(cluster)) cluster = cluster[order_x,,drop=FALSE]
  }

  kernel = tolower(kernel)
  #bwselect = toupper(bwselect)
  vce = tolower(vce)
  
  x_min = min(x);  x_max = max(x)
  N = length(x)
  range = x_max-x_min
  quant = -qnorm(abs((1-(level/100))/2))
  
  if (deriv>0 & p<=deriv) {
   p = deriv + 1
   q = p+1
  }


  #####################################################   CHECK ERRORS
  exit=0
    if (kernel!="uni" & kernel!="uniform" & kernel!="tri" & kernel!="triangular" & kernel!="epa" & kernel!="epanechnikov" & kernel!="" ){
      print("kernel incorrectly specified")
      exit = 1
    }
    
  if  (bwselect!="mse" & bwselect!="cer" & bwselect!=""){
    print("bwselect incorrectly specified")  
    exit = 1
  }
  
  if (vce!="nn" & vce!="" & vce!="hc1" & vce!="hc2" & vce!="hc3" & vce!="hc0"){ 
    print("vce incorrectly specified")
    exit = 1
  }
    
    if (c<=x_min | c>=x_max){
      print("c should be set within the range of x")
      exit = 1
    }
    
    if (p<0 | q<0 | deriv<0 | nnmatch<=0 ){
      print("p,q,deriv and matches should be positive integers")
      exit = 1
    }
    
    if (p>=q & q>0){
      print("p should be set higher than q")
      exit = 1
    }
    
    if (deriv>p & deriv>0 ){
      print("deriv should be set higher than p")
      exit = 1
    }
    
    p_round = round(p)/p;  q_round = round(q)/q;  d_round = round(deriv+1)/(deriv+1);  m_round = round(nnmatch)/nnmatch
    
    if (p_round!=1 | q_round!=1 | d_round!=1 | m_round!=1 ){
      print("p,q,deriv and matches should be integer numbers")
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
    if (!is.null(h) & is.null(rho) & is.null(b)) {
      rho = 1
      b = h
    }
    if (!is.null(h) & !is.null(rho) ) b = h/rho
    
  if (kernel=="epanechnikov" | kernel=="epa") {
    kernel_type = "Epanechnikov"
  }   else if (kernel=="uniform" | kernel=="uni") {
    kernel_type = "Uniform"
  }   else  {
    kernel_type = "Triangular"
  }

  ############################################################################################
    #print("Preparing data.") 
    if (is.null(h)) {
      #lpbws=lpbwselect(y=y, x=x,  covs=covs, cluster=cluster, c=c, deriv=deriv, p=p, q=q, vce=vce, bwselect=bwselect, kernel=kernel, scaleregul=scaleregul)
      lpbws=lpbwselect(y=y, x=x,  c=c, deriv=deriv, p=p, q=q, vce=vce, bwselect=bwselect, kernel=kernel, scaleregul=scaleregul)
      h = lpbws$bws[1]
      b = lpbws$bws[2]
    }
  ####################################################
  
  w_h = W.fun((x-c)/h, kernel)/h
  w_b = W.fun((x-c)/b, kernel)/b
  ind_h = w_h>0
  ind_b = w_b>0
  N_h = sum(ind_h);  N_b = sum(ind_b)
  
  #if (N_h_l<5 | N_h_r<5 | N_b_l<5 | N_b_r<5){
  #  stop("Not enough observations to perform calculations")
  #  exit(1)
  #}
  
  
  ind = ind_b
  if (h>b) ind = ind_h   
  
  eN = sum(ind)
  eY  = y[ind,,drop=FALSE]
  eX  = x[ind,,drop=FALSE]
  W_h = w_h[ind]
  W_b = w_b[ind]
  
  edups = edupsid = 0	
  if (vce=="nn") {
    for (i in 1:eN) {
      edups[i]=sum(eX==eX[i])
    }
    for (i in 1:eN) {
      edupsid[i:(i+edups[i]-1)]=1:edups[i]
      i=i+edups[i]-1
    }
  }          
          
  u = (eX-c)/h
  R_q = matrix(NA,eN,(q+1))
  for (j in 1:(q+1))  {
    R_q[,j] = (eX-c)^(j-1)
  }
  R_p = R_q[,1:(p+1)]

    #display("Computing RD estimates.")
  L = crossprod(R_p*W_h,u^(p+1)) 
  invG_q  = qrXXinv((sqrt(W_b)*R_q))
  invG_p  = qrXXinv((sqrt(W_h)*R_p))
  e_p1 = matrix(0,(q+1),1); e_p1[p+2]=1
  e_v  = matrix(0,(p+1),1); e_v[deriv+1]=1
  Q_q = t(t(R_p*W_h) - h^(p+1)*(L%*%t(e_p1))%*%t(t(invG_q%*%t(R_q))*W_b))
  D = eY
  eC=eZ=NULL
  dZ = dT = g = 0
  
  if (!is.null(covs)) {
    dZ = ncol(covs)
    Z  = covs
    eZ = Z[ind,,drop=FALSE]
    D  = cbind(D,eZ)
    U_p = crossprod(R_p*W_h,D)
  }
              
  if (!is.null(cluster)) {
    C  = cluster
    eC  = C[ind]
    g = length(unique(eC))
  }
                                                     
  beta_p = invG_p%*%crossprod(R_p*W_h,D); beta_q = invG_q%*%crossprod(R_q*W_b,D); beta_bc = invG_p%*%crossprod(Q_q,D) 


  if (is.null(covs)) {	
  tau_cl = tau_Y_cl = factorial(deriv)*beta_p[(deriv+1),1]
  tau_bc = tau_Y_bc = factorial(deriv)*beta_bc[(deriv+1),1]
  s_Y = 1
  } else {	
    ZWD_p  = crossprod(eZ*W_h,D)
    colsZ = (2+dT):max(c(2+dT+dZ-1,(2+dT)))
    UiGU_p =  crossprod(U_p[,colsZ],invG_p%*%U_p) 
    ZWZ_p = ZWD_p[,colsZ] - UiGU_p[,colsZ] 
    ZWY_p = ZWD_p[,1:(1+dT)] - UiGU_p[,1:(1+dT)] 
    gamma_p = chol2inv(chol(ZWZ_p))%*%ZWY_p
    s_Y = c(1 ,  -gamma_p[,1])
    tau_cl = t(s_Y)%*%beta_p[(deriv+1),]
    tau_bc = t(s_Y)%*%beta_bc[(deriv+1),]
  }

#*display("Computing variance-covariance matrix.")
  
  hii=predicts_p=predicts_q=0
  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    predicts_p=R_p%*%beta_p
    predicts_q=R_q%*%beta_q
    if (vce=="hc2" | vce=="hc3") {
      hii=matrix(NA,eN,1)	
      for (i in 1:eN) {
        hii[i] = R_p[i,]%*%invG_p%*%(R_p*W_h)[i,]
      }
    }
  }
  						
	res_h = lprobust_res(eX, eY, eZ, predicts_p, hii, vce, nnmatch, edups, edupsid, p+1)
	if (vce=="nn") {
			res_b = res_h
	} 	else {
			res_b = lprobust_res(eX, eY, eZ, predicts_q, hii, vce, nnmatch, edups, edupsid, q+1)
  }
			                       
	V_Y_cl = invG_p%*%lprobust_vce(dT+dZ, s_Y, R_p*W_h, res_h, eC)%*%invG_p
	V_Y_bc = invG_p%*%lprobust_vce(dT+dZ, s_Y, Q_q,       res_b, eC)%*%invG_p
	V_tau_cl = factorial(deriv)^2*V_Y_cl[deriv+1,deriv+1]
	V_tau_rb = factorial(deriv)^2*V_Y_bc[deriv+1,deriv+1]
	se_tau_cl = sqrt(V_tau_cl);	se_tau_rb = sqrt(V_tau_rb)

  tau = c(tau_cl, tau_bc, tau_bc)
  se  = c(se_tau_cl,se_tau_cl,se_tau_rb)
  t  =  tau/se
  pv = 2*pnorm(-abs(t))
  ci = matrix(NA,nrow=3,ncol=2)
  rownames(ci)=c("Conventional","Bias-Corrected","Robust")
  colnames(ci)=c("Lower","Upper")
  ci[1,] = c(tau[1] - quant*se[1], tau[1] + quant*se[1])
  ci[2,] = c(tau[2] - quant*se[2], tau[2] + quant*se[2])
  ci[3,] = c(tau[3] - quant*se[3], tau[3] + quant*se[3])
    
  #print("Estimation Completed.") 
    coef=matrix(tau,3,1)
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
    
  tabl1.str=matrix(NA,4,1)
  rownames(tabl1.str)=c("Number of Obs", "NN Matches", "BW Type", "Kernel Type")
  dimnames(tabl1.str) <-list(c("Number of Obs", "NN Matches", "BW Type", "Kernel Type"), rep("", dim(tabl1.str)[2]))
  tabl1.str[1,1]=N
  tabl1.str[2,1]=nnmatch
  tabl1.str[3,1]=bwselect
  tabl1.str[4,1]=kernel_type
  
  tabl2.str=matrix(NA,6,1)
  colnames(tabl2.str)=c("Left")
  rownames(tabl2.str)=c("Number of Obs","Order Loc Poly (p)","Order Bias (q)","BW Loc Poly (h)","BW Bias (b)","rho (h/b)")
  tabl2.str[1,]=formatC(c(eN),digits=0, format="f")
  tabl2.str[2,]=formatC(c(p),digits=0, format="f")
  tabl2.str[3,]=formatC(c(q),digits=0, format="f")
  tabl2.str[4,]=formatC(c(h),digits=4, format="f")
  tabl2.str[5,]=formatC(c(b),digits=4, format="f")
  tabl2.str[6,]=formatC(c(h/b),digits=4, format="f")
  
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

  out=list(tabl1.str=tabl1.str,tabl2.str=tabl2.str,tabl3.str=tabl3.str,coef=coef,bws=bws,se=se,z=z,pv=pv,ci=ci,p=p,q=q,h=h,b=b,rho=rho,N=N)
  out$call <- match.call()
  class(out) <- "lprobust"
  return(out)
}

#lprobust <- function(y,x, ...) UseMethod("lprobust")

#lprobust.default <- function(y,x,  ...){
#  est <- lprobustEst(y,x, ...)
#  est$call <- match.call()
#  class(est) <- "lprobust"
#  est
#}

print.lprobust <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nSummary:\n")
  print(x$tabl1.str,quote=F)  
  cat("\n")
  print(x$tabl2.str,quote=F)
  cat("\nEstimates:\n")
  print(x$tabl3.str,quote=F)
}

summary.lprobust <- function(object,...) {
  TAB <- cbind(Estimate    =object$coef,
               "Std. Error"=object$se,
               "z"         =object$z,
               "Pr(>|z|)"  =object$pv,
               "95% CI"    =object$ci)
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.lprobust"
  res
}

#print.summary.lprobust <- function(x, ...){
#  cat("Call:\n")
#  print(x$call)
#  cat("\n")
#  printCoefmat(x$coef)
#  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
#}