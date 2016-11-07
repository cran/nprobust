### version 0.0.1  06Nov2016

lpbwselect = function(y, x, c, p=1, q=2, deriv=0, 
                      kernel="epa", bwselect="mse", scaleregul=1, vce="nn",  nnmatch=3, all=FALSE, subset = NULL){
  
  covs=NULL
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
  
  x = as.matrix(x[na.ok])
  y = as.matrix(y[na.ok])
  
  if (!is.null(covs))    covs    = as.matrix(   covs[na.ok])
  if (!is.null(cluster)) cluster = as.matrix(cluster[na.ok])

  if (deriv>0 & p==deriv) {
    p = deriv + 1
    q = p+1
  }

  if (kernel=="epanechnikov" | kernel=="epa") {
    kernel_type = "Epanechnikov"
    C_c=2.34
  }  else if (kernel=="uniform" | kernel=="uni") {
    kernel_type = "Uniform"
    C_c=1.843
  }   else  {
    kernel_type = "Triangular"
    C_c=2.576
  }
  
  if (!is.null(covs)) {
    dZ = ncol(covs)
    Z  = covs
  }
   
  if (!is.null(cluster)) {
    C  = cluster
    g = length(unique(C))
  }
  
  
    if (vce=="nn") {
      order_x = order(x)
      x = x[order_x,,drop=FALSE]
      y = y[order_x,,drop=FALSE]
      if (!is.null(covs))    covs    =    covs[order_x,,drop=FALSE]
      if (!is.null(cluster)) cluster = cluster[order_x,,drop=FALSE]
    }
    
    ### reescaling
    #y_sd = sd(y)
    #y = y/y_sd
    #x_sd = sd(x)
    #x = x/x_sd
    #c_orig = c
    #c = c/x_sd
    x_iq = quantile(x,.75) - quantile(x,.25)
    
    X = x   
    x_min = min(X)
    x_max = max(X)
    range = x_max-x_min
    Y = y
    N = length(X)

    dZ=Z=C=Cind=g=dups=dupsid=covs=cluster=NULL  
    if (vce=="nn") {
      for (i in 1:N) {
        dups[i]=sum(X==X[i])
      }
      for (i in 1:N) {
        dupsid[i:(i+dups[i]-1)]=1:dups[i]
        i=i+dups[i]-1
      }
    }
    
    c_bw = C_c*min(c(1,x_iq/1.349))*N^(-1/5)
    C_d = lprobust_bw(Y, X, Z, C, c=c, o=q+1, nu=q+1, o_B=q+2, h_V=c_bw, h_B=range, 0, vce, nnmatch, kernel, dups, dupsid)
    d_bw = (C_d$V/C_d$B^2)^C_d$rate
    C_b  = lprobust_bw(Y, X, Z, C, c=c, o=q, nu=p+1, o_B=q+1, h_V=c(c_bw), h_B=c(d_bw), scaleregul, vce, nnmatch, kernel, dups, dupsid)
    b_bw = (C_b$V/(C_b$B^2 + scaleregul*C_b$R))^C_b$rate
    C_h  = lprobust_bw(Y, X, Z, C, c=c, o=p, nu=deriv, o_B=q, h_V=c(c_bw), h_B=c(b_bw), scaleregul, vce, nnmatch, kernel, dups, dupsid)
    h_bw = (C_h$V/(C_h$B^2 + scaleregul*C_h$R))^C_h$rate
    #h_mse = x_sd*h_bw
    #b_mse = x_sd*b_bw
    h_mse = h_bw
    b_mse = b_bw
    cer_h = N^(-(p/((3+p)*(3+2*p))))
    cer_b = N^(-(q/((3+q)*(3+2*q))))
  	h_rot = h_mse*cer_h
		b_rot = b_mse*cer_b
		
		results = matrix(NA,1,2)
		colnames(results)=c("h","b")
		rownames(results)=bwselect
		if  (bwselect=="mse" | bwselect=="") results[1,] = c(h_mse, b_mse)
		if  (bwselect=="cer") results[1,] = c(h_rot,  b_rot)
  

  if  (bwselect=="rbc") { 
    h=h_mse
    b=h_mse
    N=length(y)
    rho = h/b
    X.h = (x-c)/h
    X.b = (x-c)/b
    K.h = (0.75*(1-X.h^2)*(abs(X.h)<=1))
    L.b = (0.75*(1-X.b^2)*(abs(X.b)<=1))
    ##############################
    ind.h = K.h>0
    #ind.b = L.b>0
    ind = ind.h
    #if (h>b) ind=ind.h   
    N = sum(ind)
    y=y[ind]
    x=x[ind]
    X.h=X.h[ind]
    X.b=X.b[ind]
    K.h=K.h[ind]
    L.b=L.b[ind]
    #################################
    
    W.p = K.h/h
    W.q = L.b/b
    
    R.p.2 = matrix(NA,N,(p+3))
    for (j in 1:(p+3))  R.p.2[,j] = X.h^(j-1)
    R.p.1 = R.p.2[,1:(p+2)]
    R.p = R.p.2[,1:(p+1)]
    R.q = matrix(NA,N,(q+1))
    for (j in 1:(q+1))  R.q[,j] = X.b^(j-1)
    
    L.p.1 = crossprod(R.p*W.p, X.h^(p+1))/N 
    L.p.2 = crossprod(R.p*W.p, X.h^(p+2))/N 
    L.p.3 = crossprod(R.p*W.p, X.h^(p+3))/N
    L.q.1 = crossprod(R.q*W.q, X.b^(q+1))/N 
    L.q.2 = crossprod(R.q*W.q, X.b^(q+2))/N
    L.q.3 = crossprod(R.q*W.q, X.b^(q+3))/N 
    
    G.p  = t(R.p)%*%(W.p*R.p)/N
    G.q  = t(R.q)%*%(W.q*R.q)/N
    invG.p = solve(G.p)
    invG.q = solve(G.q)
    
    ############################### Residuals ####################################################################
    edups = edupsid = hii = predicts = 0	
    if (vce=="nn") {
      order_x = order(x)
      x = x[order_x]
      y = y[order_x]
      for (i in 1:N) {
        edups[i]=sum(x==x[i])
      }
      for (i in 1:N) {
        edupsid[i:(i+edups[i]-1)]=1:edups[i]
        i=i+edups[i]-1
      }
    }
    
    if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
      H.q = 0
      for (j in 1:(q+1)) H.q[j] = b^(-(j-1))
      beta.q = H.q*invG.q%*%crossprod(R.q*W.q,y)/N
      r.q = matrix(NA,N,q+1)
      predicts = 0
      for (j in 1:(q+1))  r.q[,j] = (x-c)^(j-1)
      for (i in 1:N) predicts[i] = r.q[i,]%*%beta.q
      
      if (vce=="hc2" | vce=="hc3") {
        hii=matrix(NA,N,1)	
        for (i in 1:N) hii[i] = (R.p[i,]%*%invG.p%*%(R.p*W.p)[i,])/N
      }
    }
    res.q = lprobust_res(x, y, Z=NULL, as.matrix(predicts), hii, vce, nnmatch, edups, edupsid, q+1)
    ###################################################################################################
    
    ### Bias (eta term)
    #m.p.3 = lprobust(y=y, x=x, c=c, deriv=(p+3), kernel=kernel, vce=vce)$coef[1] 
    #m.p.2 = lprobust(y=y, x=x, c=c, deriv=(p+2), kernel=kernel, vce=vce)$coef[1] 
    k = p+3
    r_k = matrix(NA,N,k+3)
    for (j in 1:(k+3))  r_k[,j] = x^(j-1)
    gamma = lm(y~r_k-1)
    m.p.3 = gamma$coeff[p+4]*factorial(p+3) + gamma$coeff[p+5]*factorial(p+4)*c + gamma$coeff[p+6]*factorial(p+5)*c^2/2  
    r_k = r_k[,1:(k+2)]
    gamma = lm(y~r_k-1)
    m.p.2 = gamma$coeff[p+3]*factorial(p+2) + gamma$coeff[p+4]*factorial(p+3)*c + gamma$coeff[p+5]*factorial(p+4)*c^2/2 
    
    e.p.1 = matrix(0,(q+1),1); e.p.1[p+2]=1
    e.0   = matrix(0,(p+1),1); e.0[1]=1
    eta.bc1 = (t(e.0)%*%invG.p)%*%(   (m.p.2/factorial(p+2))*L.p.2     + (m.p.3/factorial(p+3))*L.p.3 )
    eta.bc2 = rho^(-2)*b^(q-p-1)*(t(e.0)%*%invG.p)%*%L.p.1%*%t(e.p.1)%*%invG.q%*%( (m.p.2/factorial(p+2))*L.q.1     + (m.p.3/factorial(p+3))*L.q.2 )
    #eta.bc = sqrt(N*h)*h^(p+3)*(eta.bc1-eta.bc2)
    eta.bc = (eta.bc1-eta.bc2)
    ############################################################################################################
    
    ### q_rbc terms
    q_terms = hrbcpp( y = y, x = x,  K=K.h, L=L.b, res=res.q, c = c, p=p, q=q, h=h, b=b)
    q1.rbc=q_terms$q1rbc
    q2.rbc=q_terms$q2rbc
    q3.rbc=q_terms$q3rbc
    ###Solving for optimal rbc bandwidth
    H.bc = function(H) {abs(H^(-1)*q1.rbc + H^(1+2*(p+3))*eta.bc^2*q2.rbc + H^(p+3)*eta.bc*q3.rbc)}
    h.bc <- optimize(H.bc , interval=c(0, 10))
    h_dpi = h.bc$minimum*N^(-1/(p+4))
    b_dpi = b_mse
    results[1,] = c(h_dpi,  b_dpi)
  }
  
  #if (all=="FALSE"){
  #  results = matrix(NA,1,2)
  #  colnames(results)=c("h","b")
  #  rownames(results)=bwselect
  #  if  (bwselect=="mse" | bwselect=="") results[1,] = c(h_mse, b_mse)
  #  if  (bwselect=="cer") results[1,] = c(h_rot, b_rot)
  #}

  if (all=="TRUE"){
    bwselect="All"
    results = matrix(NA,3,2)
    colnames(results)=c("h","b")
    rownames(results)=c("mse","rbc.rot","rbc.dpi") 
    results[1,] =c(h_mse, b_mse)
    results[2,] =c(h_rot, b_rot)
    results[3,] =c(h_dpi, b_dpi)
  }
  
  tabl1.str=matrix(NA,4,1)
  dimnames(tabl1.str) <-list(c("BW Selector", "Number of Obs", "NN Matches", "Kernel Type"), rep("", dim(tabl1.str)[2]))
  tabl1.str[1,1]=bwselect
  tabl1.str[2,1]=N
  tabl1.str[3,1]=nnmatch
  tabl1.str[4,1]=kernel_type
  
  tabl2.str=matrix(NA,3,1)
  colnames(tabl2.str)=c("Left")
  rownames(tabl2.str)=c("Number of Obs","Order Loc Poly (p)","Order Bias (q)")
  tabl2.str[1,]=formatC(c(N),digits=0, format="f")
  tabl2.str[2,]=formatC(c(p),digits=0, format="f")
  tabl2.str[3,]=formatC(c(q),digits=0, format="f")

  bws=results
  out = list(tabl1.str=tabl1.str,tabl2.str=tabl2.str,bws=bws,bws,bwselect=bwselect,kernel=kernel_type,p=p,q=q)
  out$call <- match.call()
  class(out) <- "lpbwselect"
  return(out)
}


#lpbwselect <- function(y,x, ...) UseMethod("lpbwselect")

#lpbwselect.default <- function(y,x, ...){
#  est <- lpbwselectEst(y,x, ...)
#  est$call <- match.call()
#  class(est) <- "lpbwselect"
#  est
#}

print.lpbwselect <- function(x,...){
  cat("Call:\n")
  print(x$call)
  print(x$tabl1.str,quote=F)  
  cat("\n")
  print(x$tabl2.str,quote=F) 
  cat("\n")
  print(x$bws)  
}

summary.lpbwselect <- function(object,...) {
  TAB <- object$bws
  res <- list(call=object$call, coefficients=TAB)
  class(res) <- "summary.lpbwselect"
  res
}

#print.summary.lpbwselect <- function(x, ...){
#  cat("Call:\n")
#  print(x$call)
#  cat("\n")
#  printCoefmat(x$coefficients, P.values=FALSE, has.Pvalue=FALSE)
#}