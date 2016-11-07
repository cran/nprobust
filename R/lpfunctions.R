### version 0.0.1  06Nov2016

W.fun = function(u,kernel){
  if (kernel=="epanechnikov" | kernel=="epa") w = 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uniform"      | kernel=="uni") w =          0.5*(abs(u)<=1)
  if (kernel=="triangular"   | kernel=="tri") w =   (1-abs(u))*(abs(u)<=1)
  return(w)  
}

k.fun = function(u){
  if (kernel=="epanechnikov" | kernel=="epa") w = 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uniform"      | kernel=="uni") w =          0.5*(abs(u)<=1)
  if (kernel=="triangular"   | kernel=="tri") w =   (1-abs(u))*(abs(u)<=1)
  return(w)  
}

qrXXinv = function(x, ...) {
  #tcrossprod(solve(qr.R(qr(x, tol = 1e-10)), tol = 1e-10))
  #tcrossprod(solve(qr.R(qr(x))))
  chol2inv(chol(crossprod(x)))
}


lprobust_res = function(X, y, Z, m, hii, vce, matches, dups, dupsid, d) {
  n = length(y)
  dZ=0
  if (!is.null(Z)) dZ = ncol(Z)
  res = matrix(NA,n,1+dZ)  	
  
  if (vce=="nn") {
    for (pos in 1:n) {
      rpos = dups[pos] - dupsid[pos]
      lpos = dupsid[pos] - 1
      while (lpos+rpos < min(c(matches,n-1))) {
        if (pos-lpos-1 <= 0) rpos = rpos + dups[pos+rpos+1]
        else if (pos+rpos+1>n) lpos = lpos + dups[pos-lpos-1]
        else if ((X[pos]-X[pos-lpos-1]) > (X[pos+rpos+1]-X[pos])) rpos = rpos + dups[pos+rpos+1]
        else if ((X[pos]-X[pos-lpos-1]) < (X[pos+rpos+1]-X[pos])) lpos = lpos + dups[pos-lpos-1]
        else {
          rpos = rpos + dups[pos+rpos+1]
          lpos = lpos + dups[pos-lpos-1]
        }
      }
      ind_J = (pos-lpos):min(c(n,(pos+rpos)))
      y_J   = sum(y[ind_J])-y[pos]
      Ji = length(ind_J)-1
      res[pos,1] = sqrt(Ji/(Ji+1))*(y[pos] - y_J/Ji)
      if (!is.null(Z)) {
        for (i in 1:dZ) {
          Z_J = sum(Z[ind_J,i])-Z[pos,i]
          res[pos,1+i] = sqrt(Ji/(Ji+1))*(Z[pos,i] - Z_J/Ji)
        }
      }
    }		
  }
  else {
    if (vce=="hc0") w = 1
    else if (vce=="hc1") w = sqrt(n/(n-d))
    else if (vce=="hc2") w = sqrt(1/(1-hii))
    else                 w =      1/(1-hii)
    res[,1] = w*(y-m[,1])
    if (dZ>0) {
      for (i in 1:dZ) {
        res[,1+i] = w*(Z[,i]-m[,1+i])
      }
    }
  }
  return(res)
}


lprobust_bw = function(Y, X, Z, C, c, o, nu, o_B, h_V, h_B, scale, vce, nnmatch, kernel, dups, dupsid){
  dZ = dC = eC = 0
  w = W.fun((X-c)/h_V, kernel)/h_V
  ind_V = w> 0; eY = Y[ind_V];eX = X[ind_V];eW = w[ind_V]
  n_V = sum(ind_V)
  D_V = eY
  R_V = matrix(NA,n_V,o+1)
  for (j in 1:(o+1)) R_V[,j] = (eX-c)^(j-1)
  invG_V = qrXXinv(R_V*sqrt(eW))
  e_v = matrix(0,(o+1),1); e_v[nu+1]=1
  s = 1
  eC=eZ=NULL
  if (!is.null(Z)) {
    dZ = ncol(Z)
    eZ = Z[ind_V,,drop=FALSE]
    D_V = cbind(D_V,eZ)
    U = crossprod(R_V*eW,D_V)
    ZWD  = crossprod(eZ*eW,D_V)
    colsZ = (2):max(c(2+dZ-1,(2)))
    UiGU =  crossprod(U[,colsZ],invG_V%*%U) 
    ZWZ = ZWD[,colsZ] - UiGU[,colsZ] 
    ZWY = ZWD[,1] - UiGU[,1] 
    gamma = chol2inv(chol(ZWZ))%*%ZWY
    s = c(1 , -gamma[,1])
  }
  if (!is.null(C)) {
    dC = 1
    eC =  C[ind_V] 
  }
  beta_V = invG_V%*%crossprod(R_V*eW,D_V)	
  dups_V=dupsid_V=predicts_V=0
  
  if (vce=="nn") {
    dups_V   = dups[ind_V]
    dupsid_V = dupsid[ind_V]
  }
  
  if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
    predicts_V=R_V%*%beta_V
    if (vce=="hc2" | vce=="hc3") {
      hii=matrix(NA,n_V,1)	
      for (i in 1:n_V) {
        hii[i] = R_V[i,]%*%invG_V%*%(R_V*eW)[i,]
      }
    }
  }	
  res_V = lprobust_res(eX, eY, eZ, predicts_V, hii, vce, nnmatch, dups_V, dupsid_V, o+1)
  V_V = (invG_V%*%lprobust_vce(dZ, s, R_V*eW, res_V, eC)%*%invG_V)[nu+1,nu+1]
  v = crossprod(R_V*eW,((eX-c)/h_V)^(o+1))
  Hp = diag(c(1,poly(h_V,degree=o,raw=TRUE)))
  BConst = (Hp%*%(invG_V%*%v))[nu+1]
  
  w = W.fun((X-c)/h_B, kernel)
  ind = w> 0 
  n_B = sum(ind)
  eY = Y[ind];eX = X[ind];eW = w[ind]
  D_B = eY
  R_B = matrix(NA,n_B,o_B+1)
  for (j in 1:(o_B+1)) R_B[,j] = (eX-c)^(j-1)
  invG_B = qrXXinv(R_B*sqrt(eW))
  eC=eZ=NULL
  if (!is.null(Z)) {
    eZ = Z[ind,,drop=FALSE]
    D_B = cbind(D_B,eZ)
  }
  if (!is.null(C)) {
    eC=C[ind]
  }	
  beta_B = invG_B%*%crossprod(R_B*eW,D_B)	
  BWreg=0
  if (scale>0) {
    e_B = matrix(0,(o_B+1),1); e_B[o+2]=1
    dups_B=dupsid_B=hii=predicts_B=0
    if (vce=="nn") {
      dups_B   = dups[ind]
      dupsid_B = dupsid[ind]
    }
    if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
      predicts_B=R_B%*%beta_B
      if (vce=="hc2" | vce=="hc3") {
        hii=matrix(NA,n_B,1)	
        for (i in 1:n_B) {
          hii[i] = R_B[i,]%*%invG_B%*%(R_B*eW)[i,]
        }
      }
    }	
    res_B = lprobust_res(eX, eY, eZ, predicts_B, hii, vce, nnmatch, dups_B, dupsid_B,o_B+1)
    V_B = (invG_B%*%lprobust_vce(dZ, s, R_B*eW, res_B, eC)%*%invG_B)[o+2,o+2]
    BWreg = 3*BConst^2*V_B
  }
  B =  sqrt(2*(o+1-nu))*BConst%*%(t(s)%*%(beta_B[o+2,]))
  V = (2*nu+1)*h_V^(2*nu+1)*V_V
  R = scale*(2*(o+1-nu))*BWreg
  rate = 1/(2*o+3)
  output = list(V=V,B=B,R=R,rate=rate)
  return(output)
}

lprobust_vce = function(d, s, RX, res, C) {	
  k = ncol(RX)
  M = matrix(0,k,k)
  n  = length(C)
  if (is.null(C)) {
    w = 1
    if (d==0){
      M  = crossprod(c(res)*RX)
    }
    else {
      for (i in 1:(1+d)) {
        SS = res[,i]*res
        for (j in 1:(1+d)) {
          M = M + crossprod(RX*(s[i]*s[j])*SS[,j],RX)
        }
      }
    }
  }
  else {	
    clusters = unique(C)
    g     = length(clusters)
    w=((n-1)/(n-k))*(g/(g-1))
    if (d==0){
      for (i in 1:g) {
        ind=C==clusters[i]
        Xi = RX[ind,,drop=FALSE]
        ri = res[ind,,drop=FALSE]
        M = M + crossprod(t(crossprod(Xi,ri)),t(crossprod(Xi,ri)))
      }
    }
    else {
      for (i in 1:g) {
        ind=C==clusters[i]
        Xi = RX[ind,,drop=FALSE]
        ri = res[ind,,drop=FALSE]
        for (l in 1:(1+d)) {	
          for (j in 1:(1+d)) {
            M = M + crossprod(t(crossprod(Xi,s[l]*ri[,l])),t(crossprod(Xi,s[j]*ri[,j])))
          }	
        }					
      }
    }
  }
  return(w*M)		
}





