### version 0.0.1  06Nov2016

kern.fun = function(x,r,d,kernel){
  if (kernel=="norm")  {
    if (d==0)  k = function(u)         dnorm(u)
    if (d==1)  k = function(u)      -u*dnorm(u)
    if (d==2)  k = function(u) (u^2-1)*dnorm(u)
  }
  else if (kernel=="unif")  {
    if (d==0) {
      if (r==2) k = function(u) (abs(u)<=1)*(1/2)
      if (r==4) k = function(u) (abs(u)<=1)*(3/8)*(3-5*u^2)
    }
    if (d==2)  k = function(u) (abs(u)<=1)*(15/4)*(3*u^2-1)
  }
  else if (kernel=="epan")  {
    if (d==0) {
      if (r==2) k = function(u) (abs(u)<=1)*(3/4)*(1-u^2)
      if (r==4) k = function(u) (abs(u)<=1)*(15/32)*(7*u^4-10*u^2+3)
    }
    if (d==2) k = function(u) (abs(u)<=1)*(105/16)*(6*u^2-5*u^4-1)
  }
  else if (kernel=="tria")  {
    k = function(u) (1-abs(u))*(abs(u)<=1)
  }
  
  k.x = k(x)
  m.k = integrate(function(u) (-1)^r*u^(r)*k(u)/factorial(r), -Inf,Inf)$value
  v.k = integrate(function(u) (k(u))^r,   -Inf,Inf)$value
  
  out=list(k.x=k.x, m.k=m.k, v.k=v.k)
  return(out)
}


h.cer.den = function(x,x0,h,b,r,d,kernel) {
  n = length(x)
  u = (x-x0)/h
  rho = h/b
  m.K = kern.fun(u, r=r,d=0, kernel=kernel)$m.k
  f.r.2 =  mean((u^4 -6*u^2 + 3)*dnorm(u))/h^(r+3)
  M.fun = function(u) kern.fun(u, r=r,d=0, kernel=kernel) - 0.5*rho^(1+r)*kern.fun(rho*u,r=r+2,d=r,kernel=kernel)*m.K
  v.fun = function(p,kern) integrate(function(x)   (kern(x))^p,-Inf,Inf)$value
  m.fun = function(m,kern) integrate(function(x) x^(m)*kern(x),-Inf,Inf)$value
  v.M.2 = v.fun(2,M.fun)
  v.M.3 = v.fun(3,M.fun)
  v.M.4 = v.fun(4,M.fun)
  m.M.4 = m.fun(4,M.fun)
  
  z = qnorm(0.975)
  q1 = v.M.4*(z^2-3)/6 - v.M.3^2*(z^4-4*z^2+15)/9
  q2 = f.r.2^2 * m.M.4^2 * v.M.2*z
  q3 = f.r.2   * m.M.4   * v.M.3*(2*z^2)/3

  h.rbc.star <- optimize(function(H) {(H^(-1)*q1 - H^(1+2*(r+2))*q2 + H^(r+2)*q3)^2} , interval=c(0,10))
  h.cer = h.rbc.star$minimum*n^(-1/(r+3))
  return(h.cer)
}







