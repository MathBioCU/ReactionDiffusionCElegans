# Helper functions for reaction diffusion model

#Called by nleqslv to calculate the steady state solution (u0,v0). Not sure how to pass parameters to this function so the parameter values need to be defined when dslnex is called.
dslnex = function(x) {
  y = numeric(2)
  p = x[2]/(km*(1+x[1]/ki)+x[2])
  y[1] = 1 - a*p^r/(h^r+p^r) - x[1]
  y[2] = 1 - b*p^r/(h^r+p^r) - c*x[2]
  y #Return y
}

#Calculates the derivate of the reaction fucntions at (u0,v0)
findDerivatives = function(a,b,c,h,Km,Ki,r,u0,v0) {
  Z = v0^r+h^r*(Km*(1+u0/Ki) + v0)^r
  fu = r*Km*h^r*a*v0^r*(Km*(u0/Ki+1)+v0)^(r-1)/Ki/Z^2 - 1
  fv = -a*r*v0^(r-1)/Z + a*v0^r*(r*v0^(r-1)+r*h^r*(v0+Km*(u0/Ki+1))^(r-1))/Z^2
  gu = r*h^r*Km*b*v0^r*(Km*(u0/Ki+1)+v0)^(r-1)/Ki/Z^2
  gv = -b*r*v0^(r-1)/Z + b*v0^r*(r*v0^(r-1)+r*h^r*(v0+Km*(u0/Ki+1))^(r-1))/Z^2 - c
  return(c(fu,fv,gu,gv))
}

#Calculate the dispersion relation (i.e., lambda as function of wavenumber k)
calculateLambda = function(fu,fv,gu,gv,k,gamma,d) {
  a = 1
  b = k^2*(1+d) - gamma*(fu+gv)
  c = d*k^4 - gamma*(d*fu+gv)*k^2 + gamma^2*(fu*gv - fv*gu)
  root1 = (-b+sqrt(b^2-4*a*c))/(2*a)
  root2 = (-b-sqrt(b^2-4*a*c))/(2*a)
  return(c(root1,root2))
}

#Find all modes that are unstable (i.e., dispersion relation at the specified wavenumber is greater than zero). Depends on calculateLambda(). Currently finds up to mode n=10 and m=10.
findModes = function(fu,fv,gu,gv,gamma,d,p,q) {
  k12 = gamma/2/d*((d*fu+gv) - ((d*fu+gv)^2 - 4*d*(fu*gv-fv*gu))^(1/2))
  k22 = gamma/2/d*((d*fu+gv) + ((d*fu+gv)^2 - 4*d*(fu*gv-fv*gu))^(1/2))
  nm_df = data.frame(n=c(),m=c())
  if (is.na(k12) | is.na(k22)) {
    return(nm_df)
  }
  for (n in seq(0,5,1)) {
    for (m in seq(0,1,1)) {
      k2 = pi^2*(n^2/p^2 + 4*m^2/q^2)
      if ((k2 > k12) & (k2 < k22)) {
        lambdaReal = max(Re(calculateLambda(fu,fv,gu,gv,sqrt(k2),gamma,d)))
        nm_df = rbind(nm_df,data.frame(n=n,m=m,k2=k2,gamma=gamma,lambdaReal=lambdaReal))
      }
    }
  }
  #nm_df$k = sqrt(pi^2*(nm_df$n^2/p^2 + 4*nm_df$m^2/q^2))
  #nm_df$w = 2*pi/nm_df$k
  return(nm_df)
}

