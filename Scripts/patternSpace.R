patternSpace = function(km,ki,h,r,d,gamma){
  library(nleqslv)
  setwd('C:/Users/Jacqui/Documents/CU_Boulder/MathBio/Aging/ReactionDiffusionModel2016/')
  source('Scripts/RDSystemFunctions.R')
  
  # Need to make parameters global for the nleqslv function to work
  assign("km", km, envir = .GlobalEnv)
  assign("ki", ki, envir = .GlobalEnv)
  assign("h", h, envir = .GlobalEnv)
  assign("r", r, envir = .GlobalEnv)
  assign("a", 1, envir = .GlobalEnv)
  assign("b", 1, envir = .GlobalEnv)
  assign("c", 1, envir = .GlobalEnv)
  
  l1 = 1
  l2 = .2

  # Find steady state solution
  ss = nleqslv(c(1,1), dslnex)
  u0 = ss$x[1]
  v0 = ss$x[2]
  
  # Calculate partial derivatives
  derivs = findDerivatives(a,b,c,h,km,ki,r,u0,v0)
  fu = derivs[1]
  fv = derivs[2]
  gu = derivs[3]
  gv = derivs[4]
  
  # Determine the value of h(k2) at the minimum
  #k2 = gamma*(d*fu+gv)/(2*d)
  #hk2 = d*k2^2 - gamma*(d*(fu+gv))*k2 + gamma^2*(fu*gv-fv*gu)
  
  mode_df = findModes(fu,fv,gu,gv,gamma,d,l1,l2)
  return(nrow(mode_df))
}