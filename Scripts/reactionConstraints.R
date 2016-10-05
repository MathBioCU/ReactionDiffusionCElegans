reactionConstraints = function(km,ki,h,r){
  library(nleqslv)
  setwd('C:/Users/Jacqui/Documents/CU_Boulder/MathBio/Aging/ReactionDiffusionModel2016/')
  source('Scripts/RDSystemFunctions.R')

  assign("km", km, envir = .GlobalEnv)
  assign("ki", ki, envir = .GlobalEnv) #Set ki equal to km
  assign("h", h, envir = .GlobalEnv)
  assign("r", r, envir = .GlobalEnv)
  assign("a", 1, envir = .GlobalEnv)
  assign("b", 1, envir = .GlobalEnv)
  assign("c", 1, envir = .GlobalEnv)
  
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
  
  # Return list containing values of fu + gv and fu*gv-fv*gu
  pattern = 0
  if ((fu < 0 & gv > 0) | (fu > 0 & gv < 0)) {
    if (((fu+gv)<0) & ((fu*gv - fv*gu)>0)) {
      pattern = 1
    }
  }
  return(pattern)
}