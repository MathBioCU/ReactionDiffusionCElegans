# Perform Reaction Diffusion Model Simulations
#
# 5/17/2017
# Jacqui Wentz

# Set working directory
setwd('C:/Users/Jacqui/Documents/CU_Boulder/MathBio/ReactionDiffusionCElegans/')

# Loading required libraries and R scripts
library(ReacTran)
library(nleqslv)
#library(magick)
source('Scripts/RDSystemFunctions.R')
#source("Scripts/findGammaForTwoStableModes.R")

# Setting up reaction diffusion model reactions.
# f represents the change in u with time.
# g represents the change in v with time.
f = function(u,v) 1 - a*p(u,v)^r/(h^r+p(u,v)^r) - u
g = function(u,v) 1 - b*p(u,v)^r/(h^r+p(u,v)^r) - c*v
p = function(u,v) v/(km*(1+u/ki)+v)

# Function for running pde on a 2D grid
pde2D <- function(t, y, parms) {
  # Extracting concentrations from last time step and setting them up on a grid
  u <- matrix(nrow = Nx, ncol = Ny, data = y[1:(Nx*Ny)])
  v <- matrix(nrow = Nx, ncol = Ny, data = y[(Nx*Ny+1) : (2*Nx*Ny)])
  # Calculating the change in u and v using zero flux boundary conditions (BC) in x direction and periodic BC in y direction. 
  # For periodic BC the following equation for flux is used: volume fraction * diffusion * difference in concentration/ grid cell length (volume fraction is equal to one).
  du = gamma*f(u,v) + 
    tran.2D(C = u, D.x = 1, D.y = 1,dx = Gridx, dy = Gridy, 
            flux.x.up = 0, 
            flux.x.down = 0,
            flux.y.up = -(u[,Ny] - u[,1])/Gridx$dx[1],
            flux.y.down = -(u[,1] - u[,Ny])/Gridx$dx[1])$dC 
  dv = gamma*g(u,v) +
    tran.2D(C = v, D.x = d, D.y = d,dx = Gridx, dy = Gridy,
            flux.x.up = 0, 
            flux.x.down = 0,
            flux.y.up = -d*(v[,Ny] - v[,1])/Gridx$dx[1],
            flux.y.down = -d*(v[,1] - v[,Ny])/Gridx$dx[1])$dC 
  # Return calculated changes
  list(c(du, dv))
}

# Function for running a single simulation 
runSimulation = function(func,y_ini,time,time_interval) {
  times <- seq(0,time,time_interval)
  print(system.time(
    out <- ode.2D(y = y_ini, parms = NULL, func = pde2D,
                  nspec = 2, dimens = c(Nx, Ny), times = times,
                  lrw = 200000000, names=c("u", "v"),maxsteps=100000)
  ))
  return(out)
}

# Function for saving the simulation results as an image at each timepoint (currently outputs activator, u, concentrations)
saveSimulationResults = function(out,fn,timeInterval) {
  frames = nrow(out)
  for (i in 1:frames) {
    name = paste(fn,i,'.png',sep='')
    png(name,height=570,width=900)
    par(mfrow=c(2,1),cex=2,mar=c(2, 4, 1, 1),oma=c(2,0,1,0))
    image(1:Nx,1:Ny,matrix(out[i,2:(Nx*Ny+1)],Nx,Ny),col=colorRampPalette(c("black", "green"))(100),xlab='',ylab='Circumference')  
    title(main=paste("t*=",(i-1)*timeInterval,sep=""),outer=TRUE)
    plot(1:Nx,out[i,2:(Nx+1)],xlab='Length',ylab='Activator Magnitude')
    mtext('Length', side = 1, outer=T,cex=2)
    dev.off()             
  }       
}

# Convert simulation results to a video
makeVideo = function(input_fn,output_fn) {
  # Load frames
  frames = list()
  for (i in 1:201) {
    frames[[i]] <- image_read(paste0(input_fn,i,".png"))
  }
  # Turn frames into animation
  animation <- image_animate(image_join(frames))
  print(animation)
  
  # Save as GIF
  image_write(animation, output_fn)
}

# Function to set up initial values on the grid. 
initialValues = function(type,seed) {
  set.seed(seed)
  if (type=='Gaussian') {
    u_ini = matrix(nrow = Nx, ncol = Ny, data = c(rep(u0,Nx*Ny)),byrow=T)
    v_ini = matrix(nrow = Nx, ncol = Ny, data = c(rep(v0,Nx*Ny)),byrow=T)
    u_ini = u_ini + matrix(nrow = Nx, ncol = Ny, data = rnorm(Nx*Ny,0,u0*.1))
    v_ini = v_ini + matrix(nrow = Nx, ncol = Ny, data = rnorm(Nx*Ny,0,v0*.1))
  } else if (type=="LeftConstant") {
    u_ini <- matrix(nrow = Nx, ncol = Ny, data = c(rep(u0,Ny),rep(0,(Nx-1)*Ny)),byrow=T)
    v_ini <- matrix(nrow = Nx, ncol = Ny, data = c(rep(v0,Ny),rep(0,(Nx-1)*Ny)),byrow=T)
  } 
  return(c(u_ini,v_ini))
}

# Function used to find the stable model the simulations enter
simulateStableModes = function(outputFile,endTime,drange,gammarange,IC,seed){
  #Find derivatives
  derivs = findDerivatives(a,b,c,h,km,ki,r,u0,v0)
  #Initalize domain
  y_ini = initialValues(IC,seed)
  #Intialize the output file
  write.table(file=outputFile,data.frame(d='d',gamma='gamma',dominateMode='mode'),row.names=F,quote= FALSE, sep=",",col.names=F)
  #Iterate over values of d and gamma to find the mode the simulation stabilizes to
  for (d in drange) {
    print(d)
    assign("d", d, envir=globalenv())
    for (gamma in gammarange) {
      assign("gamma", gamma, envir=globalenv())
      print(gamma)
      #Run the simulation
      out = runSimulation(pde2D,y_ini,endTime,endTime/10)
      saveSimulationResults(out,paste('Output/Simulation/',IC,'/',seed,'/image_',d,'_',gamma,'_',sep=""),endTime/10)
      #Find all stable modes
      modes_df = findModes(derivs[1],derivs[2],derivs[3],derivs[4],gamma,d,L1,L2)
      #See which mode fits the data best (this assumes a 1D pattern, i.e. m=0)
      data1D = out[nrow(out),2:Nx]
      A = (max(data1D)-min(data1D))/2
      data1D_adjusted = abs(data1D-(max(data1D)+min(data1D))/2)
      minsum = 1000
      nfit = 0
      for (n in modes_df$n) {
        tryCatch({
          fit = nls(y ~ abs(A*cos(n*pi*x/L1)), data=data.frame(x=seq(1,99)/99,y=data1D_adjusted), start=list(n=n))
          sfit = summary(fit)
          resid_sum = sum(sfit$residuals^2)
          if (resid_sum < minsum) {
            minsum = resid_sum 
            nfit = n
          }
        },
        error=function(e) {
        }
        )
      }
      print(nfit)
      write.table(file=outputFile,data.frame(d=d,gamma=gamma,dominateMode=nfit),append=T,row.names=F,na="NA",quote= FALSE, sep=",", col.names=F)
    }
  }
}

# Set up global parameters ------------------------------------------------
a = 1
b = 1
c = 1
km = 0.01
ki = 0.01
h = 0.4
r = 10
L1 = 1
L2 = .2
Nx = 100
Ny = 20
Gridx = setup.grid.1D(x.up = 0, x.down = L1, N = Nx)
Gridy = setup.grid.1D(x.up = 0, x.down = L2, N = Ny)

ss = nleqslv(c(1,1), dslnex)
u0 = ss$x[1]
v0 = ss$x[2]
runTime = 2 #Length of simulation in dimensionless time (tstar)

# Find mode of steady-state solution --------------------------------------

drange = seq(4,10,by=1)
gammarange = seq(26,34,by=1)
simulateStableModes('Output/Simulation/dominateModes_7821.csv',runTime,drange,gammarange,"Gaussian",7821) #Output corresponds to Figure 7a.
simulateStableModes('Output/Simulation/dominateModes_2381.csv',runTime,drange,gammarange,"Gaussian",2381) #Output corresponds to Figure 7b.

#Some simulations didn't reach equilibrium. Rerun these with tstar=3
runTime = 3 #Length of simulation in dimensionless time (tstar)

gammarange = c(31,32)
drange = c(4,7,9)
simulateStableModes('Output/Simulation/dominateModes_7821_extraTime.csv',runTime,drange,gammarange,"Gaussian",7821) #Output corresponds to Figure 7a.

gammarange = c(29,30)
drange = c(4,6,7,9,10)
simulateStableModes('Output/Simulation/dominateModes_2381_extraTime.csv',runTime,drange,gammarange,"Gaussian",2381) #Output corresponds to Figure 7a.

# FIGURE 8, 9, and 10: Generate movies for specified simulations ----------
timeInterval = 0.01
runTime = 2
gamma = 30.2
d = 5

#Figure 8
out7821 = runSimulation(pde2D,initialValues("Gaussian",7821),runTime,timeInterval)
saveSimulationResults(out7821,'Output/Simulation/Gaussian/7821/Movie/image',timeInterval)
makeVideo("Output/Simulation/Gaussian/7821/Movie/image","Output/Simulation/Gaussian/7821/Movie/movie.gif")

#Figure 9
out2381 = runSimulation(pde2D,initialValues("Gaussian",2381),runTime,timeInterval)
saveSimulationResults(out2381,'Output/Simulation/Gaussian/2381/Movie/image',timeInterval)
makeVideo("Output/Simulation/Gaussian/2381/Movie/image","Output/Simulation/Gaussian/2381/Movie/movie.gif")

#Figure 10
gamma = 37
d = 4
outLeftConstant = runSimulation(pde2D,initialValues("LeftConstant",9999),runTime,timeInterval)
saveSimulationResults(outLeftConstant,'Output/Simulation/LeftConstant/9999/Movie/image',timeInterval)
makeVideo("Output/Simulation/LeftConstant/9999/Movie/image","Output/Simulation/LeftConstant/9999/Movie/movie.gif")

