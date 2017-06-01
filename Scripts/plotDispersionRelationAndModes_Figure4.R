# Examining Dispersion Relations and Finding Modes Over Cylindrical domain
#
# 5/17/2017
# Jacqui Wentz

# Set up ------------------------------------------------------------------
setwd('C:/Users/Jacqui/Documents/CU_Boulder/MathBio/ReactionDiffusionCElegans/')
source('Scripts/RDSystemFunctions.R')
library(nleqslv)
library(plot3D)
library(rgl)

# Set parameter values, determine steady state, find derivatives ----------
a = 1
b = 1
c = 1
h = .4
km = 0.01
ki = 0.01
r = 10
L1 = 1
L2 = .2

#Calculate steady state values
ss = nleqslv(c(1,1), dslnex) #Need to define parameters globally before calling nleqslv function.
u0 = ss$x[1]
v0 = ss$x[2]

#Calculate partial derivatives at steady state values
derivs = findDerivatives(a,b,c,h,km,ki,r,u0,v0)
fu = derivs[1]
fv = derivs[2]
gu = derivs[3]
gv = derivs[4]

# Determine critical d value
dc = max(-gv/fu,Re(polyroot(c(gv^2,2*(2*fv*gu - fu*gv),fu^2))))

# FIGURE 4a: 3D Dispersion Relation ---------------------------------------

#Function for plotting 3D dispersion relation
plot3DDispersionRelation = function(d,file) {
  k2 = seq(0,260,.1) #range of k^2 values to include in plot
  gammas = seq(1,100,.1) #range of gamma values to include in plot
  kg = mesh(sqrt(k2),gammas)
  
  #For each value of k2 and gamma calculate the corresponding value of lambda
  lambdaReals = matrix(0,nrow(kg$x),ncol(kg$x))
  for (i in 1:nrow(kg$x)) {
    for (j in 1:ncol(kg$x)) {
      lambdas = calculateLambda(fu,fv,gu,gv,kg$x[i,j],kg$y[i,j],d)
      lambdaReals[i,j] = max(Re(lambdas))
    }
  }
  
  # Find possible modes for each value of gamma
  nm_df_all = data.frame()
  for (gamma in gammas) {
    nm_df = findModes(fu,fv,gu,gv,gamma,d,L1,L2)
    nm_df_all = rbind(nm_df_all,nm_df)
  }
  
  #Remove unnecessarily small values from lambdaReals
  lambdaReals[lambdaReals<(-50)]=NA
  
  #Create color map for the plot
  zlim <- range(lambdaReals,na.rm=T)
  zlen <- zlim[2] - zlim[1] + 1
  colorlut <- rainbow(zlen)
  col <- colorlut[ lambdaReals - zlim[1] + 1 ]
  
  #Create plane at lambda=0 for the plot
  zeroPlane = matrix(0,nrow(kg$x),ncol(kg$x))
  
  um = matrix(c(0.94271821, -0.3187330, 0.0984453,0,
                0.08448815,  0.5136152, 0.8538507,0,
                -0.32271343, -0.7966233, 0.5111236, 0,
                0.00000000,  0.0000000, 0.0000000,1),byrow=T,ncol=4)
  open3d(userMatrix=um)
  plot3d(0,0,0,zlim=c(-40,40),box=FALSE,axes=FALSE,xlab="",ylab="",zlab="",aspect=FALSE)
  surface3d(kg$x^2,kg$y,lambdaReals,alpha=.7,col=col)
  surface3d(kg$x^2,kg$y,zeroPlane,alpha=.7,col='black')
  axis3d('x', pos = c(NA,0,0),cex=1.5)
  axis3d('y++', pos = c(260,NA,0),cex=1.5)
  axis3d('z', pos = c(0,100,NA),cex=1.5)
  title3d(paste("d=",d,sep=""),cex=2)
  plotmath3d(140,10,-45, expression(k^2),cex=1.7)
  plotmath3d(310,60,10, expression(gamma),cex=1.9)
  plotmath3d(-42,100,0, expression(Re(lambda)),cex=1.7)
  points3d(nm_df_all$k2,nm_df_all$gamma,nm_df_all$lambdaReal,col='orange',size=8)
  par3d(windowRect = c(0, 23, 1067, 560),zoom=.7)
  snapshot3d(paste('Output/DispersionRelation/',file,'_d=',d,'.png',sep=""),top=T)
  rgl.close()
  #Use to extract parameters from 3d object (used to find value of userMatrix specified above).
  #userMatrix<-par3d()$userMatrix
  #windowRect<-par3d()$windowRect
}

#Plot dispersion relation for different values of d
for (d in c(3,4,5,10)) {
  plot3DDispersionRelation(d,'dispersionRelation')
}

# FIGURE 4b: Possible Modes -----------------------------------------------

# Calculate the value of w = [u-u0,v-v0] for each mode across the domain
res = 200
modes = 5
gamma = 100
d = 10
nm_df = findModes(fu,fv,gu,gv,gamma,d,L1,L2)
w = array(data = 0, dim = c(L2*res+1,L1*res+1,modes))
for (m in 1:modes) {
  for (i in seq(0,L2*res,1)) {
    for (j in seq(0,L1*res,1))
      w[i+1,j+1,m] = cos(nm_df[m,1]*pi*j/res/L1)*(
        sin(2*nm_df[m,2]*pi*i/res/L2) +
          cos(2*nm_df[m,2]*pi*i/res/L2)
      )
  }
}

# Create meshes for specifying cylindrical surface on which to plot w[u-u0,v-v0]
height = 1
R = 0.2/pi/2
theta = seq(0, 2*pi, length.out=res*.2)
x = R*cos(theta)
y = R*sin(theta)
h = seq(0,height, length.out=res+1)
M1 = mesh(x, h)
M2 = mesh(y, h)

# Plot the five cylindrical surfaces on the same plot for comparison
pdf('Output/Modes/modesOnCylindricalSurface.pdf',height=5,width=7)
offset = c(.3,.15,0,-.15,-.3)
for (m in 1:modes) {
  if (m==1) {add = F} else {add = T}
  surf3D(
    x = M1$x, y = M1$y, z = M2$x+offset[m],
    colvar = unlist(w[1:(res*.2),1:(res+1),m]), colkey = F,
    col = rgb(0,seq(0,1,.01),0), #Use shades of green to color surface
    theta = 60, phi = 0,
    lighting = T, scale = F, add = add,
    zlim = c(-0.3,0.3) #Needs to be specified to make sure there is room for additional surfaces
  )
}
dev.off()

