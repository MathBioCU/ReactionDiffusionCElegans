# Generate figure comparing the predicted dominate mode (calculated here) based on the maximum lambda(k) with the mode the simulated dominate mode (calculated using PDE_Simulation.R).

# Set up ------------------------------------------------------------------

library(nleqslv)
library(ggplot2)
setwd('C:/Users/Jacqui/Documents/CU_Boulder/MathBio/ReactionDiffusionCElegans/')
source('Scripts/RDSystemFunctions.R')

# Set up parameter values -------------------------------------------------

a = 1
b = 1
c = 1
km = 0.01
ki = 0.01
h = 0.4
r = 10
L1 = 1
L2 = .2

# Calculated predicted domiante modes -------------------------------------

res = 300
dominateModes_n = matrix(nrow=res,ncol=res)
dominateModes_m = matrix(nrow=res,ncol=res)
d_vector = seq(3,10,length.out=res)
gamma_vector = seq(26,34,length.out=res)
df = data.frame()
for (i in seq(1,res,1)) {
  for (j in seq(1,res,1)) {
    d = d_vector[i]
    gamma = gamma_vector[j]
    ss = nleqslv(c(1,1), dslnex)
    u0 = ss$x[1]
    v0 = ss$x[2]
    derivs = findDerivatives(a,b,c,h,km,ki,r,u0,v0)
    modes_df = findModes(derivs[1],derivs[2],derivs[3],derivs[4],gamma,d,L1,L2)
    if (nrow(modes_df)==0) {
      dominateModes_n[i,j] = 0
      dominateModes_m[i,j] = 0
      df = rbind(df,data.frame(d=d,gamma=gamma,n=0,m=0))
    } else {
      #if (max(modes_df$m)>0) {
      #  print(paste('m mode exists; d=',d,' gamma=',gamma,sep=""))
      #}
      mode_n = subset(modes_df,lambdaReal==max(modes_df$lambdaReal))$n
      mode_m = subset(modes_df,lambdaReal==max(modes_df$lambdaReal))$m
      dominateModes_n[i,j] = mode_n
      dominateModes_m[i,j] = mode_m
      df = rbind(df,data.frame(d=d,gamma=gamma,n=mode_n,m=mode_m))
    }
  }
}

#Save contour plot results
write.csv(df,'Output/Simulation/PredictedModes.csv')

# Graph results -----------------------------------------------------------

#Load contour plot results
df = read.csv('Output/Simulation/PredictedModes.csv')

#Gaussian simulation (seed 7821)
simulatedModes1a = read.csv('Output/Simulation/dominateModes_7821.csv')
simulatedModes1b = read.csv('Output/Simulation/dominateModes_7821_extraTime.csv')
simulatedModes1 = merge(simulatedModes1a,simulatedModes1b,by=c('d','gamma'),all.x=T)
simulatedModes1[is.na(simulatedModes1$mode.y),]$mode.y=simulatedModes1[is.na(simulatedModes1$mode.y),]$mode.x
simulatedModes1$Simulations = as.factor(simulatedModes1$mode.y)
df[["Linear Model"]] = as.factor(df$n)
df[df[["Linear Model"]]==0,][["Linear Model"]]=NA
#Plot predicted and simulate modes zoomed in plot
p1 <- ggplot(data=df) +
  geom_contour(col='black',aes(d, gamma, z = n),binwidth=1) +
  geom_tile(aes(d, gamma, fill=`Linear Model`,alpha=`Linear Model`))+#,binwidth=1) +
  geom_point(data=simulatedModes1,aes(d,gamma,color=Simulations,shape=Simulations),size=2) +
  scale_alpha_manual(values=c(.05,.2))+
  #annotate("text",x=c(9.75,9.75),
  #         
  #         y=c(28.5,33.5),
  #         label=c(1,2),size=10) +
  scale_x_continuous(limits=c(2,10)) +
  scale_y_continuous(limits=c(26,34))+
  ylab(expression(gamma))+
  theme_bw(base_size=22)
p1
ggsave('Output/Simulation/dominateModePredictionZoomed_Gaussian7821.pdf',height=5,width=8)

# Gaussian simulation (seed 2381)
simulatedModes2a = read.csv('Output/Simulation/dominateModes_2381.csv')
simulatedModes2b = read.csv('Output/Simulation/dominateModes_2381_extraTime.csv')
simulatedModes2 = merge(simulatedModes2a,simulatedModes2b,by=c('d','gamma'),all.x=T)
simulatedModes2[is.na(simulatedModes2$mode.y),]$mode.y=simulatedModes2[is.na(simulatedModes2$mode.y),]$mode.x
simulatedModes2$Simulations = as.factor(simulatedModes2$mode.y)
#Plot predicted and simulate modes zoomed in plot
p2 <- ggplot(data=df) +
  geom_contour(col='black',aes(d, gamma, z = n),binwidth=1) +
  geom_tile(aes(d, gamma, fill=`Linear Model`,alpha=`Linear Model`))+#,binwidth=1) +
  geom_point(data=simulatedModes2,aes(d,gamma,color=Simulations,shape=Simulations),size=2) +
  #annotate("text",x=c(9.75,9.75),
  #         y=c(28.5,33.5),
  #         label=c(1,2),size=10) +
  scale_alpha_manual(values=c(.05,.2))+
  scale_x_continuous(limits=c(2,10)) +
  scale_y_continuous(limits=c(26,34))+
  ylab(expression(gamma))+
  theme_bw(base_size=22)
p2
ggsave('Output/Simulation/dominateModePredictionZoomed_Gaussian2381.pdf',height=5,width=8)
