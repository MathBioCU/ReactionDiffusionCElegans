### Script used to perform projection analysis of data presented in Mendenhall 2015.
### 5/17/2017
### Jacqui Wentz (jacqueline.wentz@colorado.edu)

# Load in required packages and set working directory ---------------------
library(ggplot2)
library(reshape2)
library(fpc) #For clustering data/bootstrap analysis

# Set working directory. Data should be within Data folder of this directory.
setwd('C:/Users/Jacqui/Documents/CU_Boulder/MathBio/ReactionDiffusionCElegans/')

# Load in and visualize raw data ------------------------------------------
data = read.csv('Data/singleReportCopyCElegansData_Mendenhall2015.csv',header=F,row.names=1)
N = 28 #number of worms

#Find total hsp expression in each worm
hspTotal = colSums(data)

#Rearrange data for plotting and analyzing
data_t = as.data.frame(t(data))
data_t$Worm = as.factor(1:N) #Label each worm with a numbers
dataPlot = melt(data_t,id="Worm") #Make a long dataframe where each row corresponds to an intestinal cell

#Label the rows of the dataPlot with the cell position
dataPlot$position = 1
p=1
for (v in unique(dataPlot$variable)) {
  dataPlot[dataPlot$variable==v,]$position=p
  p=p+1
}

#Plot orginial data
ggplot(data=dataPlot)+
  geom_line(aes(position,value,color=Worm)) +
  ylab('HSP Expression')+
  xlab('Intestinal Cell (anterior to posterior)')+
  theme_bw()+
  theme(text=element_text(size=12))
#ggsave('Output/DataVisualization/OriginalData.png',height=6,width=10)

# Average, smooth, and normalize data -------------------------------------

# 1) Average across cells in the same ring
data_t_cellRing = data.frame(cellRing1=rowMeans(data_t[,1:4]))
for (cellRing in seq(2,9)) {
  data_t_cellRing[[paste('cellRing',cellRing,sep="")]] = rowMeans(data_t[,(cellRing*2+1):(cellRing*2+2)])
}

# 2) Take rolling average to smooth data
data_t_cellRing_rollavg = data.frame(cellRing1=character(N))
for (cellRing in seq(1,9)) {
  if (cellRing==1) {
    data_t_cellRing_rollavg[[paste('cellRing',cellRing,sep="")]] = rowMeans(data_t_cellRing[,cellRing:(cellRing+1)])
  } else if (cellRing==9) {
    data_t_cellRing_rollavg[[paste('cellRing',cellRing,sep="")]] = rowMeans(data_t_cellRing[,(cellRing-1):(cellRing)])
  } else {
    data_t_cellRing_rollavg[[paste('cellRing',cellRing,sep="")]] = rowMeans(data_t_cellRing[,(cellRing-1):(cellRing+1)])
  }
}

# 3) Normalize data by mean value for each worm
data_t_cellRing_norm = t(apply(data_t_cellRing_rollavg, 1, function(x) x/mean(x)))
#uncommnet to skip normalization
#data_t_cellRing_norm = as.matrix(data_t_cellRing_rollavg)


# FIGURE 6a: Calculate and plot cosine interpolation ----------------------

# Specify the max mode. The expansion coefficinets (ECs) will be calculated for modes from 0 to the max mode.
maxMode = 8

# Initialize matrix that will contain ECs.
EC_matrix = data.frame(w1=integer(maxMode+1))

# Iterate through each worm to calculate the ECs.
for (w in 1:nrow(data_t_cellRing_norm)) {
  # Get data for current worm
  d1 = data_t_cellRing_norm[w,]
  # Initialize a datafame with evenly spaced values from 1 to 9. 
  df = data.frame(x=seq(1,9,length.out=(maxMode*100)+1))
  # Calculate cos(i pi (x-1)/8) at each of these values where i is the mode number.
  for (i in 0:maxMode) {
    df[[paste('y',i,sep='')]] = cos(i*pi*(df$x-1)/8) 
  }
  # Find the cosine value that corresponds to the discrete spatial data from x=1 to x=9. These correspond to the values for each intestinal cell ring
  df_subset = subset(df,x %in% 1:9)
  # Add the worm data to this dataframe
  df_subset$avgdata = d1
  # Find expansion coefficients for each cosine mode
  EC = c()
  for (i in 0:maxMode) {
    modeData = df_subset[[paste('y',i,sep="")]]
    #Calculate expansion coefficient, depends on mode number
    if (i %in% c(0,8)) {
      EC_current = (d1[1]+d1[9]*modeData[9])/16 + 1/8*sum(modeData[2:8]*df_subset$avgdata[2:8])
    } else {
      EC_current = (d1[1]+d1[9]*modeData[9])/8 + 1/4*sum(modeData[2:8]*df_subset$avgdata[2:8])
    }
    #Add current expansion coefficient to EC vector
    EC = c(EC,EC_current)
    #Calculate the cosine mode weighted by the expansion coefficient at high resolution
    df[[paste('ynew',i,sep="")]]=EC_current*cos(i*pi*(df$x-1)/8)
  }

  # Calculate the sum of all the cosine mode components
  df$ytot = rowSums(df[, grep("*new*", names(df))])
  # Generate and save plot of data with cosine interpolation
  g=ggplot(data=df_subset) +
    geom_point(aes(x,avgdata),color='red') +
    geom_line(data=df,aes(x,ytot)) +
    xlab('Cell Ring in Worm') +
    ylab('Normalized Expression Data') +
    scale_x_continuous(breaks=1:9) +
    theme_bw()
  ggsave(paste('Output/ProjectionAnalysis/worm',w,'.pdf',sep=""),height=3,width=4)
  print(g)
  # Note Figure 6a in paper represents worm 22.
  
  # Add the EC vector to a matrix for all the worm data
  EC_matrix[[paste('w',w,sep="")]] = EC
  
  # Add pause to look at output
  cat ("Press [enter] to continue")
  line <- readline()
}

# Cluster data with kmeans and do bootstrap analysis ----------------------
EC_data = as.data.frame(t(EC_matrix))
colnames(EC_data)=c('EC0','EC1','EC2','EC3','EC4','EC5','EC6','EC7','EC8')
EC_data$hspTotal = hspTotal

#Clustering
EC_cluster = kmeans(EC_data[,1:9],2,nstart=40)

#Boot Strap Analysis
km.boot = clusterboot(EC_data[,1:9], B=100, bootmethod="boot",
                    clustermethod=kmeansCBI,
                    krange=2,seed=992)


# Look for correlation between clusters and total HSP levels --------------

# 1) Visually
hspdf=data.frame(clust=EC_cluster$cluster,hspLevels=c(hspTotal))
plot(hspdf$clust,hspdf$hspLevels)

# 2) Using a t-test
clust1 = subset(hspdf,clust==1)
clust2 = subset(hspdf,clust==2)
t.test(clust1$hspLevels,clust2$hspLevels)

# 3) Calculate mean values
mean1 = mean(clust1$hspLevels)
mean2 = mean(clust2$hspLevels)
#No Significant Difference Detected.


# FIGURE 6b: Plot of expansion coefficients for each cluster --------------
clusterCenters=EC_cluster$centers
ecData = as.data.frame(t(clusterCenters))
ecData$numb = 0:8
ecData=melt(ecData,id="numb")
colnames(ecData)=c("Mode","Cluster","Value")
ecData$Cluster = as.factor(ecData$Cluster)
levels(ecData$Cluster) = c("1","2","A","B")
ecData[ecData$Cluster=="2",]$Cluster = "A"
ecData[ecData$Cluster=="1",]$Cluster = "B"
ggplot(data=ecData) +
  geom_point(aes(Mode,Value,shape=Cluster,color=Cluster),size=2) +
  geom_hline(aes(yintercept=0),linetype=5) +
  scale_x_continuous(breaks=0:8) +
  ylab('Cluster Centers') +
  theme_bw() +
  theme(legend.position=c(.8,.8))#,text=element_text(size=30))
ggsave('Output/ProjectionAnalysis/ExpansionCoefficients.pdf',height=3,width=4)


# FIGURE 6c: Projection results overlaid on individual data ---------------

#Calculate the sum of the different modes for each cluster
int = 8
x1 = seq(1,int+1,.01)
y1 = clusterCenters[1,1]
for (i in 2:(int+1)) {
  y1 = y1 + clusterCenters[1,i]*cos((i-1)*pi/int*(x1-1))
}
y1b = clusterCenters[1,1]+clusterCenters[1,2]*cos(pi/int*(x1-1))+clusterCenters[1,3]*cos(2*pi/int*(x1-1))
y2 = clusterCenters[2,1]
for (i in 2:9) {
  y2 = y2 + clusterCenters[2,i]*cos((i-1)*pi/int*(x1-1))
}
y2b = clusterCenters[2,1]+clusterCenters[2,2]*cos(pi/int*(x1-1))+clusterCenters[2,3]*cos(2*pi/int*(x1-1))

#Organize data for plotting
wormData = as.data.frame(data_t_cellRing_norm)
wormData$worm = as.factor(seq(1,N,1))
wormData$clust = EC_cluster$cluster
dataPlot = melt(data=wormData,id=c("worm","clust"))
dataPlot$position = 1
i=1
for (r in unique(dataPlot$variable)) {
  dataPlot[dataPlot$variable==r,]$position =i
  i=i+1
}
dataPlot$Cluster = 3
dataPlot[dataPlot$clust==1,]$Cluster='B'
dataPlot[dataPlot$clust==2,]$Cluster='A'

#Plot results
ggplot(data=dataPlot) +
  geom_line(aes(x=position,y=value,fill=as.factor(worm),color=Cluster),linetype=5) +
  geom_line(data=data.frame(x=x1,y=y1),aes(x,y,color="B"),linetype=2,size=2) +
  geom_line(data=data.frame(x=x1,y=y1b),aes(x,y,color="B"),size=2) +
  geom_line(data=data.frame(x=x1,y=y2),aes(x,y,color="A"),linetype=2,size=2) +
  geom_line(data=data.frame(x=x1,y=y2b),aes(x,y,color="A"),size=2) +
  theme_bw() +
  theme(legend.position=c(.6,.8))+#,text=element_text(size=30)) +
  scale_x_continuous(breaks=1:9) +
  xlab('Cell Ring') +
  ylab('Normalized HSP Expression')
ggsave('Output/ProjectionAnalysis/FinalPlotWithIndividaulWorms.pdf',height=3,width=4)


# FIGURE 6d: Projection results overlaid on mean/se data ------------------

#Calculate mean and se of hsp values at each cell ring for each cluster
clust1Data = subset(wormData,clust == 1)
clust2Data = subset(wormData,clust == 2)
clust1Means = colMeans(clust1Data[1:9])
clust2Means = colMeans(clust2Data[1:9])
clust1SE = sapply(clust1Data[1:9],function(x) sd(x)/sqrt(length(x)))
clust2SE = sapply(clust2Data[1:9],function(x) sd(x)/sqrt(length(x)))
meandf1 = data.frame(ring=1:9,hsp=clust1Means[1:9],se=clust1SE[1:9],clust=1)
meandf2 = data.frame(ring=1:9,hsp=clust2Means[1:9],se=clust2SE[1:9],clust=2)
meandf = rbind(meandf1,meandf2)
meandf$Cluster = 3
meandf[meandf$clust==2,]$Cluster='A'
meandf[meandf$clust==1,]$Cluster='B'

#Plot results
ggplot(data=dataPlot) +
  geom_errorbar(data=meandf,aes(ring,ymin=hsp-se,ymax=hsp+se,color=Cluster),width=.5,size=1) +
  geom_point(data=meandf,aes(ring,hsp,color=Cluster),size=2) +
  geom_line(data=data.frame(x=x1,y=y1b),aes(x,y,color="B"),size=2) +
  geom_line(data=data.frame(x=x1,y=y2b),aes(x,y,color="A"),size=2) +
  theme_bw() +
  theme(legend.position=c(.6,.8)) +#,text=element_text(size=30)) +
  scale_x_continuous(breaks=1:9) +
  xlab('Cell Ring') +
  ylab('Normalized HSP Expression')
ggsave('Output/ProjectionAnalysis/FinalPlotWithMaxMode.pdf',height=3,width=4)
