###########################################################
#   PLOT RESULTS OF A SIMULATION SCENARIO
#
# Generate the results file for a runned simulation saved in 
# the folderSim (usually it contains 400 csv files)
# The mean and sd are computed and saved in a csv file with 
# the same name. Then it plots the results.
# 
# You have to change the name of the file with the simulation
# parameters scenario (present: ParamCrohnSc1.csv) and the name
# of the folder containing the simulation results.
###########################################################
library(ggplot2)

fileParamScen = 'ParamCrohnSc1.csv' #Param scenario file
folderSim ='ItrKmeansParamCrohnSc1'  #folder of the simulation files results
fileout =paste(folderSim,'.csv',sep="");
fout = file(fileout,"w")
xx= '"mean_Sebal","mean_clr","mean_coda","sd_Sebal","sd_clr","sd_coda"'
writeLines(xx, fout)

for (k in (1:4)) {
  fname=paste('./',folderSim,'/iparam',toString(k),'*.csv', sep="")
  filenames <- Sys.glob(fname)
  
  data=NULL;
  for (i in (1:length(filenames))){
    fileAct = filenames[i];
    dataAUC = read.csv(file= fileAct, sep=",",header=FALSE);
    data=rbind(data,dataAUC);
  }
  
  mm=c(mean(data[[1]]),mean(data[[2]]),mean(data[[3]]));
  ss=c(sd(data[[1]]),sd(data[[2]]),sd(data[[3]]));
  writeLines(toString(c(mm,ss)), fout)
}
close(fout)


# -------------- plot results

  paramDir ='./simulParam';
  dataParam = read.csv(file=paste(paramDir,fileParamScen,sep="/"), sep=",", header=T, comment.char="\\");
  numPos=dataParam$ntrue[1];
  numIndv=dataParam$numObs[1];
  titleAct=sprintf("num positives = %d,  num individuals = %d", numPos, numIndv)
  dataPlot = read.csv(file=fileout, sep=",", header=T, comment.char="\\");
  nc=ncol(dataPlot);
  nr=nrow(dataPlot);
  
  numVar=dataParam$nfalse;
  xplot=rep(numVar,3);
  plotMethod=c(rep("selbal",4),rep("clr lasso",4),rep("coda lasso",4));
  ymean=c(as.matrix(dataPlot[,1:3]));
  ysd=c(as.matrix(dataPlot[,4:6]));
  dataTot <- data.frame(plotMethod,xplot,ymean,ysd)
  
  
  pdog=0.3*min(numVar);
  pd <- position_dodge(pdog) # move them .05 to the left and right
  
  ggplot(dataTot, aes(x=xplot, y=ymean, colour=plotMethod)) + 
    geom_errorbar(aes(ymin=ymean-ysd, ymax=ymean+ysd), width=2, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd) + 
    expand_limits(y=c(0.4,1.1)) +
    labs(x = "number of Negatives", y="mean AUC") +
    labs(colour = "Method") +
    labs(title = titleAct)
  
  lengS=nchar(fileout)
  newFile=substr(fileout,1,lengS-4)
  
  ggsave(paste(newFile,"png",sep='.'),device = "png")
  


  