############################################################################
#   SIMULATION BASED ON CROHN MICROBIOME DATA (with log-contrast)
#
#  You have to choose the name of the scenario parameter file 
#  (default:  namefsimul = 'ParamCrohnSc1')
#  and the name of the output folder for the simulation results
#  (default: ItrLogContrast+namefsimul)
# 
############################################################################
rm(list=ls())

#------------------------------------------------------------------------------#
# Required packages
#------------------------------------------------------------------------------#

# List of required packages
list.of.packages <- c("glmnet","doParallel","foreach","ggplot2","grid","gridExtra","gtable","plyr","pROC","qdapRegex","zCompositions","phyloseq")
# Not installed packages from the list.of.packages
new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
# Install not installed ones
if(length(new.packages)) install.packages(new.packages)


#libraries and source code
suppressMessages(library(glmnet))
source("functions_coda_penalized_regression_simulations.R")
source("functions_selbal_simulations.R")

# activate (=1) or deactivate (=0) one of the methods
SEBAL=1;
Coda_lasso=1;
Clrlasso=1;

#show computation time
printTIME=1; # (=0) to not showing
start_time_TOT<- Sys.time() #time control


# Load Simulation parameters
namefsimul = 'ParamCrohnSc1';
simul_data <- read.csv(paste('./simulParam/',namefsimul,'.csv',sep=""));
dirWrite=paste('./ItrLogContrast',namefsimul,sep="");
dir.create(dirWrite);

meanMethods=NULL;
sdMethods=NULL;

  
  #load DATA
  # # Crohn's disease dataset
  load("Crohn.rda")
  D <- Crohn[,1:48]
  # Data transformation
  p<-ncol(D)   # p: number of covariates after filtering
  # zero imputation
  D = D + 1
  
    
    
  #############################################################################
  # simulation according to read parameters
  #############################################################################
  
  nr=nrow(simul_data);
  nc=ncol(simul_data);
  cat(sprintf("simulation param: nrows= %d ncol= %d \n", nr,nc));
  
  for (iParam in (1:nr)){
    
    params=as.numeric(simul_data[iParam,]);
    niter=params[1];
    ntrue=params[2];
    nfalse=params[3];
    numObs =params[4];
    coeffTrueModel=params[5:nc];
      
    set.seed(12345)
    
    areaAUC_CodaLasso = NULL;
    areaAUC_clrLasso = NULL;
    areaAUC_selbal = NULL;
    
    
  for(k in (1:niter)){ 
      fileout=paste(dirWrite,'/iparam',toString(iParam),'iter',toString(k),'.csv',sep = "")
      fout = file(fileout,"w")
    
      cat(sprintf('\n'));
      cat(sprintf('nPar= %d, niter= %d \n',iParam,k));
      start_iter<- Sys.time() #time control
      
      # DATA simulation
        # At each iteration we generate a new dataset with ntrue positive and nfalse negatives:
        Pos<-sample(1:ncol(D),ntrue)
        Neg<-setdiff(1:ncol(D),Pos)
        Neg<-Neg[sample(1:length(Neg), nfalse)]
        
        # Selection of a submatrix of D
        X<-D[sample((1:nrow(D)),numObs),c(Pos,Neg)]   # we may reduce the number of individuals
    
        # Matrix of proportions
        X<-X/rowSums(X)
        
        # CLR transformation Z=log(X)
        z1<-log(X)
        clrz<- apply(z1,2,function (x) x-rowMeans(z1))
  
        # BETA simulations
        beta_i<-coeffTrueModel;  # coefficients of the true model 
        score1z<-as.matrix(log(X)[,(1:ntrue)])%*%(beta_i)   # linear predictor Z*beta
        
        # Y simulation
        prob1z<-1/(1+exp(-(score1z-mean(score1z))))  #logistic model
        #prob1z<-prob1z-mean(prob1z)+0.5
        y<-rep(0,nrow(z1))
        y<-as.numeric(runif(nrow(z1))<prob1z);
  
        if(SEBAL==1){    # METHOD selbal
        colnames(X)<-(1:ncol(X))
        Xselbal=X;
          start_selbal<- Sys.time() #time control
          y_selbal<-as.factor(y)
          nulldev0 <<- glm(y~1, family=binomial())[[10]];
          selbal_sim<-selbal_simulations(x = Xselbal, y = y_selbal, logt=T, logit.acc="Dev", maxV=ncol(Xselbal), th.imp=0, draw=F)
          D_selbal_X<-selbal_sim[[6]]
  
          denom=as.numeric(D_selbal_X$DENOMINATOR)
          denom[is.na(denom)]=0;
          numerator=as.numeric(D_selbal_X$NUMERATOR)
          numerator[is.na(numerator)]=0;
          
         
         selected_selbal = t(c(numerator[1],denom[1]));
          for (j in (2:length(numerator))){
            selected_selbal = cbind(selected_selbal,max(numerator[j],denom[j]));
          }
          
          matTP=NULL;
          matFP=NULL;
          matTP[1]=0;
          matFP[1]=0;
          for (j in (2:length(c(selected_selbal)))){
            sel_selbal=selected_selbal[1:j];
            matTP[j]<-sum(sel_selbal %in% (1:(ntrue)))/ntrue;   # TP
            matFP[j]<-sum(sel_selbal %in% ((ntrue+1):ncol(X)))/(ncol(X)-ntrue);  # FP
          }
          matTP[length(c(selected_selbal))+1]=1;
          matFP[length(c(selected_selbal))+1]=1;
          
          indXSort = order(matFP);
          sortX=matFP[indXSort];
          sortY=matTP[indXSort];
          
          plot(sortX, sortY, xlim=c(0,1),ylim=c(0,1),type = "l");
          
          areaAUC_selbal = trapezInteg(sortX,sortY);
          end_selbal<- Sys.time()
          cat(sprintf(" Time selbal = %f \n", end_selbal-start_selbal));
        }
        if(Clrlasso == 1){    # METHOD: clr - lasso
          start_crl = Sys.time()
          ROC<-NULL   # ROC matrix 
          lambda_i<-0
          k_i<-0
          seqlambdas<-seq(0,1, by=0.01)   # sequence of lambdas
          M<-matrix(0,ncol=2,nrow=length(seqlambdas))
          
          for (j in seqlambdas){  
            lambda_i<-j
            k_i<-k_i+1
            clrlasso <- glmnet(clrz, as.numeric(y), standardize=FALSE, alpha=1,family="binomial", lambda=lambda_i);
            selected2<-which(as.numeric(abs(coef(clrlasso, lambda=lambda_i)))>0)[-1];  # indices of selected variables (abs(coef)>0)
            
            M[k_i,1]<-sum(selected2 %in% (2:(ntrue+1)))/ntrue
            M[k_i,2]<-sum(selected2 %in% ((ntrue+2):(ncol(X)+1)))/(ncol(X)-ntrue)
            
            ROC<-rbind(ROC,c(lambda_i,M[k_i,1],M[k_i,2]))
          }
      
          colnames(ROC)<-c("lambda","TP_clr_lasso","FP_clr_lasso" )
          ROC<-as.data.frame(ROC)
          # Area for clr lasso
          indSortX = sort(ROC$FP_clr_lasso, index.return=T);
          sortX=ROC$FP_clr_lasso[indSortX[[2]]];
          sortY=ROC$TP_clr_lasso[indSortX[[2]]];
          sortX=c(0,sortX);
          sortY=c(0,sortY);
          sortX=c(sortX,1);
          sortY=c(sortY,1);
          
          plot(sortX, sortY, col =2, xlim=c(0,1),ylim=c(0,1),type = "l")
          
          areaAUC_clrLasso = trapezInteg(sortX,sortY);
          end_crl<- Sys.time()
          cat(sprintf(" Time crl = %f \n", end_crl-start_crl));
        }
        
        if (Coda_lasso == 1){      
        # METHOD: coda - lasso
          start_coda = Sys.time()
          ROC<-NULL   # ROC matrix 
          lambda_i<-0
          k_i<-0
          seqlambdas<-seq(1,0, by=-0.01)   # sequence of lambdas
          M<-matrix(0,ncol=2,nrow=length(seqlambdas))
          for (j in seqlambdas){  
            
            lambda_i<-j
            k_i<-k_i+1
            
            # Model coda logistic lasso
            coda_lasso<-coda_logistic_lasso(y,X,lambda=lambda_i,maxiter=200)
            
            selected1<-coda_lasso[[3]][-1]   # indices of selected variables (abs(coef)>0)
            
            M[k_i,1]<-sum(selected1 %in% (1:ntrue))/ntrue    # TP
            M[k_i,2]<-sum(selected1 %in% ((ntrue+1):ncol(X)))/(ncol(X)-ntrue)  # FP
            ROC<-rbind(ROC,c(lambda_i,M[k_i,1],M[k_i,2]))
          }
          colnames(ROC)<-c("lambda","TP_Coda_lasso","FP_Coda_lasso")
          ROC<-as.data.frame(ROC)
          # compute AUC for each iteration
          # Area for Coda lasso
          indSortX = sort(ROC$FP_Coda_lasso, index.return=T);
          sortX=ROC$FP_Coda_lasso[indSortX[[2]]];
          sortY=ROC$TP_Coda_lasso[indSortX[[2]]];
          sortX=c(0,sortX);
          sortY=c(0,sortY);
          sortX=c(sortX,1);
          sortY=c(sortY,1);
          
          plot(sortX, sortY, xlim=c(0,1),ylim=c(0,1),type = "l")
          
          areaAUC_CodaLasso = trapezInteg(sortX,sortY);
          end_coda =Sys.time()
          cat(sprintf(" Time coda - lasso: %f \n ", as.double(end_coda-start_coda), units = "secs"));
        }
          end_iter=Sys.time();
          cat(sprintf("Time iter = %f \n", as.double(end_iter-start_iter), units = "secs"));
          xx = c(areaAUC_selbal,areaAUC_clrLasso,areaAUC_CodaLasso);
          # dataTot=c(xx);
          writeLines(toString(xx), fout)
          close(fout)
          
    }
  }
  end_timeTOT<-Sys.time()
  cat(sprintf("total simulation time: \n"))
  as.difftime(end_timeTOT- start_time_TOT)
  
   