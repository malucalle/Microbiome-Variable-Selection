#library(selbal)
  
# Define the function selbal_simulations
  selbal_simulations <- function(x, y, th.imp = 0, covar = NULL, logit.acc="AUC",
                     logt=T, col = c("steelblue1", "tomato1"), tab=T,
                     draw=T, maxV = 1e10, zero.rep = "bayes"){
    
    nulldev0<<-glm(y~1, family=binomial())[[10]]
    
    
  #----------------------------------------------------------------------------#
  # STEP 0: load libraries and extract information
  #----------------------------------------------------------------------------#
    
    #----------------------------------------------#
    # 0.1: information about the response variable
    #----------------------------------------------#
    
    # Class of the response variable
      classy <- class(y)
    # Family for the glm (default: gaussian)
      f.class <- "gaussian"
    # numy to be y and it will be modified if y is a factor
      numy <- y
    # If y is a factor, compute the number of levels of y
      if (classy == "factor"){
        ylev <- levels(y)
        numy <- as.numeric(y) - 1
        f.class <- "binomial"
      }
    
    #------------------------------------------------------------------#
    # 0.2: information and transformation of the independent variables
    #------------------------------------------------------------------#
    
    # Load library
      suppressMessages(library(zCompositions))
    
    # Variables name
      var.nam <- rem.nam <- colnames(x)
    
    # Build a table with the response variable and covariates for correction
      if (!is.null(covar)){ dat <- data.frame(cbind(numy, covar))
      } else { dat <-data.frame(numy)}
    
    # The logCounts (with zero replacement)
      # if (logt == F){ logCounts <- x
      # } else{
        logCounts <- log(cmultRepl2(x, zero.rep = zero.rep))
      # }

    
    #--------------------------------------------------------------------------#
    # 0.3: auxiliar functions
    #--------------------------------------------------------------------------#
    
    #--------------------------------------------------------------------------#
    # Auxiliar function to compute the first balance
    #--------------------------------------------------------------------------#
    
    first.bal<- function(logCounts, Y, covar=NULL){
      
      #------------------------------------------------------------------------#
      # STEP 0: extract information
      #------------------------------------------------------------------------#
      
      # Number and name of variables
        n <- ncol(logCounts)
        nam <- colnames(logCounts)
      
      
      #------------------------------------------------------------------------#
      # STEP 1: build the output matrix
      #------------------------------------------------------------------------#
      
      # The matrix
        if (classy=="factor"){ M<-matrix(0, nrow=n, ncol=n)
        }else{ M<-matrix(1e99, nrow=n, ncol=n)}
      # Row.names and colnames
        row.names(M)<-colnames(M)<-nam
      
      #------------------------------------------------------------------------#
      # STEP 2: complete the matrix
      #------------------------------------------------------------------------#
      
      if(classy=="factor"){
        
        # Solve the problem with libraries
          # suppressWarnings(suppressMessages(library("CMA")))
          # suppressMessages(detach("package:CMA", unload=TRUE))
          suppressMessages(library(pROC))
        
          for (i in 2:n){
            for (j in 1:(i-1)){
            # Build a table with the information
              TAB <- data.frame(cbind(Y,logCounts[,i]-logCounts[,j],covar))
            # Fit the regression model
              FIT <- glm(Y ~ .,data=TAB, family = f.class)
            
            # Add the value into the matrix
              ifelse(FIT$coefficients[2]>0,
                     M[i,j] <- logit.cor(FIT, Y, logit.acc),
                     M[j,i] <- logit.cor(FIT, Y, logit.acc))
            
            } # End j
          } # End i
        
        # Indices for the highest logit.cor value
          r <- which(M == max(M), arr.ind = T)
        
        } else {
          for (i in 2:n){
            for (j in 1:(i-1)){
            # Build a table with the information
              TAB <- data.frame(cbind(Y,logCounts[,i]-logCounts[,j],covar))
            # Fit the regression model
              FIT <- glm(Y ~ .,data=TAB, family = f.class)
            # Complete the matrix
              ifelse(FIT$coefficients[2]>0,
                     M[i,j] <- mean(FIT$residuals^2),
                     M[j,i] <- mean(FIT$residuals^2))
            } # End j
          } # End i
        
        # Indices for the lowest MSE value
          r <- which(M == min(M), arr.ind = T)
        }
      
      
      # Return the row and column of the maximum value
        return(r)
      }
    
    #--------------------------------------------------------------------------#
    
    #--------------------------------------------------------------------------#
    # Auxiliar function to compute the "association value" when adding a new 
    # variable into the balance
    #--------------------------------------------------------------------------#
    
    balance.adding.cor <- function(x, LogCounts, POS, NEG, numy, covar=NULL){
      
      #----------------------------------------#  
      # If x added into the numerator, . .
      #----------------------------------------#
      
      # The "numerator"
        S1.pos <- rowM(LogCounts[,c(POS,x)]); s1 <- length(POS) + 1
      # The "denominator"  
        S2.pos <- rowM(LogCounts[,NEG])     ; s2 <- length(NEG)
      # The balance
        BAL <- sqrt((s1*s2)/(s1+s2))*(S1.pos - S2.pos)
      
      # Data.frame with the variables
        D.pos <- data.frame(cbind(numy, BAL, covar))
      
      # Regression model
        FIT.pos <- glm(numy~., data=D.pos, family=f.class)
      # The MSE or the corresponding value for dichotomous responses  
        if(classy=="numeric"){ C.pos <- mean(FIT.pos$residuals^2)
        }else{ C.pos <- logit.cor(FIT.pos,numy,logit.acc)}
      
      #----------------------------------------#  
      # If x added into the numerator, . .
      #----------------------------------------#
      
      # The numerator    
        S1.neg <- rowM(LogCounts[,POS])       ; s1 <- length(POS)
      # The denominator  
        S2.neg <- rowM(LogCounts[,c(NEG,x)])  ; s2 <- length(NEG) + 1
      # The balance
        BAL <- sqrt((s1*s2)/(s1+s2))*(S1.neg - S2.neg)
      
      # Data.frame with the variables  
        D.neg <- data.frame(cbind(numy, BAL, covar))
      
      # Regression model
        FIT.neg <- glm(numy~., data=D.neg, family=f.class)
      # The MSE or the corresponding value for dichotomous responses
        if(classy=="numeric"){ C.neg <- mean(FIT.neg$residuals^2)
        }else{ C.neg <- logit.cor(FIT.neg,numy,logit.acc)}
      
      # Correlation values
        COR <- c(C.pos, C.neg)
      # Return the values
        return(COR)
      }
    
#------------------------------------------------------------------------------#
    
    
    #--------------------------------------------------------------------------#
    # STEP 1: depending on the response variable class, . . .
    #--------------------------------------------------------------------------#
    
    # Define the first balance
      A1 <- first.bal(logCounts, Y = numy, covar=covar)
    # Variables taking parti into the first balance
      POS <- colnames(x)[A1[1,1]]
      NEG <- colnames(x)[A1[1,2]]
    
    # Included variables in the model
      INC.VAR <- c(POS, NEG)
    # Delete these variables from rem.nam
      rem.nam <- setdiff(rem.nam, INC.VAR)
    
    # Define the initial balance (B_1)
      S1 <- logCounts[,POS]
      S2 <- logCounts[,NEG]
    # Assign the values to B
      B <- sqrt(1/2)*(S1 - S2)
    
    #--------------------------------------------------------------------------#
    # Information about the ACC for the Balance values
    #--------------------------------------------------------------------------#
    # Build a new data.frame
      dat.ini <- cbind(dat, B)
    # Fit the regression model
      FIT.initial <- glm(numy ~ .,data=dat.ini, family = f.class)
    
    # Solve the problem with libraries
      # suppressWarnings(suppressMessages(library("CMA")))
      # suppressMessages(detach("package:CMA", unload=TRUE))
      suppressMessages(library(pROC))
    
    # Define the initial "accuracy" or "association" value
      if(classy=="numeric"){ ACC.Bal <- mean(FIT.initial$residuals^2)
      }else{ ACC.Bal <- logit.cor(FIT.initial, numy, logit.acc)}
    
    
  #----------------------------------------------------------------------------#
    
    # ACC reference
    ACC.ref <- ACC.Bal
    
    #------------------------------------------------------------------------#
    # Improve the balances
    #------------------------------------------------------------------------#
    # Define some parameters
    # The p.value to compare 2 balances (one of them with an additional
    # variable)
      ACC.set <- ACC.ref
    
    # Generate a data.frame with the relevant information
      if(tab){
      # Build a data.frame with the information
        EVOL <- data.frame(POS,NEG,ACC.ref, ACC.ref)
        colnames(EVOL) <- c("NUMERATOR","DENOMINATOR", "ACC", "Increase")
      
      # Index of factor (character) columns
        f.indx <- sapply(EVOL, is.factor)
      }
    # Index of the number of variables for the balance
      nV <- 2
    
    
    #------------------------------#
    # For numeric responses, . . .
    #------------------------------#
    
      if (classy=="numeric"){  
      
      # While there is an improvement and the maximun number of variables 
      # has not been reached, . . .
        # while (ACC.set >= ACC.ref && length(rem.nam)!=0 && nV<maxV){     
        #while (abs(ACC.set - ACC.ref)>1.e-4 && length(rem.nam)!=0 && nV<maxV){     #  *** canvi ***
         while (length(rem.nam)!=0 && nV<maxV){   # canvi
            

        # The new p.bal.ref is the p.set of the previous step
          ACC.ref <- ACC.set
        
        # Function to extract the p-value
         # add2bal.ACC <- matrix(, nrow = 0, ncol = 2)
          add2bal.ACC <- matrix(0, nrow = 0, ncol = 2)  #  *** canvi ***
        
        # Solve the problem with libraries
          # suppressWarnings(suppressMessages(library("CMA")))
          # suppressMessages(detach("package:CMA", unload=TRUE))
          suppressMessages(library(pROC))
        
        
        # Extract the p-values
          add2bal.ACC <- t(sapply(rem.nam, function(x)
            balance.adding.cor(x, LogCounts = logCounts, POS, NEG, numy = numy,
                               covar = covar)))
        # Add names to the rows
          row.names(add2bal.ACC) <- rem.nam
        
        
        # Extract which is the variable (only the first row)
          ACC.opt <- which(add2bal.ACC==min(add2bal.ACC),arr.ind = T)[1,]
        # Modify p.set
          ACC.set <- min(add2bal.ACC)
        
        # If there is an improvement, . . .
          if (abs(ACC.set - ACC.ref) > th.imp){
            INC.VAR <- c(INC.VAR, rem.nam[ACC.opt[1]])
            ACC.Bal <- c(ACC.Bal, ACC.set)
            nV <- nV + 1
            if (ACC.opt[2]==1){
              POS <- c(POS, rem.nam[ACC.opt[1]])
              if(tab){
                EVOL[f.indx] <- lapply(EVOL[f.indx], as.character)
                EVOL <- rbind(EVOL, c(rem.nam[ACC.opt[1]], "-", ACC.set,
                                    ACC.set - ACC.ref))}
            } else if (ACC.opt[2]==2){
              NEG <- c(NEG, rem.nam[ACC.opt[1]])
              if(tab){
                EVOL[f.indx] <- lapply(EVOL[f.indx], as.character)
                EVOL <- rbind(EVOL, c("-", rem.nam[ACC.opt[1]], ACC.set,
                                    ACC.set - ACC.ref))}
          }
        } else {ACC.set <- 0 }
        
        # Remainng variables (possible to add to the balance)
          rem.nam <- rem.nam[-ACC.opt[1]]
        
      } # End while
    }else{ 
      #-----------------------------------#
      # For non-numeric responses, . . .
      #-----------------------------------#  
      
      # While there is an improvement and the maximun number of variables has not
      # been reached, . . .
       # while (ACC.set >= ACC.ref && length(rem.nam)!=0 && nV<maxV){     
         #while (abs(ACC.set - ACC.ref)>1.e-4 && length(rem.nam)!=0 && nV<maxV){     #  *** canvi ***
            
          while (length(rem.nam)!=0 && nV<maxV){   #  *** canvi ***
            
          
        # The new p.bal.ref is the p.set of the previous step
          ACC.ref <- ACC.set
        
        # Function to extract the p-value
          #add2bal.ACC <- matrix(, nrow = 0, ncol = 2)    #  *** canvi ***
          add2bal.ACC <- matrix(0, nrow = 0, ncol = 2)
          
        # Solve the problem with libraries
          # suppressWarnings(suppressMessages(library("CMA")))
          # suppressMessages(detach("package:CMA", unload=TRUE))
          suppressMessages(library(pROC))
        
        # Extract the p-values
          add2bal.ACC <- t(sapply(rem.nam, function(x)
            balance.adding.cor(x, LogCounts = logCounts, POS, NEG, numy = numy,
                               covar = covar)))
        # Add names to the rows
          row.names(add2bal.ACC) <- rem.nam
        
        # Extract which is the variable (only the first row)
          ACC.opt <- which(add2bal.ACC==max(add2bal.ACC),arr.ind = T)[1,]
        # Modify p.set
          ACC.set <- max(add2bal.ACC)
        
        # If there is an improvement, . . .
          if ((ACC.set - ACC.ref) > th.imp){
            # Add the included variable
              INC.VAR <- c(INC.VAR, rem.nam[ACC.opt[1]])
            # Add the Accuracy value  
              ACC.Bal <- c(ACC.Bal, ACC.set)
            # A new variable into the balance  
              nV <- nV + 1
            # Variable included into NUMERATOR or DENOMINATOR?  
              if (ACC.opt[2]==1){
                POS <- c(POS, rem.nam[ACC.opt[1]])
                if(tab){
                  EVOL[f.indx] <- lapply(EVOL[f.indx], as.character)
                  EVOL <- rbind(EVOL, c(rem.nam[ACC.opt[1]], "-", ACC.set,
                                      ACC.set - ACC.ref))}
              } else if (ACC.opt[2]==2){
                NEG <- c(NEG, rem.nam[ACC.opt[1]])
                if(tab){
                  EVOL[f.indx] <- lapply(EVOL[f.indx], as.character)
                  EVOL <- rbind(EVOL, c("-", rem.nam[ACC.opt[1]], ACC.set,
                                        ACC.set - ACC.ref))}
              }
          } else {ACC.set <- 0 }
        
        # Remainng variables (possible to add to the balance)
          rem.nam <- rem.nam[-ACC.opt[1]]
        
        }
      }
    
    # K1 and k2
      k1 <- length(POS); k2 <- length(NEG)
    
    # The final Balance
      FINAL.BAL <- sqrt((k1*k2)/(k1+k2))*
        (rowM(logCounts[,POS])- rowM(logCounts[,NEG]))
  #   
  # #----------------------------------------------------------------------------#
  # # GRAPHICAL REPRESENTATION
  # #----------------------------------------------------------------------------#
  #   
  #   # Load library
  #     library(ggplot2)
  #   
  #   #-----------------------------------------#
  #   # FIRST: The names of included variables
  #   #-----------------------------------------#
  #   
  #   # Variables included
  #     T1 <- c("NUMERATOR", POS)
  #     T2 <- c("DENOMINATOR",NEG)
  #   # Parameter to specify the limits for writting
  #     yl <- max(length(T1), length(T2)) + .5
  #   
  #   # Empty plot with text
  #     df <- data.frame()
  #     Imp.table <- ggplot(df) + xlim(0, 100) + ylim(-0.5, yl) + theme_void() +
  #       annotate("text",
  #                x = 75,
  #                y = floor(yl):ceiling(yl-length(T1)),
  #                label = T1,
  #                colour = c("royalblue1",rep("black",length(T1)-1)),
  #                fontface = 2) +
  #       annotate("text",
  #                x = 25,
  #                y = floor(yl):ceiling(yl-length(T2)),
  #                label = T2,
  #                colour = c("royalblue1",rep("black",length(T2)-1)),
  #                fontface = 2)
  #   
  #   #-----------------------------------------#
  #   # SECOND: The representation of the plots
  #   #-----------------------------------------#
  #   
  #   # Auxiliar data.frame for graphical representation
  #     U <- data.frame(dat, FINAL.BAL)
  #     colnames(U)[ncol(U)] <- "V1"
  #   # Regression model
  #     FIT.final <- glm(numy~., data=U, family = f.class)
  #   
  #   # The plot depending of the class of the response variable
  #     if (classy=="factor"){
  #     
  #     #---------------------------------------------------#
  #     # The composition of the final plot
  #     #---------------------------------------------------#
  #     
  #     # BOXPLOT
  #       BoxP <-  ggplot(U, aes(x=y, y=V1, fill=y)) +
  #         geom_boxplot(color="black", size=1, ) +
  #         scale_fill_manual(values=col) +
  #         theme_bw() +
  #         ylab("Balance") +
  #         xlab("Factor") +
  #         theme(legend.position = "none")
  #     # Density plot for the balance
  #       ydensity <- ggplot(U, aes(V1, fill=y)) +
  #                   geom_density(alpha=.5, size=1.25) +
  #                   scale_fill_manual(values = col) +
  #                   theme_bw() + xlab("") + ylab("") +
  #                   theme(legend.position = "none",
  #                         axis.text.x=element_blank(),
  #                         axis.ticks.x=element_blank()) +
  #                   coord_flip()
  #     
  #     # ROC - curve
  #       library(pROC)
  #     # Build ROC curve
  #       A<-roc(response = U$numy,predictor = FIT.final$fitted.values)
  #     # Extract the sensitivity and specificiti values
  #       ROC.TAB <- data.frame(x=1-A$specificities, y=A$sensitivities)
  #     # Order them for a correct representation
  #       ROC.TAB <- ROC.TAB[order(ROC.TAB$y),]
  #     # AUC value
  #       auc.val <- round(logit.cor(FIT.final,Y = U$numy, logit.acc = logit.acc),3)
  #     # The plot
  #       ROC.plot <- ggplot(data=ROC.TAB, aes(x=x, y=y)) +
  #                   geom_line() +
  #                   ggtitle("ROC curve") +
  #                   xlab("FPR") + ylab("TPR") +
  #                   geom_step() +
  #                   annotate("text", x = .75, y = .2,
  #                             label = paste("AUC-ROC \n",auc.val), size = 2.5) +
  #                   theme_bw() +
  #                   theme(plot.title = element_text(hjust = 0.5))
  #     
  #     # Load libraries
  #       library("gridExtra")
  #       library("grid")
  #       FINAL.P <- arrangeGrob(Imp.table, ROC.plot, BoxP, ydensity,
  #                              ncol=2, nrow=2, widths=c(5,1.25), heights=c(2, 5),
  #                              vp=viewport(width=0.8, height=0.8))
  #     
  #     } else {
  #     
  #       # Fit the regression model
  #         FIT.p <- glm(y ~ FINAL.BAL, family = f.class)
  #     
  #         PLOT.G <- ggplot(U, aes(V1, y)) +
  #                   geom_point(colour = "black", size = 3) +
  #                   geom_abline(intercept=FIT.p$coefficients[1],
  #                               slope=FIT.p$coefficients[2], lwd=3, col="blue") +
  #                   theme_bw() +
  #                   xlab("Balance value") + ylab("Response variable") +
  #                   ggtitle("") +
  #                   theme(strip.text.x = element_text(size=12, angle=0,
  #                                                     face="bold",colour="white"),
  #                         strip.text.y = element_text(size=12, face="bold"),
  #                         strip.background = element_rect(colour="black",
  #                                                         fill="black"),
  #                   plot.title = element_text(size=20, vjust=2.25, hjust= .5,
  #                                             face = "bold"),
  #                   legend.title = element_text(face="bold"),
  #                   legend.text = element_text(face="bold"))
  #     
  #     # Load libraries
  #       library(grid)
  #       library(gridExtra)
  #     
  #       FINAL.P <- arrangeGrob(Imp.table, PLOT.G, nrow=2,
  #                              heights=c(0.2,0.5),vp=viewport(width=0.8,
  #                                                             height=0.8))
  #     
  #     }
  #   
  # # Draw the plot if draw == T  
  #   if (draw){grid.draw(FINAL.P)}
    
    # Round the values
      if(tab){
        EVOL[,3]<-round(as.numeric(EVOL[,3]),5)
        EVOL[,4]<-round(as.numeric(EVOL[,4]),5)
      #   L <- list(FINAL.BAL, POS, NEG, INC.VAR, ACC.Bal, EVOL, FINAL.P,
      #             FIT.final)
      # } else{
      #   L <- list(FINAL.BAL, POS, NEG, INC.VAR, ACC.Bal, FINAL.P,
      #             FIT.final)
        L <- list(FINAL.BAL, POS, NEG, INC.VAR, ACC.Bal, EVOL)
      } else{
        L <- list(FINAL.BAL, POS, NEG, INC.VAR, ACC.Bal)
      }
    
      return(L)
  }
  
  
  ####################################################################################################
  # New   Deviance     #  *** canvi ***
  
  
  
  
  # Define the function logit.cor
  logit.cor <- function(FIT, y, logit.acc){
    if (logit.acc == "AUC"){
      d <- as.numeric(auc(y, FIT$fitted.values))
      # d<-1-(deviance(FIT)/nulldev0)   # proportion of explained deviance
    } else if (logit.acc == "Rsq"){
      d <- cor(y, FIT$fitted.values)^2
    } else if (logit.acc == "Tjur"){
      d <- mean(FIT$fitted.values[y==1]) - mean(FIT$fitted.values[y==0])
    }
    else if (logit.acc == "Dev"){
      d<-1-(deviance(FIT)/nulldev0)  # proportion of explained deviance
    }
    #Return the value
    return(d)
  }
    #------------------------------------------------------------------------------#
  # Auxiliar function in order to replace zeros if necesary
  #------------------------------------------------------------------------------#

  #' Zero replacement for compositional data
  #'
  #' \code{cmultRepl2} replaces the zeros for a matrix where each row is
  #' compositional
  #'
  #'
  #' @param x a \code{matrix} object with the information of variables
  #' (\emph{columns}) for each sample (\emph{rows}).
  #' @param zero.rep if \emph{"bayes"} the Bayesian - Multiplicative treatment
  #' implemented in \code{zCompositions} is applied. If \emph{"one"}, a
  #' pseudocount of 1 is added to the whole matrix.
  #'
  #'
  #'
  #' @return The initial matrix after the zero - replacement normalized so that
  #' each sample's composition sums one.
  #'
  #'
  #' @examples
  #'
  #' # Load the count matrix (with zeros)
  #'   x <- HIV[,1:60]
  #' # Zero replacement
  #'   x.non0 <- cmultRepl2(x, zero.rep = "bayes")
  #'
  #'
  #' @export cmultRepl2



  cmultRepl2 <- function(x, zero.rep = "bayes"){

    # Load library
    library(zCompositions)
    # If there are zeros, use cmultRepl
    if(sum(x==0)!=0){
      if (zero.rep =="bayes"){
        new.x <- cmultRepl(x, suppress.print = T)
      } else if (zero.rep =="one"){
        new.x <- x + 1
      }
    }else { new.x <- x}
    # Return new.x
    return(new.x)
  }


  ################################################################################
  # FUNCTION: rowM
  ################################################################################

  #' Calculates the mean of each row of a matrix (though having only one column).
  #'
  #'
  #'
  #' @param x a \code{matrix} object.
  #'
  #' @return A \code{vector} with the mean of each row of \code{x}, even if the
  #' matrix only contains one column.
  #'
  #' @examples
  #' # Build a matrix with one column
  #'   M <- matrix(rnorm(10), nrow=1)
  #' # rowM (resulting on the same matrix M)
  #'   rowM(M)
  #'
  #' # Build a matrix
  #'   M <- matrix (runif(100), nrow=10)
  #' # Apply rowM function
  #'   rowM(M)
  #'
  #' @export rowM

  # Define the function rowM
  rowM <- function(x){
    if(is.vector(x)) {u <- x
    } else { u <- rowMeans(x)}
    return(u)
  }
