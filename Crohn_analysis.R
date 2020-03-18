#####################################################################
#
# Variable selection in microbiome compositional data analysis
# Antoni Susin, Yiwen Wang, Kim-Anh Lê Cao, M.Luz Calle (2020)
#
# Methods: selbal, clr-lasso, coda-lasso
#
# Crohn_analysis.R
#
#####################################################################

cat("\014") #clean console window

dev.off() #clean plots window

#------------------------------------------------------------------------------#
# Required packages
#------------------------------------------------------------------------------#

# List of required packages
list.of.packages <- c("VennDiagram","glmnet","selbal")
# Not installed packages from the list.of.packages
new.packages <- list.of.packages[!(list.of.packages %in%
                                     installed.packages()[,"Package"])]
# Install not installed ones
if(length(new.packages)) install.packages(new.packages)

# Load library
suppressMessages(library(VennDiagram))
suppressMessages(library(glmnet))
suppressMessages(library(selbal))


#https://github.com/UVic-omics/CoDA-Penalized-Regression/blob/master/R/
source("./../../main_code/functions_coda_penalized_regression.R")

#https://github.com/UVic-omics/Microbiome-Variable-Selection
source("./../../main_code/functions_tutorial.R")


# Compute Methods

use_Selbal =1
use_ClrLasso =1
use_CodaLasso =1


#------------------------------------------------------------------------------#
#
# Read data "Crohn_data.rda"
#
#------------------------------------------------------------------------------#


load("Crohn_data.rda")


#------------------------------------------------------------------------------#
#
# selbal: Selection of balances
#
#------------------------------------------------------------------------------#

if (use_Selbal == 1) {

  y<-y_Crohn  # for binary outcomes (logistic regression) selbal requires y to be factor,
  #                             if y is numeric selbal implements linear regression

  x<-x_Crohn
  colnames(x)<-(1:ncol(x))  # (12 variables, 27% explained variance)

  selbal_Crohn<-selbal(x = x, y = y, logt=T, maxV=12, logit.acc="Dev") #optimization criteria Deviance

  dev.off() # clean plots window
  grid.draw(selbal_Crohn $global.plot2)

  CD.results_selbal <- selbal_wrapper(Y = y_Crohn, X = x_Crohn, maxV = 12, logit.acc = 'Dev')

  CD.results_selbal$numVarSelect

  CD.results_selbal$varSelect

}



#------------------------------------------------------------------------------#
#
# clr-lasso: CLR + LOGISTIC LASSO
#
#------------------------------------------------------------------------------#

if (use_ClrLasso == 1) {

X<-x_Crohn
y<-as.numeric(y_Crohn)   # requires y to be numeric

# CLR transformation Z=log(X)
z<-log(X)
z<- apply(z,2,function (x) x-rowMeans(z))
CD.clrx<-z
#rowMeans(z)

glmnet(z, y, family="binomial")
glmnet(z, y, standardize=T, alpha=1,family="binomial", lambda=seq(0.03, 0.05, 0.001))

lambda_clr=0.043  # 12 variables, 18,85% explained variance

clrlasso_Crohn <- glmnet(z, y, standardize=T, alpha=1,family="binomial", lambda=lambda_clr)

CD.lambda_clr = 0.043
CD.results_clrlasso <- glmnet_wrapper(Y = y_Crohn_numeric, X = CD.clrx, family = 'binomial',
                                      lambda = CD.lambda_clr)
CD.results_clrlasso$numVarSelect

CD.results_clrlasso$varSelect

CD.results_clrlasso$varIndex


}



#------------------------------------------------------------------------------#
#
# coda-lasso: CODA LOGISTIC LASSO
#
#------------------------------------------------------------------------------#

if (use_CodaLasso == 1) {

  # Main function:
  #   coda_logistic_lasso(y,x,lambda)

  # y is the binary outcome, can be numerical (values 0 and 1), factor (2 levels) or categorical (2 categories)
  #
  # x is the matrix of microbiome abundances, either absolute abundances (counts) or relative abundances (proportions)
  #   the rows of x are individuals/samples, the columns are taxa

  # Imputation of zeros:
  # The user can provide x without zeros (after implementation of an external zero imputation method) and no additional transformation will be applied
  # If x contains zeros,

  # Log-transformation:
  # x should not be the matrix of log(counts) or log(proportions). The method itself performs the log-transformation of the abundances.

  # lambda is the penalization parameter: the larger the value of lambda the fewer number of variables will be selected.
  #
  # function rangLambda(y,x,numVar, lambdaIni=1) provides a rang of lambda values corresponding to a given number of variables to be selected (numVar). The default initial lambda is lambdaIni=1.

  y<-y_Crohn

  x<-x_Crohn

  lambdaRange_codalasso(y,x)

  # lambda=0.15  (12 variables, 21,3% explained variance)
  results_codalasso<-coda_logistic_lasso(y,x,lambda=0.15)

  results_codalasso$`number of iterations`
  results_codalasso$`number of selected variables`
  results_codalasso$`indices of selected variables`
  results_codalasso$`name of selected variables`
  results_codalasso$`beta non-zero coefficients`
  results_codalasso$`proportion of explained deviance`
  results_codalasso$betas



}
#######################################################################################################
#
# Concordance of variables selected by the three methods
#
#######################################################################################################



d_selbal<-read.csv("results_Selbal_Crohn.csv", sep=",", header = T)
d_clrlasso<-read.csv("results_clrlasso_Crohn.csv", sep=",", header = T)
d_codalasso<-read.csv("results_codalasso_Crohn.csv", sep=",", header = T)

taxa_selected<-list(d_clrlasso[,2], d_codalasso[,2], d_selbal[,2])
taxa.id_selected<-list(d_clrlasso[,3], d_codalasso[,3], d_selbal[,3])


venn.plot <- venn.diagram(taxa.id_selected , NULL, fill=c("magenta", "blue", "lightblue"),
                          alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4,
                          category.names=c("clr+lasso", "coda-lasso", "selbal"),
                          main="Concordance of selected taxa for Crohn data")
grid.draw(venn.plot)


