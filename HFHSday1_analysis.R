#####################################################################
#
# Variable selection in microbiome compositional data analysis
# Antoni Susin, Yiwen Wang, Kim-Anh Lê Cao, M.Luz Calle (2020)
#
# Methods: selbal, clr-lasso, coda-lasso
#
# HFHSday1_analysis.R
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

use_CodaLasso =1
use_ClrLasso =1
use_Selbal =1



#------------------------------------------------------------------------------#
#
# Read data "HFHSday1.rda"
#
#------------------------------------------------------------------------------#


load("HFHSday1.rda")

#------------------------------------------------------------------------------#
#
# selbal: Selection of balances
#
#------------------------------------------------------------------------------#

if (use_Selbal == 1) {

  x <- X_day1

  y <- as.factor(y_day1)


  selbal_HFHS <- selbal(x = x, y = y, logt=T, maxV=2, logit.acc="Dev")

  dev.off() # clean plots window
  grid.draw(selbal_HFHS$global.plot2)


  HFHS.results_selbal <- selbal_wrapper(selbal_HFHS, X_day1)
  HFHS.results_selbal$numVarSelect

  HFHS.results_selbal$varSelect

  selected<- HFHS.results_selbal$varSelect

  HFHS.tax_selbal <- taxonomy_tab[which(rownames(taxonomy_tab)%in% selected),]

  HFHS.tax_selbal[,2:6]
}



#------------------------------------------------------------------------------#
#
# clr-lasso: CLR + LOGISTIC LASSO
#
#------------------------------------------------------------------------------#

if (use_ClrLasso == 1) {

X <- X_day1
y <- as.numeric(y_day1)   # requires y to be numeric

# CLR transformation Z=log(X)
z <- log(X)
z <- apply(z,2,function (x) x-rowMeans(z))

HFHS.clrx<-z
HFHS.y.num<-y
#rowMeans(z)


HFHS.test_clrlasso <- glmnet(x = HFHS.clrx, y = HFHS.y.num, family = 'binomial')
plot(HFHS.test_clrlasso, xvar = 'lambda', label = T)

HFHS.test_clrlasso
#cv.clrlasso <- cv.glmnet(z, y, standardize=F, alpha=1,family="binomial", lambda=seq(0.00005,0.03,0.001))
cv.clrlasso <- cv.glmnet(z, y, standardize=T, alpha=1,family="binomial", lambda=seq(0.00005,0.04,0.001))

plot(cv.clrlasso)
print(cv.clrlasso)

#
HFHS.lambda_clr = 0.03

results_clrlasso <- glmnet(x = HFHS.clrx, y = HFHS.y.num,  family = 'binomial',
                                        lambda = HFHS.lambda_clr)

HFHS.results_clrlasso <- glmnet_wrapper(results_clrlasso, HFHS.clrx)

HFHS.results_clrlasso$numVarSelect

HFHS.results_clrlasso$varIndex

HFHS.results_clrlasso$varSelect

HFHS.results_clrlasso$explained_deviance_proportion

selected<- HFHS.results_clrlasso$varSelect

HFHS.tax_clrlasso <- taxonomy_tab[which(rownames(taxonomy_tab)%in% selected),]

HFHS.tax_clrlasso[ ,2:6]


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


  x <- X_day1

  y <- y_day1

  lambdaRange_codalasso(y,x, seq(0,2,0.01))

  HFHS.results_codalasso <- coda_logistic_lasso(y,x,lambda=1.47)  # 7 variables, %Dev 1

  HFHS.results_codalasso$`number of iterations`
  HFHS.results_codalasso$`number of selected variables`
  HFHS.results_codalasso$`indices of selected variables`
  HFHS.results_codalasso$`name of selected variables`
  HFHS.results_codalasso$`beta non-zero coefficients`
  HFHS.results_codalasso$`proportion of explained deviance`
  HFHS.results_codalasso$betas

  selected<- HFHS.results_codalasso$`name of taxa with non-zero coeff`
  HFHS.tax_codalasso <- taxonomy_tab[which(rownames(taxonomy_tab)%in% selected),]

  HFHS.tax_codalasso[ ,2:6]


  selected <- results_codalasso[[4]]

  column_selected <- results_codalasso[[3]][-1]-1

  tax <- taxonomy_tab[which(rownames(taxonomy_tab)%in% selected),]
  selected_codalasso <- results_codalasso[[4]][-1]
  #
  columns_selected_codalasso <- results_codalasso[[3]][-1]

  results_codalasso$`number of iterations`
  results_codalasso$`number of selected taxa`
  results_codalasso$`indices of taxa with non-zero coeff`
  results_codalasso$`name of taxa with non-zero coeff`
  results_codalasso$`beta non-zero coefficients`
  results_codalasso$`proportion of explained deviance`
  results_codalasso$betas

  write.csv(data.frame(column_selected,selected_codalasso),"results_codalasso_BEMEday1.csv")

}


#######################################################################################################
#
# Concordance of variables selected by the three methods
#
#######################################################################################################


d_selbal <- read.csv("results_selbal_BEMEday1.csv", sep=",", header = T)
d_clrlasso <- read.csv("results_clrlasso_BEMEday1.csv", sep=",", header = T)
d_codalasso <- read.csv("results_codalasso_BEMEday1.csv", sep=",", header = T)

taxa_selected <- list(d_clrlasso[,2], d_codalasso[,2], d_selbal[,2])
taxa.id_selected <- list(d_codalasso[,3], d_clrlasso[,3], d_selbal[,3])


venn.plot <- venn.diagram(taxa.id_selected , NULL, fill=c("blue", "magenta", "lightblue"),
                          alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4,
                          category.names=c("coda+lasso", "clr-lasso", "selbal"), main="Concordance of selected taxa for BEME day1")
dev.off() #to clean plots window
grid.draw(venn.plot)

