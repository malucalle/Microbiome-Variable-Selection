#######################################################################
#                                                                     #
#      functions_microbiome_variable_selection_tutorial:              #
#                                                                     #
#######################################################################

#######################################################################
#                   source coda_lasso function
#######################################################################

# copy the repository from https://github.com/UVic-omics/CoDA-Penalized-Regression
system('git clone https://github.com/UVic-omics/CoDA-Penalized-Regression')

# fetch the last modified repository from https://github.com/UVic-omics/CoDA-Penalized-Regression
# when you have already git clone the repository
# system('git pull https://github.com/UVic-omics/CoDA-Penalized-Regression')

# source the required functions
source(file = 'CoDA-Penalized-Regression/R/functions_coda_penalized_regression.R')


#######################################################################
#                         coda_lasso_wrapper
#######################################################################

coda_lasso_wrapper  <-  function(result, X){
  
  # result: result from coda_logistic_lasso()
  # X: matrix of k covariates (positive values, taxa abundances in counts, 
  #    proportions, intensities, ...), matrix of dimension n by k.
  
  coefficientsSelect <- result$`beta non-zero coefficients`[-1]
  names(coefficientsSelect) <- result$`name of selected variables`
  coefficientsSelect <- coefficientsSelect[order(abs(coefficientsSelect), decreasing = T)]
  
  coefficients <- result$betas[-1]
  names(coefficients) <- colnames(X)
  
  varSelect = names(coefficientsSelect)
  
  varIndex <- match(varSelect, colnames(X))
  
  # choose the desired output from 'result'
  out = list(
    call = match.call(),
    iter = result$`number of iterations`,
    numVarSelect = result$`number of selected variables`,
    varSelect = varSelect,
    varIndex = varIndex,
    coefficientsSelect = coefficientsSelect,
    posCoefSelect = coefficientsSelect[which(coefficientsSelect > 0)],
    negCoefSelect = coefficientsSelect[which(coefficientsSelect < 0)],
    coefficients = coefficients,
    explained_deviance_proportion = result$`proportion of explained deviance`
  )
  
  return(invisible(out))
  
}


#######################################################################
#                           glmnet_wrapper 
#######################################################################

glmnet_wrapper  <-  function(result, X){
  
  # result: result from glmnet()
  # X: clr transformed matrix of k covariates, matrix of dimension n by k.
  
  coefficientsSelect <- result$beta[which(result$beta[,1] != 0),]
  coefficientsSelect <- coefficientsSelect[order(abs(coefficientsSelect), decreasing = T)]
  
  varSelect = names(coefficientsSelect)
  
  varIndex <- match(varSelect, colnames(X))
  
  # choose the desired output from 'result'
  out = list(
    call = match.call(),
    numVarSelect = result$df,
    varSelect = varSelect,
    varIndex = varIndex,
    coefficientsSelect = coefficientsSelect,
    posCoefSelect = coefficientsSelect[which(coefficientsSelect > 0)],
    negCoefSelect = coefficientsSelect[which(coefficientsSelect < 0)],
    coefficients = result$beta,
    explained_deviance_proportion = result$dev.ratio
  )
  
  return(invisible(out))
  
}

#######################################################################
#                           selbal_wrapper
#######################################################################

selbal_wrapper  <-  function(result, X){
  
  # result: result from selbal()
  # X: a matrix object with the information of variables (columns) for each sample (rows).
  
  # choose the desired output from 'result'
  out = list(
    call = match.call(),
    finalBal = result[[1]],
    posVarSelect = result[[2]],
    negVarSelect = result[[3]],
    numVarSelect = length(result[[4]]),
    varSelect = result[[4]],
    varIndex = match(result[[4]], colnames(X)),
    steps = length(result[[5]]),
    ACC_eachStep = result[[5]]
  )
  
  return(invisible(out))
  
}

#######################################################################
#                       Additional functions


#######################################################################
#                        Transparent color
#######################################################################
color <- c('#388ECC', '#F68B33', '#C2C2C2', '#009E73', '#CC79A7', '#F0E442', '#6073B1', 'black', 
           '#D55E00', '#999999', '#E69F00', '#56B4E9', '#994F00', '#40B0A6', '#5D3A9B', '#E1BE6A', 
           '#005AB5')

t_col <- function(color, percent = 50, name = NULL) {
  #	color = color name
  #	percent = % transparency
  #	name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255, alpha = (100 - percent) * 255 / 100,
               names = name)
  ## Save the color
  invisible(t.col)
}

#######################################################################
#                           selbal plot
#######################################################################

selbal_like_plot <- function(pos.names, neg.names, Y, X = NULL, selbal = FALSE, FINAL.BAL = NULL, OTU = F, taxa = NULL){
  # pos.names: the name of variables selected with positive coef # selbal: numerator 
  # neg.names: the name of variables selected with negative coef # selbal: denominator
  # Y: the response variable, should be binomial
  # X: a matrix object with the information of variables (columns) for each sample (rows).
  # selbal: logical value determining if the results are from selbal function
  # FINAL.BAL: if selbal = TRUE, the FINAL.BAL should not be NULL.
  # OTU: logical value determining if the input positive and negative names are OTUs
  # taxa: rownames are OTUs. if OTU = T, taxa should not be NULL
  
  require(ggplot2)
  require(gridExtra)
  
  #-----------------------------------------#
  # FIRST: The names of included variables
  #-----------------------------------------#
  
  if(OTU == TRUE){
    if(is.null(taxa)){
      stop("'taxa' should not be NULL")
    }
    
    pos.names.taxa = c()
    for(e in pos.names){
      index = ncol(taxa)
      taxa.name = ''
      while(taxa.name == ''){
        taxa.name <- taxa[e,index]
        index <- index - 1		
      }
      pos.names.taxa = c(pos.names.taxa, taxa.name)
    }
    pos.names.print <- paste0(pos.names, ': ', pos.names.taxa)
    
    
    neg.names.taxa = c()
    for(e in neg.names){
      index = ncol(taxa)
      taxa.name = ''
      while(taxa.name == ''){
        taxa.name <- taxa[e,index]
        index <- index - 1		
      }
      neg.names.taxa = c(neg.names.taxa, taxa.name)
    }
    neg.names.print <- paste0(neg.names, ': ', neg.names.taxa)
    
  }else{
    pos.names.print <- pos.names
    neg.names.print <- neg.names
  }
  
  
  ###################
  
  pos.var <- c(pos.names.print, 'A')
  neg.var <- c(neg.names.print, 'B')
  
  # Parameter to specify the limits for writting
  yl <- max(length(pos.var), length(neg.var))*1.05
  escal <- 3
  bot <- 5*escal*yl/100;
  
  # Empty plot with text
  ndiv = max(3,floor(yl))+1;
  lineh <- 0 
  df <- data.frame()
  Imp.table <- ggplot(df) + xlim(0, 100) + ylim(-bot, ndiv) + theme_void() +
    geom_segment(aes(x = 10, y = lineh, xend = 90, yend = lineh),color='brown3', size=1.3) +
    geom_segment(aes(x = 50-escal, y = lineh-bot, xend = 50+escal, yend = lineh-bot),color='brown3', size=1) +
    geom_segment(aes(x = 50+escal, y = lineh-bot, xend = 50+escal, yend = lineh),color='brown3', size=1) +
    geom_segment(aes(x = 50-escal, y = lineh-bot, xend = 50-escal, yend = lineh),color='brown3', size=1) +
    annotate('text',
             x = 75,
             y = seq(bot, floor(yl), length.out = ndiv),
             label = c(pos.var,rep('',ndiv-length(pos.var))),
             colour = c(rep('black',length(pos.var)-1),rep('brown3',ndiv-length(pos.var)+1)),
             fontface = 2) +
    annotate('text',
             x = 25,
             y = seq(bot, floor(yl), length.out = ndiv),
             label = c(neg.var,rep('',ndiv-length(neg.var))),
             colour = c(rep('black',length(neg.var)-1),rep('brown3',ndiv-length(neg.var)+1)),
             fontface = 2)
  
  
  #-----------------------------------------#
  # SECOND: The representation of the plots
  #-----------------------------------------#
  
  if(selbal == TRUE){
    if(is.null(FINAL.BAL)){
      stop("'FINAL.BAL' should not be NULL")
    }
    Final.score <- FINAL.BAL
  }else{
    if(is.null(X)){
      stop("'X' should not be NULL")
    }
    data.log <- log(X)
    pos.index <- match(pos.names, colnames(X))
    neg.index <- match(neg.names, colnames(X))
    FINAL.mean.pos <- apply(data.log, 1, FUN = function(x, pos.index){mean(x[pos.index])}, pos.index)
    FINAL.mean.neg <- apply(data.log, 1, FUN = function(x, neg.index){mean(x[neg.index])}, neg.index)
    Final.score <- FINAL.mean.pos - FINAL.mean.neg
  }
  
  Y = as.factor(Y)
  
  U <- data.frame(Y, Final.score)
  colnames(U)[ncol(U)] <- 'V1'
  
  # The plot depending of the class of the response variable
  #---------------------------------------------------#
  # The composition of the final plot
  #---------------------------------------------------#
  # BOXPLOT
  col = c('#388ECC', '#F68B33')
  if(selbal == TRUE){
    BoxP <-  ggplot(U, aes(x = Y, y = V1, fill = Y)) +
      geom_boxplot(color = 'black', size = 1) +
      scale_fill_manual(values = col) +
      theme_bw() +
      ylab('Balance') +
      xlab('') +
      theme(legend.position = 'none') + 
      coord_flip()
  }else{
    BoxP <-  ggplot(U, aes(x = Y, y = V1, fill = Y)) +
      geom_boxplot(color = 'black', size = 1) +
      scale_fill_manual(values = col) +
      theme_bw() +
      ylab('Log mean difference between NEG and POS variables') +
      xlab('') +
      theme(legend.position = 'none') +
      coord_flip()
    }
  # Density plot for the balance
  ydensity <- ggplot(U, aes(V1, fill = Y)) +
    geom_density(alpha = .5, size = 1) +
    scale_fill_manual(values = col) +
    theme_bw() + xlab('') + ylab('') +
    theme(legend.position = 'none') 
  
  grid.arrange(Imp.table, BoxP, ydensity, ncol = 1)

}


#######################################################################
#                           plotLoadings
#######################################################################

plotcoefficients <- function(coef, data, Y, method = 'mean', contrib.method = 'max', title = 'Coefficients', 
                             cex.x = 0.8, cex.y = 0.6, cex.legend = 0.8, OTU = F, taxa = NULL){
  # coef: the coefficients that are going to be plotted, should have variable names
  # data: the data of variables with non-zero coefficients, samples in the row, and variables in the column
  # Y: sample group info
  # method: the way we calculate, can be median or mean, the default is mean
  # contrib.method: the way we use to decide the contributing group, can be min or max, the default is max
  # OTU: logical value determining if the input coef names are OTUs
  # taxa: rownames are OTUs. if OTU = T, taxa should not be NULL
  
  color <- c('#388ECC', '#F68B33', '#C2C2C2', '#009E73', '#CC79A7', '#F0E442', '#6073B1', 'black', 
             '#D55E00', '#999999', '#E69F00', '#56B4E9', '#994F00', '#40B0A6', '#5D3A9B', '#E1BE6A', 
             '#005AB5')
  
  if(is.null(names(coef))){
    stop('The coefficients should have variable names')
  }
  
  if(OTU == TRUE){
    if(is.null(taxa)){
      stop("'taxa' should not be NULL")
    }
    
    coef.taxa = c()
    for(e in names(coef)){
      index = ncol(taxa)
      taxa.name = ''
      while(taxa.name == ''){
        taxa.name <- taxa[e,index]
        index <- index - 1		
      }
      coef.taxa = c(coef.taxa, taxa.name)
    }
    coef.name.taxa <- paste0(names(coef), ': ', coef.taxa)
  }else{
    coef.name.taxa <- names(coef)
  }
  
  coef <- coef[order(abs(coef), decreasing = T)]
  data <- data[ ,as.character(names(coef))]
  
  Y = as.factor(Y)
  which.contrib = data.frame(matrix(FALSE, ncol = nlevels(Y) + 1, nrow = length(coef), 
                                    dimnames = list(names(coef), c(paste0('Contrib.', levels(Y)), 'GroupContrib'))))
  
  method.group <- list()
  for(k in 1:ncol(data))
  {
    method.group[[k]] = tapply(data[, k], Y, method, na.rm = TRUE) # method is either mean or median
    # determine which group has the highest mean/median
    which.contrib[k, 1:nlevels(Y)] = (method.group[[k]]) == get(contrib.method)(method.group[[k]]) # contrib.method is either min or max
  }
  
  col.ties = 'white'
  
  which.contrib$color = apply(which.contrib, 1, function(x)
  {
    if (length(which(x)) > 1)
    {
      return(col.ties)
    } else { # otherwise we use legend color provided
      return(color[1:nlevels(Y)][which(x)])
    }
  })
  
  
  which.contrib$GroupContrib = apply(which.contrib[ ,1:(nlevels(Y))], 1, function(x)
  {
    if (length(which(x)) > 1)
    {
      return('tie')
    } else {
      return(levels(Y)[which(x)])
    }
  })
  
  par(mar=c(5.1, 8.1, 4.1, 8.1), xpd=TRUE)
  p <- barplot(coef, horiz = T, col = which.contrib$color, yaxt = 'n',main = title, 
               cex.axis = cex.x, xlim = c(round((min(coef) * 1.2),3), round((max(coef) * 1.3),3)))  
  axis(side = 2, at = p, labels = coef.name.taxa, las = 2, tick = F, cex.axis = cex.y)
  legend('right', inset = c(-0.3,0), legend = levels(Y), pch = 16, title = 'Outcome', col = color[1:nlevels(Y)], cex = cex.legend)
  par(mar = c(5.1,4.1,4.1,2.1))
  
  results <- data.frame(coef, which.contrib)
  return(invisible(results))
}


#######################################################################
#                         Trajectory plot
#######################################################################

TRAJ_plot <- function(selectVar_coef_method1, selectVar_coef_method2, selectMethods, OTU = F, taxa = NULL){
  # selectVar_coef_method1: the vector of coefs with variable names as the vector names selected by the first method
  # selectVar_coef_method2: the vector of coefs with variable names as the vector names selected by the second method
  # selectMethods: methods used to select variables
  # OTU: logical value determining if the input coef names are OTUs
  # taxa: rownames are OTUs. if OTU = T, taxa should not be NULL
  
  require(ggplot2)
  require(ggforce)
  
  # order the absolute coefs of the first method
  selectVar_coef_method1 <- selectVar_coef_method1[order(abs(selectVar_coef_method1), decreasing = T)]
  
  selectVar_method1 <- names(selectVar_coef_method1)
  selectVar_method2 <- names(selectVar_coef_method2)
  unionVar <- union(selectVar_method1, selectVar_method2)
  
  plot.matrix <- matrix(0, nrow = length(unionVar), ncol = 2, dimnames = list(unionVar, selectMethods))
  plot.matrix[match(selectVar_method1, unionVar),1] <- selectVar_coef_method1
  plot.matrix[match(selectVar_method2, unionVar),2] <- selectVar_coef_method2
  
  # abs
  plot.matrix.abs <- abs(plot.matrix)
  
  # proportion
  plot.matrix.prop <- apply(plot.matrix.abs, 2, function(x){x / sum(x)})
  
  # rank
  plot.matrix.rank <- apply(plot.matrix.abs, 2, rank, ties.method = 'max')
  
  # increase the diff between variables
  plot.matrix.rank <- plot.matrix.rank * 5
  
  matrix.names <- paste0(rep(rownames(plot.matrix.rank), each = 2), 1:2)
  plot.matrix.rank_prop <- data.frame(matrix(0, nrow = 2 * nrow(plot.matrix.rank), ncol = 2, 
                                             dimnames = list(matrix.names, selectMethods)))
  
  for(i in 1:nrow(plot.matrix.rank)){
    plot.matrix.rank_prop[(2 * i-1), ] <- plot.matrix.rank[i, ] + plot.matrix.prop[i, ] * 5 # * 10 / 2
    plot.matrix.rank_prop[2 * i, ] <- plot.matrix.rank[i, ] - plot.matrix.prop[i, ] * 5
  }
  
  # same dist between each two sequential variables
  plot.matrix.rank_prop.nospace <- plot.matrix.rank_prop
  plot.matrix.rank_prop.nospace$methodVar2 <- plot.matrix.rank_prop.nospace$methodVar1 <- rep(rownames(plot.matrix.rank), each = 2)
  method1.null <- setdiff(selectVar_method2, selectVar_method1)
  method2.null <- setdiff(selectVar_method1, selectVar_method2)
  
  plot.matrix.rank_prop.nospace$methodVar1[c(match(method1.null, plot.matrix.rank_prop.nospace$methodVar1),
                                             match(method1.null, plot.matrix.rank_prop.nospace$methodVar1) + 1)] = 'Not selected'
  
  plot.matrix.rank_prop.nospace$methodVar2[c(match(method2.null, plot.matrix.rank_prop.nospace$methodVar2),
                                             match(method2.null, plot.matrix.rank_prop.nospace$methodVar2) + 1)] = 'Not selected'
  
  for(j in 1:2){
    data <- plot.matrix.rank_prop.nospace[ ,j]
    names(data) <- rownames(plot.matrix.rank_prop.nospace)
    data.order <- sort(data, decreasing = T)
    
    p <- length(data.order)
    data.res <- data.order
    for(i in 2:nrow(plot.matrix.rank)){
      data.res[(2 * i-1):p] <- data.res[(2 * i-1):p] + data.res[2 * (i-1)] - data.res[(2 * i-1)] - 0.2
    }
    plot.matrix.rank_prop.nospace[names(data.res),j] = data.res
  }
  
  method1.index.null <- which(plot.matrix.rank_prop.nospace$methodVar1 == 'Not selected')
  plot.matrix.rank_prop.nospace[method1.index.null,1] = plot.matrix.rank_prop.nospace[method1.index.null[1],1]
  
  method2.index.null <- which(plot.matrix.rank_prop.nospace$methodVar2 == 'Not selected')
  plot.matrix.rank_prop.nospace[method2.index.null,2] = plot.matrix.rank_prop.nospace[method2.index.null[1],2]
  
  
  method1.loc <- tapply(plot.matrix.rank_prop.nospace[,1], plot.matrix.rank_prop.nospace$methodVar1, mean)
  
  method2.loc <- tapply(plot.matrix.rank_prop.nospace[,2], plot.matrix.rank_prop.nospace$methodVar2, mean)
  
  # plot
  color <- c('#388ECC', '#F68B33', '#C2C2C2', '#009E73', '#CC79A7', '#F0E442', '#6073B1', 'black', 
             '#D55E00', '#999999', '#E69F00', '#56B4E9', '#994F00', '#40B0A6', '#5D3A9B', '#E1BE6A', 
             '#005AB5')
  
  more_color = c()
  percent = 30
  for(k in 1:2){
    for(c in color){
      new.color <- t_col(color = c, percent = percent)
      more_color <- c(more_color, new.color)
    }
    percent = percent + 30
  }
  
  color <- c(color, more_color) # in case the color is run out of
  
  y = c()
  for(i in 1:nrow(plot.matrix.rank)){
    y = c(y, unlist(plot.matrix.rank_prop.nospace[(2 * i-1),1:2]),
          unlist(plot.matrix.rank_prop.nospace[2 * i,2:1]))
    
  }
  data <- data.frame(x = rep(c(1, 6, 6, 1), nrow(plot.matrix.rank)), 
                     y = y, group = rep(1:nrow(plot.matrix.rank), each = 4))
  
  
  
  if(OTU == TRUE){
    if(is.null(taxa)){
      stop("'taxa' should not be NULL")
    }
    
    method1.taxa = c()
    for(e in names(method1.loc)){
      if(e != 'Not selected'){
        index = ncol(taxa)
        taxa.name = ''
        while(taxa.name == ''){
          taxa.name <- taxa[e,index]
          index <- index - 1		
        }
        taxa.name <- paste0(e, ': ', taxa.name)
        method1.taxa = c(method1.taxa, taxa.name)
      }else{
        taxa.name = e
        method1.taxa = c(method1.taxa, taxa.name)
      }}
    
    method2.taxa = c()
    for(e in names(method2.loc)){
      if(e != 'Not selected'){
        index = ncol(taxa)
        taxa.name = ''
        while(taxa.name == ''){
          taxa.name <- taxa[e,index]
          index <- index - 1		
        }
        taxa.name <- paste0(e, ': ', taxa.name)
        method2.taxa = c(method2.taxa, taxa.name)
      }else{
        taxa.name = e
        method2.taxa = c(method2.taxa, taxa.name)
      }}
    
  }else{
    method1.taxa <- names(method1.loc)
    method2.taxa <- names(method2.loc)
  }      
  
  
  
  
  
  ggplot(data) + geom_diagonal_wide(aes(x, y, group = group, 
                                        fill = as.factor(group)), alpha = 0.8) + 
    theme_bw() + theme(panel.border = element_blank(), panel.grid = element_blank(), 
                       axis.title = element_blank(), axis.ticks.x = element_blank()) + 
    scale_fill_manual(values = color) + 
    scale_y_continuous(breaks = method1.loc, labels = method1.taxa, 
                       sec.axis = sec_axis(~ . , breaks = method2.loc, 
                                           labels = method2.taxa)) + 
    theme(legend.position = 'none') + scale_x_continuous(breaks = c(1,6), 
                                                         labels = selectMethods)
  
}


#######################################################################
#                         GraPhlAn
#######################################################################

graphlan_annot_generation <- function(taxa_list, save_folder){
  # taxa_list: a list of taxa info selected by different methods. The list should have names for each element
  # save_folder: where to save your generated txt files
  
  # checking
  if(is.null(names(taxa_list))){
    stop('The taxa list should have names for each element.')
  }
  
  #---------------------------
  # generate taxa.txt file
  #---------------------------
  for(i in 1:length(taxa_list)){
    taxa_list[[i]] <- as.data.frame(taxa_list[[i]])
    taxa_list[[i]]$OTUs <- rownames(taxa_list[[i]])
    rownames(taxa_list[[i]]) <- NULL
  }
  
  tax.df <- do.call(rbind, taxa_list)
  tax.df.unique <- tax.df[match(unique(tax.df$OTUs),tax.df$OTUs),]
  
  # make each column to be character
  for(i in 1:ncol(tax.df.unique)){
    tax.df.unique[,i] = as.character(tax.df.unique[,i])
  }
  
  # save the taxa.txt file
  write.table(tax.df.unique, file = paste0(save_folder, '/taxa.txt'), sep = '.', 
              quote = FALSE, col.names = F, row.names = F)
  
  #----------------------------
  # generate annot.txt file
  #----------------------------
  color.vector <- c('#388ECC', '#F68B33', '#C2C2C2', '#009E73', '#CC79A7', '#F0E442', '#6073B1', 'black', 
                    '#D55E00', '#999999', '#E69F00', '#56B4E9', '#994F00', '#40B0A6', '#5D3A9B', '#E1BE6A', 
                    '#005AB5')
  
  ################################################################
  # We characterise each taxa with background color
  # except Kingdom (1st column) and OTUs (last column) 
  branch.each_taxa <- c()
  for(i in 2:(ncol(tax.df.unique) - 1)){
    branch.each_taxa <- c(branch.each_taxa, unique(tax.df.unique[ ,i]))
  }
  branch.each_taxa <- branch.each_taxa[branch.each_taxa != '']
  
  # The background color is blocks of colors according to Phylum 
  Phylum.color.list <- list()
  for(i in 1:length(unique(tax.df.unique[ ,2]))){
    phy.inx <- which(tax.df.unique[ ,2] == tax.df.unique[ ,2][i])
    Phylum.color <- c()
    for(j in 2:(ncol(tax.df.unique) - 1)){
      Phylum.color <- c(Phylum.color, unique(tax.df.unique[phy.inx,j]))
    }
    Phylum.color <- Phylum.color[Phylum.color != '']
    Phylum.color.list[[i]] <- Phylum.color 
  }
  names(Phylum.color.list) = unique(tax.df.unique[,2])
  
  # build the background color matrix
  background.color <- data.frame(matrix(NA, nrow = 3 * length(branch.each_taxa), ncol = 3))
  background.color[ ,1] <- rep(branch.each_taxa, each = 3)
  background.color[ ,2] <- c('annotation', 'annotation_background_color', 'annotation_font_size')
  background.color[ ,3] <- ifelse(background.color[,2] == 'annotation', background.color[ ,1], 5)
  color.b.index <- which(background.color[,2] == 'annotation_background_color')
  background.color.vector <- c(NA, length = length(branch.each_taxa))
  for(i in 1:length(Phylum.color.list)){
    index <- match(Phylum.color.list[[i]], branch.each_taxa)
    background.color.vector[index] <- color.vector[i + length(taxa_list)]
  }
  
  background.color[color.b.index,3] <- background.color.vector
  
  
  # make it able to merge with the later matrix
  background.color.four_column <- cbind(background.color, NA)
  colnames(background.color.four_column) <- paste0('c', 1:4)
  
  ####################################
  # build the rings outside
  # each ring represents one method
  outside.ring.global <- data.frame(matrix(NA, nrow = 2, ncol = 4))
  outside.ring.global[1, ] <- c('total_plotted_degrees', 340, NA, NA)
  outside.ring.global[2, ] <- c('start_rotation', 270, NA, NA)
  
  n_rings <- length(taxa_list) + 1
  outside.shape.all <- data.frame(matrix(NA, nrow = n_rings * 4, ncol = 4))
  outside.shape.all[ ,1] <- c('ring_internal_separator_thickness', 'ring_label',
                              'ring_label_color', 'ring_label_font_size')
  outside.shape.all[ ,2] <- rep(1:n_rings, each = 4)
  outside.shape.all[ ,3] <- ifelse(outside.shape.all[ ,1] == 'ring_internal_separator_thickness', 0.5, 7)
  ring_label_index <- which(outside.shape.all[ ,1] == 'ring_label')
  outside.shape.all[ring_label_index,3] <- c(names(taxa_list), 'outside')
  ring_label_col.index <- which(outside.shape.all[ ,1] == 'ring_label_color')
  outside.shape.all[ring_label_col.index,3] <- c(color.vector[1:length(taxa_list)], 'w')
  outside.shape.all[ ,4] <- NA
  
  tax.name.list <- read.table(file = paste0(save_folder, '/taxa.txt'))
  outside.ring.shape <- data.frame(matrix(NA, nrow = n_rings*nrow(tax.name.list), ncol = 4))
  outside.ring.shape[,1] <- tax.name.list[ ,1]
  outside.ring.shape[,2] <- 'ring_color'
  outside.ring.shape[,3] <- rep(1:n_rings, each = nrow(tax.name.list))
  outside.ring.shape[,4] <- rep(c(color.vector[1:length(taxa_list)], 'w'), each =  nrow(tax.name.list))
  
  outside.ring.alpha <- outside.ring.shape
  outside.ring.alpha[ ,2] <- 'ring_alpha'
  outside.ring.alpha[ ,4] <- 0
  alpha.index <- c()
  for(i in 1:length(taxa_list)){
    index <- match(taxa_list[[i]]$OTUs, tax.df.unique$OTUs) + nrow(tax.name.list) * (i-1)
    alpha.index <- c(alpha.index, index) 
  }
  outside.ring.alpha[alpha.index,4] <- 1
  
  
  # The OTUs that are selected only by one method will be enlarged
  OTU.size <- data.frame(matrix(NA, nrow = nrow(tax.df.unique), ncol = 4))
  OTU.size[ ,1] <- tax.name.list[ ,1]
  OTU.size[ ,2] <- 'clade_marker_size'
  OTU.sizes <- c()
  for(i in 1:nrow(tax.df.unique)){
    OTU.sizes[i] <- length(which(tax.df$OTUs == tax.df.unique$OTUs[i]))
  }
  OTU.size[ ,3] <- OTU.sizes * 25
  OTU.size[ ,4] <- NA
  
  
  label_index <- which(OTU.sizes == min(OTU.sizes))
  OTU.clade <- data.frame(matrix(NA, nrow = length(label_index), ncol = 4))
  OTU.clade[ ,1] <- tax.df.unique$OTUs[label_index]
  OTU.clade[ ,2] <- 'clade_marker_color'
  OTU.clade[ ,3] <- '#555555'
  OTU.clade[ ,4] <- NA
  
  four_column <- rbind(outside.ring.global, outside.shape.all, outside.ring.shape, 
                       outside.ring.alpha, OTU.size, OTU.clade)
  colnames(four_column) <- paste0('c', 1:4)
  annot_file <- rbind(background.color.four_column, four_column)
  
  # save the annot file
  write.table(annot_file, file = paste0(save_folder, '/annot_all.txt'), sep = '\t', quote = FALSE, col.names = F, row.names = F, na = '')
  
}





