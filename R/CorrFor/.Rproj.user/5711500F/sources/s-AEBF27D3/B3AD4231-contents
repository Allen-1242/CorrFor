#Loading libraries 
library(dplyr)
library(tidyr)
library(ranger)
library(lsr)
library(igraph)
library(data.table)
library(Matrix)
library(purrr)


#Loading the example dataset 
df <- read.csv("C:\\Users\\sunny\\OneDrive\\Documents\\GitHub\\CorFor\\default of credit card clients.csv")

#Generalized correlation function
cor2 = function(df)
{
  #If both are numeric
  if(class(df[[pos_1]]) %in% c("integar", "numeric") && 
     class(df[[pos_2]]) %in% c("integar", "numeric")){
       r <- stats::cor(df[[pos_1]], df[[pos_2]], use = 'pairwise.complete.obs')
    }

  #if one is numeric and the other is a factor or character 
  if(class(df[[pos_1]]) %in% c("integar", "numeric") && 
    class(df[[pos_2]]) %in% c("factor", "character")){
      r <- sqrt(summary(stats::lm(df[[pos_1]] ~ as.factor(df[[pos_2]])))[['r.squared']])
    }
  
  
  if(class(df[[pos_2]]) %in% c("integar", "numeric") && 
     class(df[[pos_1]]) %in% c("factor", "character")){
    r <- sqrt(summary(stats::lm(df[[pos_2]] ~ as.factor(df[[pos_1]])))[['r.squared']])
  }
  
  #If both are factors/characters
  if(class(df[[pos_1]]) %in% c("factor", "character") && 
     class(df[[pos_2]]) %in% c("factor", "character")){
    r <- lsr::cramersV(df[[pos_1]], df[[pos_2]], simulate.p.value = TRUE)
  }
}

#Main function
Main_func = function(Data, Y_var, Focus_variables, corr_cutoff, RF_info_cutoff, plot = FALSE, fast_calculation = FALSE)
{
  Data <- as.data.table(data)
  
  cor_matrix <- cor2(Data[, Y_var, with = FALSE])
  
  pairs_mat <- cor_matrix %>%
                as.data.frame() %>%
                dplyr::mutate(var1 = rownames(.)) %>%
                tidyr::gather(var2, value, -var1) %>%
                arrange(desc(value))
  
  #Sub-setting values only above 0.9
  pairs_mat = pairs_mat[which(abs(pairs_mat$value) >= corr_cutoff & (pairs_mat$var1 != pairs_mat$var2)),]
  
  
  #Creating clusters based upon graph theory 
  list1 = list()
  
  if(dim(pairs_mat)[1] != 0)
  {
    for(j in 1:nrom(pairs_mat))
    {
      list1[[j]] = c(pairs_mat[j,1], pairs_mat[j,2])
      list1[[j]] = sort(list1[[j]])
    }
    
    list2 = unique(list1)
    i = rep(1:length(list1), lengths(list1))
    j = factor(unlist(list1))
    tab = sparseMatrix(i = i , j = as.integer(j), x = TRUE, dimnames= list(NULL, levels(j)))
    connects = tcrossprod(tab, boolArith = TRUE)
    group = clusters(graph_from_adjacency_matrix(as(connects, "lsCMatrix")))$membership
    var_groups <- tapply(list1, group, function(x) sort(unique(unlist(x))))
    var_groups <- as.list(var_groups)
    
    #Condition to extract the focus variables
    for(i in 1:length(var_groups))
    {
      if(any(var_groups[[i]] %in% Focus_variables))
      {
        if(sum(which(var_groups[[i]] %in% Focus_variables >= 2)))
        {
          Intersection <- dplyr::intersect(Focus_variables, var_groups[[i]])
          var_groups[[i]] <- Intersection[1]
          
          next
        }
        
        var_groups[[i]] <- dplyr::intersect(Focus_variables, var_groups[[i]])
        
        
      }
    }
    
    #Getting every combination of variables possible
    if(fast_calculation == TRUE)
    {
      result <- purrr::map(var_groups, 1)
      result <- as.data.frame(t(unlist(result)))
    }else
    {
      result <- expand.grid(var_groups)
    }
    
    RF_list <- list()
    Data_nocor <- Data[, dplyr::union(unlist(var_groups), unlist(Focus_variables)), with = FALSE]
    
    ##Start of the RF
    for(i in 1:nrow(result))
    {
      Data_temp <- cbind(Data_nocor, Data[, unlist(result[i,]), with = FALSE])
      
      Rf <- ranger(as.formula(paste(paste(Y_var, '~'), paste(colnames(Data_temp), collapse = "+"))), data = Data_temp, mtry = ncol(Data_temp/3), importance = 'permutation')
      Rf_2 <- data.frame(Rf$variable.importance)
      RF_list[[i]] <- Rf_2
    }
  }else
  {
    RF_list <- list()
    Rf <- ranger(as.formula(paste(paste(Y_var, '~'), paste(colnames(Data_temp), collapse = "+"))), data = Data_temp, mtry = ncol(Data_temp/3), importance = 'permutation')
    Rf_2 <- data.frame(Rf$variable.importance)
    RF_list[[i]] <- Rf_2
  }
  
  #Fast Aggregation across multiple frames
  l <- lapply(RF_list, function(x) {x$Rowname <- row.names(x) ; x})
  Res <- Reduce(function(...) merge(..., by = 'RowName', all = TRUE), l)
  
  Rf_2 <- dcast(melt(setDT(Res), "RowName"),
                RowName ~ sub("\\..*","",variable),
                mean,
                na.rm = TRUE,
                value.var = "value")

  #Taking the best variable from each group
  if(dim(pairs_mat)[1] != 0)
  {
    for(bv in 1:nrow(var_groups))
    {
      comp <- var_groups[bv]
      comp <- unlist(strsplit(comp[[1]]), '"')# Change is needed?
      temp <- Rf_2[which(Rf_2$Rowname %in% comp)]
      keep_var <- Rf_2$RowName[Rf_2$Rf == max(temp$Rf)]
      rem_var <- comp[which(comp != keep_var)]
      
      #Dropping all values not needed
      Rf_2 <- Rf_2[!(Rf_2$RowName %in% unlist(rem_var))]
    }
    
    #Fitting the RF without the correlated variables 
    temp_var <- c(Rf_2$RowName, 'Y_var')
    Data_temp <- Data[,c(...temp_var)] 
    
    Rf_list = list()
    
    Rf <- ranger(as.formula(paste(paste(Y_var, '~'), paste(colnames(Data_temp), collapse = "+"))), data = Data_temp, mtry = ncol(Data_temp/3), importance = 'permutation')
    Rf_2 <- data.frame(Rf$variable.importance)
    RF_list[[1]] <- Rf_2             
    
    l <- lapply(RF_list, function(x) {x$Rowname <- row.names(x) ; x})
    Res <- Reduce(function(...) merge(..., by = 'RowName', all = TRUE), l)
    
    Rf_2 <- dcast(melt(setDT(Res), "RowName"),
                  RowName ~ sub("\\..*","",variable),
                  mean,
                  na.rm = TRUE,
                  value.var = "value")
    }

  Rf_2 <- na.omit(Rf_2)
  Rf_2 <- Rf_2[which(Rf_2$Rf >= 0)]
  Rf_2 <- Rf_2[order(-Rf)]
  
  Rf_2 <- Rf_2$Rf/sum(Rf_2$Rf)
  Rf_2 <- cumsum(Rf_2$Rf)
  temp <- which(Rf_2$Rf > Rf_info_cutoff)[1]
  
  if(plot == TRUE)
  {
    print('lo')
  }
  
}