# Summary functions for quick normalization checks
library(Polychrome)
library(lawstat)
library(tidyverse)

microbiome.raw = read.csv("./0. Raw data input/20221005_wp2/count-table.tsv", sep="\t")

subjectColours = cbind(unique(microbiome.raw$subject), createPalette(63,  c("#ff0000", "#00ff00", "#0000ff")))
colnames(subjectColours) = c("subject", "colour")
rownames(subjectColours) = subjectColours[,1]
subjectColours = subjectColours[,-1]

nicheColours.tongue =  c("#CC79A7", "#BB6B98", "#AB5E88", "#9A5079", "#89426A", "#79355A", "#68274B")
nicheColours.saliva =  c("#00E0A5", "#00C894", "#00B182", "#009971", "#00815F", "#006A4E", "#00523C")
nicheColours.up_prox = c("#89E1EF", "#75CEDC", "#61BBC9", "#4DA8B7", "#3995A4", "#258291", "#116F7E")
nicheColours.low_prox= c("#0E7FE6", "#0C71CD", "#0B63B3", "#09559A", "#074680", "#063867", "#042A4D")
nicheColours.up_ling = c("#F0E442", "#E1D639", "#D2C730", "#C4B927", "#B5AB1E", "#A69C15", "#978E0C")
nicheColours.low_ling= c("#E69F00", "#D49300", "#C28600", "#B07A00", "#9E6D00", "#8C6100", "#7A5400")
RFcolours = c("#898775", "#1100FF", "#FF0000")

rf_data = read.csv("./0. Raw data input/RFdata.csv")
colnames(rf_data) = c("subject", "id", "fotonr", "day", "group", "RFgroup", "MQH", "SPS(tm)", "Area_delta_R30", "Area_delta_Rmax", "Area_delta_R30_x_Rmax", "gingiva_mean_R_over_G", "gingiva_mean_R_over_G_upper_jaw", "gingiva_mean_R_over_G_lower_jaw")
rf_data = rf_data %>% as_tibble()

rf_data[rf_data$subject == "VSTPHZ", 1] = "VSTPH2"
rf_data[rf_data$subject == "D2VZH0", 1] = "DZVZH0"
rf_data[rf_data$subject == "DLODNN", 1] = "DLODDN"
rf_data[rf_data$subject == "O3VQFX", 1] = "O3VQFQ"
rf_data[rf_data$subject == "F80LGT", 1] = "F80LGF"
rf_data[rf_data$subject == "26QQR0", 1] = "26QQrO"

rf_data2 = read.csv("./0. Raw data input/red_fluorescence_data.csv") %>% as_tibble()
rf_data2 = rf_data2[,c(2,4,181:192)]
rf_data = rf_data %>% left_join(rf_data2)

rf = rf_data %>% select(subject, RFgroup) %>% unique()

dataOverviewScatter <- function(df, row_lab, col_lab, titleA, titleB, titleC, titleD){
  
  row_num = 1:nrow(df)
  col_num = 1:ncol(df)
  robustSum = function(x) sum(x, na.rm=TRUE)
  robustStd = function(x) sd(x, na.rm=TRUE)
  
  # Row sum
  plotA = apply(df, 1, robustSum) %>% as_tibble() %>% mutate(row_number=row_num) %>% ggplot(aes(x=row_number,y=value)) + geom_point() + theme(axis.text.x=element_blank()) + ggtitle(titleA) + xlab(row_lab) + ylab("Row sum")
  
  # Col sum
  plotB = apply(df, 2, robustSum) %>% as_tibble() %>% mutate(col_number=col_num) %>% ggplot(aes(x=col_number,y=value)) + geom_point() + theme(axis.text.x=element_blank()) + ggtitle(titleB) + xlab(col_lab) + ylab("Column sum")
  
  # Row std
  plotC = apply(df, 1, robustStd) %>% as_tibble() %>% mutate(row_number=row_num) %>% ggplot(aes(x=row_number,y=value)) + geom_point() + theme(axis.text.x=element_blank()) + ggtitle(titleC) + xlab(row_lab) + ylab("Row std")
  
  # Col std
  plotD = apply(df, 2, robustStd) %>% as_tibble() %>% mutate(col_number=col_num) %>% ggplot(aes(x=col_number,y=value)) + geom_point() + theme(axis.text.x=element_blank()) + ggtitle(titleD) + xlab(col_lab) + ylab("Column std")
  
  ggarrange(plotA, plotB, plotC, plotD, ncol=2, nrow=2, labels=LETTERS[1:4])
}

dataOverviewHist <- function(df, row_lab, col_lab, titleA, titleB, titleC, titleD){
  
  row_num = 1:nrow(df)
  col_num = 1:ncol(df)
  robustSum = function(x) sum(x, na.rm=TRUE)
  robustStd = function(x) sd(x, na.rm=TRUE)
  
  # Row sum
  plotA = apply(df, 1, robustSum) %>% as_tibble() %>% ggplot(aes(x=value)) + geom_histogram(bins=30) + ggtitle(titleA) + xlab("Row sum") + ylab("Counts")
  
  # Col sum
  plotB = apply(df, 2, robustSum) %>% as_tibble() %>% ggplot(aes(x=value)) + geom_histogram(bins=30) + ggtitle(titleB) + xlab("Column sum") + ylab("Counts")
  
  # Row std
  plotC = apply(df, 1, robustStd) %>% as_tibble() %>% ggplot(aes(x=value)) + geom_histogram(bins=30) + ggtitle(titleC) + xlab("Row std") + ylab("Counts")
  
  # Col std
  plotD = apply(df, 2, robustStd) %>% as_tibble() %>% ggplot(aes(x=value)) + geom_histogram(bins=30) + ggtitle(titleD) + xlab("Row std") + ylab("Counts")
  
  ggarrange(plotA, plotB, plotC, plotD, ncol=2, nrow=2, labels=LETTERS[1:4])
}

heteroScedastic = function(df){
  colMeans = colMeans(df, na.rm=TRUE)
  colStds  = apply(df, 2, function(x) sd(x, na.rm=TRUE))
  
  result = cbind(colMeans, colStds)
  result %>% as_tibble() %>% ggplot(aes(x=colMeans, y=colStds)) + geom_point()
}

testNormality <- function(df){
  result = vector(length=ncol(df))
  
  for(i in 1:ncol(df)){
    testResult = shapiro.test(df[,i])
    result[i] = testResult$p.value
  }
  return(result)
}

testSymmetry <- function(df){
  result = vector(length=ncol(df))
  
  for(i in 1:ncol(df)){
    testResult = symmetry.test(df[,i], boot=FALSE)
    result[i] = testResult$p.value
  }
  return(result)
}

importPARAFAC = function(path, featureNames, dayVector, featureColumnOfInterest){
  #  Output:
  #    [[1]]: Loadings in ID mode
  #    [[2]]: Loadings in feature mode
  #    [[3]]: Loadings in time mode
  #    [[4]]: Data as modeled, wide format
  #    [[5]]: Input data, wide format
  #    [[6]]: Data as modeled, long format
  #    [[7]]: Input data, long format
  #    [[8]]: Data as modeled, component 1, wide format
  #    [[9]]: Data as modeled, component 1, long format
  #   [[**]]: and so on until all components have a list item.
  
  id_mode = read.csv(paste0(path,"_individual_mode.csv"), header=FALSE)
  feature_mode = read.csv(paste0(path,"_feature_mode.csv"), header=FALSE)
  time_mode = read.csv(paste0(path,"_time_mode.csv"), header=FALSE)
  numComponents = ncol(time_mode)
  
  if(dim(id_mode)[2] == (numComponents+1)){
    colnames(id_mode) = c(paste0("Component_", 1:numComponents), "subject")
    id_mode = id_mode %>% left_join(rf)
  }
  colnames(id_mode) = c(paste0("Component_", 1:numComponents), "subject", "RFgroup")
  id_mode = as_tibble(id_mode)
  
  colnames(feature_mode) = c(paste0("Component_", 1:numComponents), featureNames)
  feature_mode = as_tibble(feature_mode)
  
  colnames(time_mode) = c(paste0("Component_", 1:numComponents))
  time_mode = as_tibble(time_mode) %>% mutate(days=dayVector)
  
  rawmodel = scan(paste0(path,"_model.csv"), sep=",")
  model_wide = matrix(0, nrow=nrow(id_mode), ncol=nrow(feature_mode)*nrow(time_mode))
  
  for(i in 1:nrow(id_mode)){
    model_wide[i,] = rawmodel[((nrow(feature_mode)*nrow(time_mode)*(i-1))+1):(nrow(feature_mode)*nrow(time_mode)*i)]
  }
  column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
  
  for(i in 2:length(dayVector)){
    column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
  }
  colnames(model_wide) = column_names
  
  model_wide = as_tibble(model_wide)
  
  rawinput = scan(paste0(path,"_input.csv"), sep=",")
  input_wide = matrix(0, nrow=nrow(id_mode), ncol=nrow(feature_mode)*nrow(time_mode))
  
  for(i in 1:nrow(id_mode)){
    input_wide[i,] = rawinput[((nrow(feature_mode)*nrow(time_mode)*(i-1))+1):(nrow(feature_mode)*nrow(time_mode)*i)]
  }
  column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
  
  for(i in 2:length(dayVector)){
    column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
  }
  colnames(input_wide) = column_names
  
  input_wide = as_tibble(input_wide)
  model_long = make_longer(model_wide, id_mode$subject)
  input_long = make_longer(input_wide, id_mode$subject)
  
  result = list(id_mode, feature_mode, time_mode, model_wide, input_wide, model_long, input_long)
  listIterator = 8
  
  for(i in 1:numComponents){
    modeled_component_raw = scan(paste0(path,"_component_", i, ".csv"), sep=",")
    modeled_component_wide = matrix(0, nrow=nrow(id_mode), ncol=nrow(feature_mode)*nrow(time_mode))
    
    for(i in 1:nrow(id_mode)){
      modeled_component_wide[i,] = modeled_component_raw[((nrow(feature_mode)*nrow(time_mode)*(i-1))+1):(nrow(feature_mode)*nrow(time_mode)*i)]
    }
    column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
    
    for(i in 2:length(dayVector)){
      column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
    }
    colnames(modeled_component_wide) = column_names
    
    modeled_component_wide = as_tibble(modeled_component_wide)
    modeled_component_long = make_longer(modeled_component_wide, id_mode$subject)
    
    result[[listIterator]] = modeled_component_wide
    listIterator = listIterator + 1
    result[[listIterator]] = modeled_component_long
    listIterator = listIterator + 1
  }
  
  return(result)
}

importNPLS = function(path, featureNames, dayVector, featureColumnOfInterest){
  id_mode = read.csv(paste0(path,"_individual_mode.csv"), header=FALSE)
  feature_mode = read.csv(paste0(path,"_feature_mode.csv"), header=FALSE)
  time_mode = read.csv(paste0(path,"_time_mode.csv"), header=FALSE)
  numComponents = ncol(time_mode)
  
  colnames(id_mode) = c(paste0("Component_", 1:numComponents), "subject", "RFgroup")
  id_mode = as_tibble(id_mode)
  
  colnames(feature_mode) = c(paste0("Component_", 1:numComponents), featureNames)
  feature_mode = as_tibble(feature_mode)
  
  colnames(time_mode) = c(paste0("Component_", 1:numComponents))
  time_mode = as_tibble(time_mode) %>% mutate(days=dayVector)
  
  model = read.csv(paste0(path,"_model.csv"), header=FALSE)
  dim(model)
  
  ypred = read.csv(paste0(path, "_ypred.csv"), header=FALSE)
  colnames(ypred) = c(paste0("Component_", 1:numComponents), "subject", "y")
  ypred = as_tibble(ypred)
  
  
  column_names = paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t1")
  for(i in 2:length(dayVector)){
    column_names = c(column_names, paste0(feature_mode[featureColumnOfInterest] %>% pull, "_t", i))
  }
  column_names
  colnames(model) = column_names
  
  model = as_tibble(model)
  
  return(list(id_mode, feature_mode, time_mode, model, ypred))
}

make_longer = function(data, subjects){
  df = data %>% mutate(subject=subjects) %>% pivot_longer(-subject, names_to=c("asv", "visit"), names_sep="_t") %>% as.data.frame()
  df$visit = as.integer(df$visit)
  df$value = as.double(df$value)
  df = df %>% as_tibble() %>% pivot_wider(id_cols=c("subject","visit"), names_from=asv, values_from=value) %>% select(-subject,-visit)
  return(df)
}