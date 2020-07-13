options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)

# Inputs
args = c("A", "-60", "60")


# Function to convert output from cut, if created with SummariseAB.R
getClass <- function(x, y){
  x = as.numeric(x)
  y = as.numeric(y)
  classNames = list()
  classNames[[paste("(-Inf,", x - 1, "]", sep = "")]] = paste("<", x, sep = '')
  for (c in seq( x, y )){
    classNames[[paste("(", c - 1, ",", c, "]", sep = "") ]] = c
  }
  classNames[[paste("(", y, ", Inf]", sep = "")]] = paste(">", y, sep = '')
  return(classNames)
}

# 
classVec <- function(x, y){
  classes = c(paste("<", args[2], sep = ''))
  for (cl in seq(as.numeric(args[2]), as.numeric(args[3]) )){
    classes =  c(classes, c(cl))
  }
  classes = c(classes, c(paste(">", args[3], sep = '')))
  return(classes)
}


# Get algorithms in the folder
getAlgorithms <- function(){
  algorithms = list.dirs(".", recursive = F)
  if (!dir.exists("PLOTS")){
    outdir = dir.create("PLOTS")
  } else {
    algorithms = algorithms[!grepl("PLOTS", algorithms)]
  }
  return(algorithms)
}

# Create value dataframe
makeValues <- function(samples, sample_dirs, algorithm){
  values = data.frame(GROUP =  classes)
  values$GROUP <- factor(values$GROUP, levels = values$GROUP)
  for (n in seq(1, nrow(samples))){
    sample = samples[n, 2]
    dataset = read.csv(list.files(sample_dirs[n], full.names = T), h=T) 
    dataset$GROUP <- factor(dataset$GROUP, levels = dataset$GROUP)
    dataset = dataset %>% select(GROUP, avgDP, N)
    colnames(dataset)[2] = paste(sample, "DP", sep = "_")
    colnames(dataset)[3] = paste(sample, "N", sep = "_")
    if (!exists("values")){
      values = dataset
    } else {
      values = merge(values, dataset, all.x = T, all.y = T, by = "GROUP")
    }
  }
  values["MEAN_DP"] = rowMeans(values[,seq(2, ncol(values), 2)], na.rm=TRUE) 
  values["MEAN_Nv"] = rowMeans(values[,seq(3, ncol(values), 2)], na.rm=TRUE) 
  values = values %>% select(GROUP, MEAN_DP, MEAN_Nv)
  colnames(values)[c(2,3)] = c( 
    paste(gsub("./","", algorithm), "DP", sep = "_"), 
    paste(gsub("./","", algorithm), "Nv", sep = "_") 
  )
  return(values)
}

# Create value dataframe
makeCounts <- function(samples, sample_dirs, algorithm, breed){
  values = data.frame(GROUP =  classes)
  values$GROUP <- factor(values$GROUP, levels = values$GROUP)
  for (n in seq(1, nrow(samples))){
    sample = samples[n, 2]
    dataset = read.csv(list.files(sample_dirs[n], full.names = T), h=T) 
    dataset$GROUP <- factor(dataset$GROUP, levels = dataset$GROUP)
    dataset = dataset %>% select(GROUP, N)
    colnames(dataset)[2] = sample
    if (!exists("values")){
      values = dataset
    } else {
      values = merge(values, dataset, all.x = T, all.y = T, by = "GROUP")
    }
  }
  return(values)
}


# Start the script
classNames = getClass( args[2], args[3] )
algorithms = getAlgorithms()


for (G in c("ALL", "NOVEL", "KNOWN") ){
  if ( !dir.exists( paste("PLOTS/", G, sep = "") ) ){
    outdir = dir.create( paste("PLOTS/", G, sep = "" ))
  } 
  for (br in c("NDama", "Angus", "Sahiwal")){
    if (exists("final")){rm(final)}
    if (exists("values")){rm(values)}
    
    classes = classVec(args[2], args[3])
    tmp_algorithms = algorithms[grepl(G,algorithms)]
    for (algorithm in tmp_algorithms){
      sample_dirs = list.dirs(algorithm, recursive = F)
      sample_dirs = sample_dirs[grepl(br,sample_dirs)]
      samples = data.frame(dir = sample_dirs) %>% 
        separate(dir, sep = "/", into = c("CWD","ALGORITHM", "SAMPLE")) %>%
        select(ALGORITHM, SAMPLE)
      
      values = makeValues(samples, sample_dirs, algorithm)
      if (!exists("final")){
        final = values
      } else {
        final = merge(final, values, by = "GROUP", all.x = T)
      }
      #rm(values)
    }
    
    to_plot_tmp1 = melt(final[,c("GROUP", colnames(final)[grepl("_DP", colnames(final))] )], id.vars = "GROUP") %>% 
      select("Length" = GROUP, "DP"=value, "Algorithm"=variable)
    to_plot_tmp1$Algorithm = gsub("_DP","", to_plot_tmp1$Algorithm)
    to_plot_tmp2 = melt(final[,c("GROUP", colnames(final)[grepl("_Nv", colnames(final))] )], id.vars = "GROUP") %>% 
      select("Length" = GROUP, "NumberVariants"=value, "Algorithm"=variable)
    to_plot_tmp2$Algorithm = gsub("_Nv","", to_plot_tmp2$Algorithm)
    to_plot = merge(to_plot_tmp1, to_plot_tmp2, by = c("Length", "Algorithm"))
    to_plot$NumberVariants = log10(to_plot$NumberVariants)
    p = to_plot %>% 
      ggplot(aes(x=Length, y=DP, color = Algorithm)) +
      geom_point() + geom_line(aes(group = Algorithm)) + 
      scale_y_continuous() +
      theme(axis.text.x = element_text(angle=90, hjust=1)) +
      labs(x="Length", y="DP (avg)", color = "Software", title = paste("DP of variant by variant size and number")) 
    ggsave(paste("PLOTS/", G, "/Number_Variants_DP_", br,"_",G,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
  }
  
  if (exists("final")){rm(final)}
  if (exists("values")){rm(values)}
  classes = classVec(args[2], args[3])
  final = data.frame(GROUP =  classes)
  final$GROUP <- factor(final$GROUP, levels = final$GROUP)
  for (algorithm in tmp_algorithms){
    sample_dirs = list.dirs(algorithm, recursive = F)
    samples = data.frame(dir = sample_dirs) %>% 
      separate(dir, sep = "/", into = c("CWD","ALGORITHM", "SAMPLE")) %>%
      select(ALGORITHM, SAMPLE)
    
    values = makeValues(samples, sample_dirs, algorithm)
    
    if (!exists("final")){
      final = values
    } else {
      final = merge(final, values, by = "GROUP", all.x = T)
    }
    rm(values)
  }
  
  to_plot_tmp1 = melt(final[,c("GROUP", colnames(final)[grepl("_DP", colnames(final))] )], id.vars = "GROUP") %>% 
    select("Length" = GROUP, "DP"=value, "Algorithm"=variable)
  to_plot_tmp1$Algorithm = gsub("_DP","", to_plot_tmp1$Algorithm)
  to_plot_tmp2 = melt(final[,c("GROUP", colnames(final)[grepl("_Nv", colnames(final))] )], id.vars = "GROUP") %>% 
    select("Length" = GROUP, "NumberVariants"=value, "Algorithm"=variable)
  to_plot_tmp2$Algorithm = gsub("_Nv","", to_plot_tmp2$Algorithm)
  to_plot = merge(to_plot_tmp1, to_plot_tmp2, by = c("Length", "Algorithm"))
  to_plot$NumberVariants = log10(to_plot$NumberVariants)
  p = to_plot %>% 
    ggplot(aes(x=Length, y=DP, color = Algorithm)) +
    geom_point() + geom_line(aes(group = Algorithm)) + 
    scale_y_continuous() +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    labs(x="Length", y="DP (avg)", color = "Software", title = paste("DP of variant by variant size and number")) 
  ggsave(paste("PLOTS/", G, "/Number_Variants_DP_WGS_",G,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
}


# Make general count plots
for (G in c("ALL", "NOVEL", "KNOWN") ){
  if (exists("tmp")){rm(tmp)}
  if (exists("final")){rm(final)}
  if ( !dir.exists( paste("PLOTS/", G, sep = "") ) ){
    outdir = dir.create( paste("PLOTS/", G, sep = "" ))
  } 
  tmp_algorithms = algorithms[grepl(G,algorithms)]
  for (br in c("NDama", "Angus", "Sahiwal")){
    classes = classVec(args[2], args[3])
    if (exists("tmp")){rm(tmp)}
    for (algorithm in tmp_algorithms){
      sample_dirs = list.dirs(algorithm, recursive = F)
      sample_dirs = sample_dirs[grepl(br,sample_dirs)]
      samples = data.frame(dir = sample_dirs) %>% 
        separate(dir, sep = "/", into = c("CWD","ALGORITHM", "SAMPLE")) %>%
        select(ALGORITHM, SAMPLE)
      
      values = makeCounts(samples, sample_dirs, algorithm, br)
      nvals = data.frame(colSums(values[,2:ncol(values)], na.rm = T))
      nvals$SAMPLE = rownames(nvals)
      nvals = nvals[,c(2,1)]
      rownames(nvals) = NULL
      colnames(nvals) = c("SAMPLE", paste(gsub("./", "", algorithm)))
      if (!exists("tmp")){
        tmp = nvals
      } else {
        tmp = merge(tmp, nvals, all.x = T, all.y = T, by= "SAMPLE")
      }
    }
    tmp = tmp[,c("SAMPLE", gsub("./", "", tmp_algorithms)) ]
    if (!exists("final")){
      final = tmp
    } else {
      final = rbind(final, tmp)
    }
  }
  to_plot = melt(final, id.vars = "SAMPLE") %>% select("Sample" = SAMPLE, "VariantsNumber"=value, "Algorithm"=variable)
  p = to_plot %>% 
    ggplot(aes(x=Sample, y=VariantsNumber, color = Algorithm)) +
    geom_point() + geom_line(aes(group = Algorithm)) + geom_hline(aes(yintercept = 0.5), linetype="dashed", colour = "red") +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    labs(x="Sample", y="Variants", color = "Software", title = paste("Number of variant by variant size")) 
  ggsave(paste("PLOTS/", G, "/Summary_NumberVariants_WGS_",G,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
}

