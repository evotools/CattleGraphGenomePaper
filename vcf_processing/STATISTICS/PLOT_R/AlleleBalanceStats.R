options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)
args = commandArgs(T)

# Inputs
if (length(args) == 0){
  args = c("A", "-60", "60")
}


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
    dataset = dataset %>% select(GROUP, AVG)
    colnames(dataset)[2] = sample
    if (!exists("values")){
      values = dataset
    } else {
      values = merge(values, dataset, all.x = T, all.y = T, by = "GROUP")
    }
  }
  values["MEAN"] = rowMeans(values[,seq(2, ncol(values))], na.rm=TRUE) 
  values = values %>% select(GROUP, MEAN)
  colnames(values)[2] = gsub("./","", algorithm)
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
      rm(values)
    }
    
    to_plot = melt(final, id.vars = "GROUP") %>% select("Length" = GROUP, "AllelicBalance"=value, "Algorithm"=variable)
    to_plot$Algorithm = gsub("_BQSR", "", gsub(paste("_",G, sep = ""), "", to_plot$Algorithm))
    p = to_plot %>% 
      ggplot(aes(x=Length, y=AllelicBalance, color = Algorithm)) +
      geom_point() + geom_line(aes(group = Algorithm)) + geom_hline(aes(yintercept = 0.5), linetype="dashed", colour = "red")  +
      scale_y_continuous(breaks=seq(0, 1, .04), limits = c(0.48, 0.88)) +
      theme(axis.text.x = element_text(angle=90, hjust=1)) +
      labs(x="Length", y="Allelic balance", color = "Software", title = paste("Allelic balance by variant size", br)) 
    ggsave(paste("PLOTS/", G, "/Allelic_Balance_", br,"_",G,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
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
    
    # Compute values
    values = makeValues(samples, sample_dirs, algorithm)
    
    if (!exists("final")){
      final = values
    } else {
      final = merge(final, values, by = "GROUP", all.x = T)
    }
    rm(values)
  }
  
  to_plot = melt(final, id.vars = "GROUP") %>% select("Length" = GROUP, "AllelicBalance"=value, "Algorithm"=variable)
  p = to_plot %>% 
    ggplot(aes(x=Length, y=AllelicBalance, color = Algorithm)) +
    geom_point() + geom_line(aes(group = Algorithm)) + geom_hline(aes(yintercept = 0.5), linetype="dashed", colour = "red") +
    scale_y_continuous(breaks=seq(0, 1, .04), limits = c(0.48, 0.88)) +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    labs(x="Length", y="Allelic balance", color = "Software", title = paste("Allelic balance by variant size")) 
  ggsave(paste("PLOTS/", G, "/Allelic_Balance_WGS_",G,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
}


for (br in c("NDama", "Angus", "Sahiwal")){
  if ( !dir.exists( paste("PLOTS/", br, sep = "") ) ){
    outdir = dir.create( paste("PLOTS/", br, sep = "" ))
  } 
  if (exists("final")){rm(final)}
  if (exists("values")){rm(values)}
  # Read algorithm of interest
  classes = classVec(args[2], args[3])
  for (algorithm in c("./VG_DIVERSITY_ALL", 
                      "./VG_DIVERSITY_KNOWN", 
                      "./VG_DIVERSITY_NOVEL", 
                      "./FB_BQSR_ALL", 
                      "./GATK_BQSR_ALL")){
    sample_dirs = list.dirs(algorithm, recursive = F)
    sample_dirs = sample_dirs[grepl(br,sample_dirs)]
    samples = data.frame(dir = sample_dirs) %>% 
      separate(dir, sep = "/", into = c("CWD","ALGORITHM", "SAMPLE")) %>%
      select(ALGORITHM, SAMPLE)
    samples$ALGORITHM = gsub("_BQSR", "", gsub("_NOVEL", "", gsub("_KNOWN", "",gsub("_ALL", "", samples$ALGORITHM))))
    
    values = makeValues(samples, sample_dirs, algorithm)
    
    if (!exists("final")){
      final = values
    } else {
      final = merge(final, values, by = "GROUP", all.x = T)
    }
    rm(values)
  }
  
  to_plot = melt(final, id.vars = "GROUP") %>% select("Length" = GROUP, "AllelicBalance"=value, "Algorithm"=variable)
  to_plot$Algorithm = gsub("_BQSR_ALL", "", to_plot$Algorithm)
  
  p = to_plot %>% 
    ggplot(aes(x=Length, y=AllelicBalance, color = Algorithm)) +
    geom_point() + geom_line(aes(group = Algorithm)) + geom_hline(aes(yintercept = 0.5), linetype="dashed", colour = "red")  +
    scale_y_continuous(breaks=seq(0, 1, .04), limits = c(0.48, 0.88)) +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    labs(x="Length", y="Allelic balance", color = "Software", title = paste("Allelic balance by variant size", br)) 
  ggsave(paste("PLOTS/", br, "/Allelic_Balance_", br,"_",G,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
}
  
if (exists("final")){rm(final)}
if (exists("values")){rm(values)}
if ( !dir.exists( paste("PLOTS/", "WGS", sep = "") ) ){
  outdir = dir.create( paste("PLOTS/", "WGS", sep = "" ))
} 
classes = classVec(args[2], args[3])
final = data.frame(GROUP =  classes)
final$GROUP <- factor(final$GROUP, levels = final$GROUP)
for (algorithm in c("./VG_DIVERSITY_ALL", 
                    "./VG_DIVERSITY_KNOWN", 
                    "./VG_DIVERSITY_NOVEL", 
                    "./FB_BQSR_ALL", 
                    "./GATK_BQSR_ALL")){
  sample_dirs = list.dirs(algorithm, recursive = F)
  samples = data.frame(dir = sample_dirs) %>% 
    separate(dir, sep = "/", into = c("CWD","ALGORITHM", "SAMPLE")) %>%
    select(ALGORITHM, SAMPLE)
  samples$ALGORITHM = gsub("_BQSR", "", gsub("_NOVEL", "", gsub("_KNOWN", "",gsub("_ALL", "", samples$ALGORITHM))))
  
  # Compute values
  values = makeValues(samples, sample_dirs, algorithm)
  
  if (!exists("final")){
    final = values
  } else {
    final = merge(final, values, by = "GROUP", all.x = T)
  }
  rm(values)
}

to_plot = melt(final, id.vars = "GROUP") %>% select("Length" = GROUP, "AllelicBalance"=value, "Algorithm"=variable)
to_plot$Algorithm = gsub("_BQSR_ALL", "", to_plot$Algorithm)

p = to_plot %>% 
  ggplot(aes(x=Length, y=AllelicBalance, color = Algorithm)) +
  geom_point() + geom_line(aes(group = Algorithm)) + geom_hline(aes(yintercept = 0.5), linetype="dashed", colour = "red") +
  scale_y_continuous(breaks=seq(0, 1, .04), limits = c(0.48, 0.88)) +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  labs(x="Length", y="Allelic balance", color = "Software", title = paste("Allelic balance by variant size")) 
ggsave(paste("PLOTS/WGS/Allelic_Balance_WGS_",G,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)


