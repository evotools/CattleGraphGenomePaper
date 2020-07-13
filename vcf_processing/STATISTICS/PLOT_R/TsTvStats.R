options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)


# Get algorithms in the folder
getAlgorithms <- function(f = "."){
  algorithms = list.dirs(f, recursive = F, full.names = T)
  if (dir.exists("PLOTS")){
    algorithms = algorithms[!grepl("PLOTS", algorithms)]
  } else {
    dir.create( paste(f, "PLOTS") )
  }
  return(algorithms)
}

# Create value dataframe
makeCounts <- function(samples, sample_dirs, algorithm, breed){
  for (n in seq(1, nrow(samples))){
    sample = samples[n, 2]
    dataset = read.table(list.files(sample_dirs[n], full.names = T), h=T) 
    dataset$MODEL <- factor(dataset$MODEL, levels = dataset$MODEL)
    colnames(dataset)[2] = sample
    if (!exists("values")){
      values = dataset
    } else {
      values = merge(values, dataset, all.x = T, all.y = T, by = "MODEL")
    }
  }
  return(values)
}


# Make general count plots
for (FILT in list.dirs(recursive = F)){
  algorithms = getAlgorithms(FILT)
  for (G in c("ALL", "NOVEL", "KNOWN") ){
    if (exists("final")){rm(final)}
    if (exists("tmp")){rm(tmp)}
    if ( !dir.exists( paste(FILT, "/PLOTS/", sep = "") ) ){
      outdir = dir.create( paste(FILT, "/PLOTS/", sep = "" ))
    } 
    if ( !dir.exists( paste(FILT,"/PLOTS/", G, sep = "") ) ){
      outdir = dir.create( paste(FILT,"/PLOTS/", G, sep = "" ))
    } 
    tmp_algorithms = algorithms[grepl(G,algorithms)]
    for (br in c("NDama", "Angus", "Sahiwal")){
      if (exists("tmp")){rm(tmp)}
      for (algorithm in tmp_algorithms){
        sample_dirs = list.dirs(algorithm, recursive = F)
        sample_dirs = sample_dirs[grepl(br,sample_dirs)]
        samples = data.frame(dir = sample_dirs) %>% 
          separate(dir, sep = "/", into = c("CWD","FILT","ALGORITHM", "SAMPLE")) %>%
          select(ALGORITHM, SAMPLE)
        
        values = makeCounts(samples, sample_dirs, algorithm, br) %>% filter(MODEL == "Ts" | MODEL == "Tv")
        tstv = values[1,seq(2, ncol(values))] / values[2, seq(2, ncol(values))]
        nvals = data.frame(Samples = colnames(values)[seq(2, ncol(values))], TsTv = as.numeric(as.vector(tstv[1,])) )
        colnames(nvals) = c("SAMPLE", paste(gsub( paste(FILT, "/", sep = "") , "", algorithm, fixed = T)))
        if (!exists("tmp")){
          tmp = nvals
        } else {
          tmp = merge(tmp, nvals, all.x = T, all.y = T, by= "SAMPLE")
        }
        rm(values)
      }
      tmp = tmp[,c("SAMPLE", gsub(paste(FILT, "/", sep = ""), "", tmp_algorithms)) ]
      if (!exists("final")){
        final = tmp
      } else {
        final = rbind(final, tmp)
      }
      to_plot = melt(tmp, id.vars = "SAMPLE") %>% select("Sample" = SAMPLE, "TsTv"=value, "Algorithm"=variable)
      p = to_plot %>% 
        ggplot(aes(x=Sample, y=TsTv, color = Algorithm)) +
        geom_point() + geom_line(aes(group = Algorithm)) + geom_hline(aes(yintercept = 2.0), linetype="dashed", colour = "red") +
        theme(axis.text.x = element_text(angle=90, hjust=1)) +
        labs(x="Sample", y="Ts/Tv", color = "Software", title = paste("Transition/Trasversion rate", br)) 
      ggsave(paste(FILT,"/PLOTS/", G, "/Summary_TsTv_WGS_",G,"_",br,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
      rm(tmp)
    }
    to_plot = melt(final, id.vars = "SAMPLE") %>% select("Sample" = SAMPLE, "TsTv"=value, "Algorithm"=variable)
    p = to_plot %>% 
      ggplot(aes(x=Sample, y=TsTv, color = Algorithm)) +
      geom_point() + geom_line(aes(group = Algorithm)) + geom_hline(aes(yintercept = 2.0), linetype="dashed", colour = "red") +
      theme(axis.text.x = element_text(angle=90, hjust=1)) +
      labs(x="Sample", y="Ts/Tv", color = "Software", title = paste("Transition/Trasversion rate")) 
    ggsave(paste(FILT,"/PLOTS/", G, "/Summary_TsTv_WGS_",G,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
  }
}

