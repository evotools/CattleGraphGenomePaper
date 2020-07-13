options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)


# Get algorithms in the folder
getAlgorithms <- function(f = "."){
  algorithms = list.dirs(f, recursive = F, full.names = T)
  if (dir.exists(paste(f, "PLOTS", sep = "/"))){
    algorithms = algorithms[!grepl("PLOTS", algorithms)]
  } else {
    dir.create( paste(f, "PLOTS", sep = '/') )
  }
  return(algorithms)
}

# Create value dataframe
makeCounts <- function(samples, sample_dirs, algorithm, breed){
  for (n in seq(1, nrow(samples))){
    sample = samples[n, 2]
    dataset = read.csv(list.files(sample_dirs[n], full.names = T), h=T)
    dataset$N_AltGT = rowSums( dataset[, c("N_RA", "N_AA_HOM", "N_AA_HET") ] )
    dataset$SIZE <- factor(dataset$SIZE, levels = dataset$SIZE)
    dataset = dataset %>% select(SIZE, N_AltGT)
    colnames(dataset)[2:ncol(dataset)] = sample
    if (!exists("values")){
      values = dataset
    } else {
      values = merge(values, dataset, all.x = T, all.y = T, by = "SIZE")
    }
  }
  return(values)
}


# Make general count plots
algorithms = getAlgorithms()
for (G in c("ALL", "NOVEL", "KNOWN") ){

  if (exists("final")){rm(final)}
  if (exists("tmp")){rm(tmp)}
  if ( !dir.exists( paste("./PLOTS/", sep = "") ) ){
    outdir = dir.create( paste("./PLOTS/", sep = "" ))
  } 
  if ( !dir.exists( paste("./PLOTS/", G, sep = "") ) ){
    outdir = dir.create( paste("./PLOTS/", G, sep = "" ))
  } 
  
  tmp_algorithms = algorithms[grepl(G,algorithms)]
  for (br in c("NDama", "Angus", "Sahiwal")){
    if (exists("tmp")){rm(tmp)}
    for (algorithm in tmp_algorithms){
      sample_dirs = list.dirs(algorithm, recursive = F)
      sample_dirs = sample_dirs[grepl(br,sample_dirs)]
      sample_dirs = gsub("//","/", sample_dirs)
      samples = data.frame(dir = sample_dirs) %>% 
        separate(dir, sep = "/", into = c("CWD","ALGORITHM", "SAMPLE")) %>%
        select(ALGORITHM, SAMPLE)
      
      values = makeCounts(samples, sample_dirs, algorithm, br) 
      values = melt(values, id.vars = "SIZE") %>% 
        select("VariantSize" = SIZE, "N_AltGT"=value, "Sample"=variable)
      values$Algorithm = gsub("./","",algorithm)
      values$Breed = br
      nvals = values
      #for (coln in c(2:ncol(values)) ){
      #  values[,coln] = values[,coln]/ sum(values[,coln])
      #}
      
      #values[gsub("./","",algorithm)] = rowMeans(values[,seq(2, ncol(values), 2)], na.rm=TRUE) 
      #nvals = values[,c(1,ncol(values))]
      if (!exists("tmp")){
        tmp = nvals
      } else {
        tmp = rbind(tmp, nvals)
      }
      rm(values)
    }

    to_plot=tmp
    #to_plot = melt(tmp, id.vars = "SIZE") %>% select("VariantSize" = SIZE, "N_AltGT"=value, "Algorithm"=variable)
    # p = to_plot %>% 
    #   ggplot(aes(x=VariantSize, y=N_AltGT, color = Algorithm)) +
    #   geom_point() + geom_line(aes(group = Algorithm)) + geom_hline(aes(yintercept = 2.0), linetype="dashed", colour = "red") +
    #   theme(axis.text.x = element_text(angle=90, hjust=1)) + 
    #   scale_y_log10 +
    #   labs(x="Variant size", y="Proportion of variants", color = "Software", title = paste("Variants rate by algorithm", br)) 
    to_plot2 = to_plot
    to_plot2$N_AltGT = log10(to_plot2$N_AltGT)
    to_plot2$N_AltGT[!is.finite(to_plot2$N_AltGT) ] = 0
    to_plot2$Algorithm = factor(to_plot2$Algorithm, levels = c("GATK_BQSR_ALL", "FB_BQSR_ALL", "VG_LINEAR_ALL", "VG_CACTUS_ALL", "VG_DIVERSITY_ALL"))
    p = to_plot %>% 
      ggplot(aes(x=VariantSize, y=N_AltGT, color = Algorithm, fill = Algorithm)) +
      geom_dotplot(binaxis='y', stackdir='center') + 
      theme(axis.text.x = element_text(angle=90, hjust=1)) + 
      labs(x="Variant size", y="Variants #", title = paste("Variants by size by algorithm", br)) #+
      facet_wrap(to_plot$Algorithm, ncol = 5) 
      #scale_y_log10() + 
    
    ggsave(paste("./PLOTS/", G, "/Summary_Sizes_WGS_",G,"_",br,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
    if (exists("final")){
      final = rbind(final, tmp)
    } else {
      final = tmp
    }
    rm(tmp)
    
  }
  write_excel_csv(final, "results.size.csv")
  p = final %>% 
    ggplot(aes(x=Breed, y=N_AltGT, color = Algorithm, fill = Algorithm)) +
    geom_dotplot(binaxis='y', stackdir='center') + 
    theme(axis.text.x = element_text(angle=90, hjust=1)) + 
    labs(x="Variant size", y="Variants #", title = paste("Variants by size by algorithm", br)) +
    facet_wrap(to_plot$VariantSize, ncol = 5) 
  ggsave(paste("./PLOTS/", G, "/Summary_Sizes_WGS_",G,"_ALL.pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
  
}


# Make general count plots
algorithms = getAlgorithms()
for (G in c("ALL", "NOVEL", "KNOWN") ){
  if (exists("final")){rm(final)}
  if (exists("tmp")){rm(tmp)}
  if ( !dir.exists( paste("./PLOTS/", sep = "") ) ){
    outdir = dir.create( paste("./PLOTS/", sep = "" ))
  } 
  if ( !dir.exists( paste("./PLOTS/", G, sep = "") ) ){
    outdir = dir.create( paste("./PLOTS/", G, sep = "" ))
  } 
  tmp_algorithms = algorithms[grepl(G,algorithms)]
  if (exists("tmp")){rm(tmp)}
  for (algorithm in tmp_algorithms){
    sample_dirs = list.dirs(algorithm, recursive = F)
    sample_dirs = sample_dirs[grepl(br,sample_dirs)]
    sample_dirs = gsub("//","/", sample_dirs)
    samples = data.frame(dir = sample_dirs) %>% 
      separate(dir, sep = "/", into = c("CWD","ALGORITHM", "SAMPLE")) %>%
      select(ALGORITHM, SAMPLE)
    
    values = makeCounts(samples, sample_dirs, algorithm, br) 
    values[paste("NAltGT", gsub("./","",algorithm), sep = "_")] = rowMeans(values[,seq(2, ncol(values), 2)], na.rm=TRUE) 
    nvals = values[,c(1,ncol(values))]
    if (!exists("tmp")){
      tmp = nvals
    } else {
      tmp = merge(tmp, nvals, all.x = T, all.y = T, by= "SIZE")
    }
    rm(values)
  }
  to_plot = melt(tmp, id.vars = "SIZE") %>% select("VariantSize" = SIZE, "N_AltGT"=value, "Algorithm"=variable)
  to_plot$Algorithm = gsub("_ALL", "", gsub("NAltGT_", "", to_plot$Algorithm))
  p = to_plot %>% 
    ggplot(aes(x=VariantSize, y=N_AltGT, color = Algorithm)) +
    geom_point() + geom_line(aes(group = Algorithm)) + geom_hline(aes(yintercept = 2.0), linetype="dashed", colour = "red") +
    theme(axis.text.x = element_text(angle=90, hjust=1)) + scale_y_log10() +
    labs(x="Variant size (bp)", y="# Alternative genotypes (log10)", color = "Software", title = "Number of variants by variant size") 
  ggsave(paste("./PLOTS/", G, "/Summary_Sizes_WGS_",G, ".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
  rm(tmp)
}






