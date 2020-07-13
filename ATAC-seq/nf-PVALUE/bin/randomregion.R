#!/usr/bin/env Rscript
library(tidyverse)
library(GenomicRanges)
library(ggpubr)
args = commandArgs(trailingOnly=TRUE)


positive = 6129
numPerm<-10000


getWinds<-function(thisChrom, thisWidth, chr, chr.sizes)
{
  sample(chr.sizes[chr==thisChrom]-thisWidth,1)
}

generateRandomPos <- function(chr, chr.sizes, widths)
{
  random_chr <- sample(x=chr,size=length(widths),prob=chr.sizes,replace=T)
  random_pos<-mapply(getWinds, random_chr, widths, MoreArgs = list(chr, chr.sizes))
  res<-cbind(random_chr, random_pos, random_pos+widths)
  return(res)
}

# Get the fai file
myfaifile = read.table("ARS_UCD1.2.fai", h=F) %>%
  filter(V1 %in% c(1:29))

lengths = as.numeric(abs(read.table("lengths.txt")[,1]))

found = c()
for (test in seq(1,numPerm)){
  regions = generateRandomPos(myfaifile$V1, myfaifile$V2, lengths)
  regions = as.data.frame(cbind(regions, paste("REGION", c(1: nrow(regions)), sep = "")))
  write.table(regions, paste("myregion",test,"bed", sep = '.'), sep = "\t", row.names = F, col.names = F, quote = F)
  system(paste("bedtools intersect -a myregion",test,"bed -b OMregions.bed | sort | uniq > intersected",test,"bed", sep = "."), wait = T) 
  tryCatch({
    nfound = nrow(read.table(paste("intersected",test,"bed", sep = '.'), h=F))
    system(paste("rm myregion", test, "bed intersected", test, "bed", sep = '.'))
    }, error = function(cond){nfound = 0})
  found = c(found, nfound)
  print(test)
  flush.console() 
}
ggplot(data.frame(found=found), aes(x=found)) + 
  geom_histogram(binwidth = 1, aes(y=..density..), colour="black", fill="white", )+
  geom_density(adjust = 5, alpha=.2, fill="#FF6666") 

Perm_mean = mean(found)
Perm_sd = sd(found)
z=(positive-Perm_mean)/Perm_sd  
z_P<-2*pnorm(-abs(z))

tot = nrow(regions)
negative = tot-positive
pvals = c()
for (ctrl_pos in found){
  ctrl_neg = tot - ctrl_pos
  csq = chisq.test(c(positive, negative), p=c(ctrl_pos/tot, ctrl_neg/tot))
  pvals = c(pvals, csq$p.value)
}
write.table(found,"matches.txt", sep = "\t", quote =F, col.names= F, row.names = F)
write.table(pvals,"pvals.txt", sep = "\t", quote =F, col.names= F, row.names = F)
final_p = write.table(z_P, "pvalue.txt")


