library(tidyverse)
library(ggpubr)
args = commandArgs(trailingOnly=TRUE)
options(scipen = 999)

positive = as.numeric(args[1])
numPerm<-as.numeric(args[2])
print(paste("Positive:", positive))
print(paste("Permutations:", numPerm))

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

lengths = as.numeric(abs(read.table("lengths.txt", h=F)[,1]))

found = c()
for (test in seq(1,numPerm)){
  regions = generateRandomPos(myfaifile$V1, myfaifile$V2, lengths)
  regions = as.data.frame(cbind(regions, paste("REGION", c(1: nrow(regions)), sep = "")))
  write.table(regions, paste("myregion",test,"bed", sep = '.'), sep = "\t", row.names = F, col.names = F, quote = F)
  system(paste("bedtools intersect -a myregion",test,"bed -b OMregions.bed | sort | uniq > intersected",test,"bed", sep = "."), wait = T) 
  tryCatch({
    nfound = nrow(read.table(paste("intersected",test,"bed", sep = '.'), h=F))
    }, error = function(cond){nfound = 0})
  if (!exists("nfound")){nfound = 0}
  system(paste("rm myregion", test, "bed intersected", test, "bed", sep = '.'))
  found = c(found, nfound)
  if (test%%500 == 0){
    print(test)
    flush.console()
  }
}
print(length(found))
p = ggplot(data.frame(found=found), aes(x=found)) + 
  geom_histogram(binwidth = 1, aes(y=..density..), colour="black", fill="white", )+
  geom_density(adjust = 5, alpha=.2, fill="#FF6666") 
ggsave(filename = "dist.pdf", device = "pdf", plot = p, width = 12, height = 8)

Perm_mean = mean(found)
Perm_sd = sd(found)
print(Perm_mean)
print(Perm_sd)
z=(positive-Perm_mean)/Perm_sd 
z_P<-2*pnorm(-abs(z))
print(z)
print(z_P)
results = data.frame(PARAM = c("Positive", "Total", "Ntests", "Mean", "StD", "Z", "P-value"), 
                     VALUES = c(positive, length(lengths), numPerm, Perm_mean, Perm_sd, z, z_P) )

write.table(found,"matches.txt", sep = "\t", quote =F, col.names= F, row.names = F)
final_p = write.table(results, "pvalue.txt")


