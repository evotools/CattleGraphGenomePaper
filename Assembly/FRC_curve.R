# Make coverage plots
args = commandArgs(T)
library(tidyverse)
library(ggpubr)
library(ggsci)
frc = read.table(args[1])
colnames(frc) = c("Feature", "Coverage")
pdf(paste(args[3], 'pdf',sep = '.'), width = 12, height = 8)
ggline(data = frc, x = "Feature", y = "Coverage", title = args[2], x.text.angle = 90)
dev.off()


