# Make coverage plots
args = commandArgs(T)
library(tidyverse)
library(ggpubr)
library(ggsci)
coverage = read.table(args[1])
colnames(coverage) = c("Coverage", "Value")
pdf(paste(args[3], 'pdf',sep = '.'), width = 12, height = 8)
ggbarplot(data = coverage, x = "Coverage", y = "Value", title = args[2])
dev.off()

pdf(paste(args[3], "trim120", 'pdf',sep = '.'), width = 12, height = 8)
ggbarplot(data = coverage[coverage$Coverage < 120, ], x = "Coverage", y = "Value", title = args[2])
dev.off()