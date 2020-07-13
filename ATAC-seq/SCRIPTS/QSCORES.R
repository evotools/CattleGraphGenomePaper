# Add Q-score as described in https://github.com/jsh58/Genrich#qvalue
args = commandArgs(T)
data = read.table(args[1])
totlength = as.numeric(args[2])
cat(paste("Input data:", args[1], "\n"))
cat(paste("Genome size:", args[2], "\n"))
data$V9 = p.adjust(10^-data$V8, method = "hochberg", n = totlength)
write.table(data, gsub(".narrowPeak",".withQ.narrowPeak",args[1]), col.names = F, row.names = F, quote = F, sep = "\t")