
# Assign chromosome number to scaffolds
options(stringsAsFactors = F)
library(tidyverse)
chromosomes = c(c(1:29), c("X", "Y", "MT"))

# Plot alignments function
plotChrom <- function(x, alignments = NULL, initialV = 0, endingV = NULL){
  if (is.null(endingV)){ endingV = max(alignments$tl) + 1 }
  if (is.null(alignments)){ return(0) }
  chr = alignments %>% filter(tn == x & (alignments$tbpi >= initialV | alignments$tbpe >= initialV) & (alignments$tbpe <= endingV | alignments$tbpi <= endingV) ) %>% arrange(ql)
  ctgs = data.frame(qn = sort(unique(chr$qn)), y1 = seq(1, length(unique(chr$qn))), y2= seq(1, length(unique(chr$qn))) )
  chr = merge(chr, ctgs, by = "qn", all.x = T) %>% arrange(ql)
  p = ggplot(chr) + 
    geom_segment(aes(x = tbpi, y = y1, xend = tbpe, yend = y2, colour = strand, label = qn)) +
    labs(x = paste("Chromosome", x), y = "Scaffolds") #+ 
    #ggrepel::geom_text_repel(x = chr$tbpi, y = chr$y1, label = paste(chr$qn, chr$qbpi, chr$qbpe, sep = "-"))
  p
}

grepAlignments <- function(alignment, tgt=NULL, qry=NULL){
  if ((is.null(tgt) & is.null(qry)) | is.null(alignment)){return(1)}
  if (is.null(qry)){
    return(alignment %>% filter(tn == tgt))
  } else if (is.null(tgt)){
    return(alignment %>% filter(qn == qry))
  } else {
    return(alignment %>% filter(tn == tgt & qn == qry))
  }
  
}

# Start processing data
alignments = read.table("alignments_cut.paf", h=F, sep = "\t", )
colnames(alignments) = c("qn","ql","qbpi","qbpe","strand",
                         "tn","tl","tbpi","tbpe","nmatches","alLength", "MQ")
summaryData = alignments %>% group_by(qn,tn) %>% summarise(aligned = sum(alLength), tLength = mean(tl), qLength = mean(ql)) %>% mutate(ratio = aligned / tLength, mappedRatio = aligned / qLength)
chromosomeAssignments = summaryData[summaryData$ratio>0.80,] %>% filter(tn %in% chromosomes)
write.csv(chromosomeAssignments, "chromosomeCodes.csv")

scaffolds = alignments %>% group_by(qn) %>% summarise(qLength = mean(ql)) %>% filter(qLength > 20000000)
lft = summaryData[!summaryData$qn %in% chromosomeAssignments$qn,] %>% filter(qLength > 20000000) %>% filter(mappedRatio > .10)
write.csv(lft, "largeScaffolds.csv")
