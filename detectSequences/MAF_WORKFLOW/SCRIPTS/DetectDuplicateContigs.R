#!/usr/bin/env Rscript
## Get duplicate sequences
# Use minimap2 alignments to detect which novel contigs are repeated
options(stringsAsFactors = F, scipen = 999)
library(tidyverse)
library(reshape2)
args = commandArgs(T)

alignments = read.table(args[1])
contigs = read.table(args[2])
# Add concordance vector
contigs$CONCORDANT = ""
contigs = contigs %>% arrange(desc(V2))
cat("Processing ", nrow(alignments), " for ", nrow(contigs), " contigs.\n" )

# Add column names
colnames(alignments) = c("qseqid", "tseqid", "ident", "length", "mismatch", "gapopen", 
                         "qstart", "qend", "tstart", "tend", "evalue", "bitscore")
# Extract non-self-alignments
alignments_filt1 = alignments %>% 
                    filter(qseqid != tseqid) %>%
                    group_by(qseqid, tseqid) %>%
                    arrange(desc(ident) )
cat("Kept ", nrow(alignments_filt1), " alignments between contigs.\n" )

# Add sequence length
alignments_filt2 = merge(alignments_filt1, data.frame(qseqid = contigs[,1], qseqlen = contigs[,2]), by = "qseqid", all.x = T)
alignments_filt2 = merge(alignments_filt2, data.frame(tseqid = contigs[,1], tseqlen = contigs[,2]), by = "tseqid", all.x = T)

# Extract highly similar contigs (>99% identical)
alignments_filt3 = alignments_filt2 %>% filter(ident > 99)
cat("Kept ", nrow(alignments_filt1), " alignments with identity >99%.\n" )

# Remove if length of alignment is less than 95% of the shortest contigs
alignments_filt4 = alignments_filt3 %>% mutate( lratio = length / pmin(qseqlen, tseqlen) ) %>% filter(lratio > 0.95)
cat("Kept ", nrow(alignments_filt1), " alignments spanning 95% of the shortest contig .\n" )

# Group by values, keeping the largest contig
cat("Initial number of contigs (n)                : ", nrow(contigs), "\n")
cat("Initial number of contigs (bp)               : ", sum(contigs$V2), "\n")
outputs = contigs[!contigs$V1 %in% unique(c(alignments_filt4$tseqid, alignments_filt4$qseqid)),]
cat("Contigs with no/low identity alignments (n)  : ", nrow(outputs), "\n")
cat("Contigs with no/low identity alignments (bp) : ", sum(outputs$V2), "\n")
toprocess = contigs[contigs$V1 %in% unique(c(alignments_filt4$tseqid, alignments_filt4$qseqid)),] %>%
  arrange(desc(V2))
while (nrow(toprocess) > 0){
  ctgname = toprocess[1, 1]
  tgtalign = alignments_filt4 %>% filter(tseqid == ctgname | qseqid == ctgname)
  if ( nrow(tgtalign) > 0 ){
    outputs = rbind(outputs, toprocess[1, ])
    toadd = unique( c(tgtalign[,"qseqid"], tgtalign[,"tseqid"]) [!grepl(ctgname, c(tgtalign[,"qseqid"], tgtalign[,"tseqid"]))] )
    outputs[outputs$V1 == ctgname, "CONCORDANT"] = paste(toadd, c(outputs[outputs$V1 == ctgname, "CONCORDANT"]), sep = '', collapse = ",")
    toprocess = toprocess[-which( toprocess$V1 %in% c(tgtalign$tseqid, tgtalign$qseqid)), ]
    cat(nrow(toprocess), ' contigs left to check\n')
  }
}
cat("Final number of contigs (n)                  : ", nrow(outputs), "\n")
cat("Final number of contigs (bp)                 : ", sum(outputs$V2), "\n")

write.table(outputs, args[3], sep = "\t", quote = F, col.names = F, row.names = F)
