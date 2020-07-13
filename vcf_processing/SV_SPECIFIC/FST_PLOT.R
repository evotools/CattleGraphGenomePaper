library(ggplot2)
library(ggrepel)
library(qqman)
library(magrittr)
library(dplyr)

fst = read.table("FSTVAL.85.windowed.weir.fst.gz", h=T)
fst$CHROM = factor(fst$CHROM, paste("hereford", c(1:29), sep = "."))
fst_filt = fst[!is.nan(fst$WEIGHTED_FST),]
fst_filt$VAR = paste(fst_filt$CHROM, fst_filt$BIN_START, sep = "_")

varOfInterest = c(fst_filt[fst_filt$WEIGHTED_FST > 0.8, "VAR"])
fst_filt$CHROM = as.numeric(gsub("hereford.", "", fst_filt$CHROM))
fst_filt$CHROM = factor(fst_filt$CHROM, levels = c(1:29))

don <- fst_filt %>% 
  
  # Compute chromosome size
  group_by(CHROM) %>% 
  summarise(chr_len=max(BIN_START)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(fst_filt, ., by=c("CHROM"="CHROM")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHROM, BIN_START) %>%
  mutate( BPcum=BIN_START+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(VAR %in% varOfInterest, "yes", "no")) %>%
  mutate( is_annotate=ifelse(WEIGHTED_FST>0.8, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(CHROM) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
ggplot(don, aes(x=BPcum, y=WEIGHTED_FST)) +
  
  # Show all points
  geom_point( aes(color=CHROM), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 29 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=VAR), size=2) +
  
  # Change label name
  labs(x = "Chromosome") +
  
  # Change y limits
  ylim(0, 1) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


