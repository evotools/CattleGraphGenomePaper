library(tidyverse)

classVec <- function(x, y){
  classes = c(paste("<", x, sep = ''))
  for (cl in seq(as.numeric(x), as.numeric(y) )){
    classes =  c(classes, c(cl))
  }
  classes = c(classes, c(paste(">", y, sep = '')))
  return(classes)
}

mydata = tibble(data.frame( FULL_PATH = list.files(pattern= '*.txt.gz', recursive=T) )) %>% 
    separate(FULL_PATH, sep = "/", into = c("FLDR", "SAMPLE", "FILE"), remove = F) %>% 
    relocate(FULL_PATH, .after = FILE) %>%
    separate(FLDR, sep = "_", into = c("GENOME", "TYPE"), remove = F) 

if (!dir.exists("PLOTS")){ dir.create('PLOTS') }

for ( var_type in unique(mydata$TYPE) ){
    tmp_files = mydata %>% filter( TYPE == var_type )
    for ( n in seq(1, nrow( tmp_files ) ) ){
        tmp_QUAL = read_tsv(as.character(tmp_files[n,] %>% pull(FULL_PATH)), col_types = cols(GROUP = col_character()) ) %>%
                    mutate(ALGORITHM = as.character(tmp_files[n,] %>% pull(FLDR)), SAMPLE = as.character(tmp_files[n,] %>% pull(SAMPLE)) ) %>%
                    relocate( SAMPLE, .before = CHROM ) %>%
                    relocate( ALGORITHM, .before = CHROM ) %>%
                    select( ALGORITHM, GROUP, zQUAL )
        if (!exists("allQUALs")){
            allQUALs = tmp_QUAL
        } else {
            allQUALs = rbind(allQUALs, tmp_QUAL)
        }
        cat("Done: ", as.character(tmp_files[n,] %>% pull(FULL_PATH)), "\n")
    }
    summaries = allQUALs %>% group_by(ALGORITHM, GROUP) %>% summarise(zQUAL=median(zQUAL))
    summaries$GROUP = factor(summaries$GROUP, level = classVec(-60, 60))
    p = summaries %>% 
            ggplot(aes(x=GROUP, y=zQUAL, color = ALGORITHM)) +
            geom_point() + geom_line(aes(group = ALGORITHM)) + geom_hline(aes(yintercept = 0.5), linetype="dashed", colour = "red")  +
            scale_y_continuous(breaks=seq(0, 1, .04), limits = c(0.48, 0.88)) +
            theme(axis.text.x = element_text(angle=90, hjust=1)) +
            labs(x="Length", y="Depth", color = "Software", title = "Depth by variant size") 
    ggsave(paste("PLOTS/zQUAL_WGS_",var_type,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
    rm(allABs)
}