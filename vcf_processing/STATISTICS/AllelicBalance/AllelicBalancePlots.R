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
if (!dir.exists("SUMMARIES")){ dir.create('SUMMARIES') }

for ( var_type in unique(mydata$TYPE) ){
    tmp_files = mydata %>% filter( TYPE == var_type )
    dat<-list()
    for ( n in seq(1, nrow( tmp_files ) ) ){
        dat[[n]] = read_tsv(as.character(tmp_files[n,] %>% pull(FULL_PATH)), col_types = cols(GROUP = col_character(), AB = col_double()) ) %>%
                    mutate(ALGORITHM = as.character(tmp_files[n,] %>% pull(FLDR)), SAMPLE = as.character(tmp_files[n,] %>% pull(SAMPLE)) ) %>%
                    select( ALGORITHM, GROUP, AB )
        cat("Done: ", as.character(tmp_files[n,] %>% pull(FULL_PATH)), "\n")
    }
    dat<-bind_rows(dat, .id="Index")
    summaries = dat %>% group_by(ALGORITHM, GROUP) %>% summarise(AB=median(AB))
    summaries$GROUP = factor(summaries$GROUP, level = classVec(-60, 60))
    write.csv2(summaries, paste("SUMMARIES/AB_WGS",var_type,".csv", sep = ""))
    p = summaries %>% 
            ggplot(aes(x=GROUP, y=AB, color = ALGORITHM)) +
            geom_point() + geom_line(aes(group = ALGORITHM)) + geom_hline(aes(yintercept = 0.5), linetype="dashed", colour = "red")  +
            scale_y_continuous(breaks=seq(0, 1, .04), limits = c(0.48, 0.88)) +
            theme(axis.text.x = element_text(angle=90, hjust=1)) +
            theme_classic(base_size = 18) +
            labs(x="Length", y="Allelic balance", color = "Software", title = "Allelic balance by variant size") +
            ggpubr::rotate_x_text(angle = 90)
    ggsave(paste("PLOTS/AB_WGS",var_type,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 18)
    rm(dat)
    gc()
}