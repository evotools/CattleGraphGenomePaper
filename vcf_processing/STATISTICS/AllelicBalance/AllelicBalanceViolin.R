library(tidyverse)

mydata = tibble(data.frame( FULL_PATH = list.files(pattern= '*.txt.gz', recursive=T) )) %>% 
    separate(FULL_PATH, sep = "/", into = c("FLDR", "SAMPLE", "FILE"), remove = F) %>% 
    relocate(FULL_PATH, .after = FILE) %>%
    separate(FLDR, sep = "_", into = c("GENOME", "TYPE"), ) 

if (!dir.exists("PLOTS")){ dir.create('PLOTS') }

for ( var_type in unique(mydata$TYPE) ){
    tmp_files = mydata %>% filter( TYPE == var_type )
    for ( n in seq(1, nrow( tmp_files ) ) ){
        tmp_AB = read_tsv(as.character(tmp_files[n,] %>% pull(FULL_PATH))) %>%
                    mutate(ALGORITHM = as.character(tmp_files[n,] %>% pull(GENOME)), SAMPLE = as.character(tmp_files[n,] %>% pull(SAMPLE)) ) %>%
                    relocate( SAMPLE, .before = CHROM ) %>%
                    relocate( ALGORITHM, .before = CHROM ) %>%
                    select( ALGORITHM, SAMPLE, AB )
        if (!exists("allABs")){
            allABs = tmp_AB
        } else {
            allABs = rbind(allABs, tmp_AB)
        }
        cat("Done: ", as.character(tmp_files[n,] %>% pull(FULL_PATH)), "\n")
    }

    p = allABs %>% 
        ggplot(aes(x=SAMPLE, y=AB, fill = ALGORITHM)) +
        geom_violin(position=position_dodge(1)) +
        scale_y_continuous() +
        theme(axis.text.x = element_text(angle=90, hjust=1)) +
        labs(x="Sample", y="Allelic Balance", fill = "Method", title = paste("Allelic balance by method used")) + 
        facet_wrap(~ALGORITHM)
    ggsave(paste("PLOTS/ABviolin_WGS_",var_type,".pdf", sep = ""), plot = p,device = "pdf", height = 9, width = 16)
    rm(allABs)
}