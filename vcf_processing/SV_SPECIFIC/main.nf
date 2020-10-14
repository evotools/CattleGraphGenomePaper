#!/usr/bin/env nextflow
 
/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options.
 * Karyotype can be specified as ranges (e.g. 1-22), single different 
 * chromosomes can be added after comma (e.g. 1-22,X,Y,Mt).
 */
params.faifile = "faifile.fai"
params.ntest = 10000
params.omregions = "OMregion.bed"
params.lengths = "lengths.txt"
params.positives = 1000
params.outdir = "OUTPUT"



process createRegion{
    tag "makeregions"

    input:
    each x from 1..params.ntest

    output:
    tuple x, "myregion.${x}.bed" into regions_ch

    script:
    $/
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
    args = commandArgs(trailingOnly=TRUE)
    options(scipen = 999)

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
    myfaifile = read.table("${params.faifile}", h=F) %>%
        filter(V1 %in% c(1:29))

    lengths = as.numeric(abs(read.table("${params.lengths}", h=F)[,1]))

    regions = generateRandomPos(myfaifile[,'V1'], myfaifile[,'V2'], lengths)
    regions = as.data.frame(cbind(regions, paste("REGION", c(1: nrow(regions)), sep = "")))
    write.table(regions, "myregion.${x}.bed", sep = "\t", row.names = F, col.names = F, quote = F)

    /$
}

process intersect {
    tag "intersect"

    input:
    tuple x, region from regions_ch

    output:
    file "intersected.${x}.bed" into intersected_ch

    script:
    """
    bedtools intersect -a ${region} -b ${params.omregions} | sort | uniq > intersected.${x}.bed
    """
}

process collect{
    tag collect
    publishDir "${params.outdir}/intersections", overwrite: true, mode: "copy"

    input:
    file intersected from intersected_ch.collect()

    output:
    file "results.txt" into numbers_ch

    script:
    """
    for fname in ${intersected}; do
        wc -l \$fname | awk '{print \$1}'
    done > results.txt
    """
}

process zscores{
    tag "zscores"
    publishDir "${params.outdir}/intersections", overwrite: true, mode: "copy"

    input:
    file result from numbers_ch

    output:
    file "dist.pdf" into plot_ch
    file "pvalue.txt" into pvalues_ch

    script:
    $/
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
    suppressPackageStartupMessages(library(ggpubr, quietly = TRUE))
    args = commandArgs(trailingOnly=TRUE)
    options(scipen = 999)

    positive = as.numeric("${params.positives}")
    print(paste("Positive:", positive))


    # Get the fai file
    myfaifile = read.table("${params.faifile}", h=F) %>%
    filter(V1 %in% c(1:29))

    lengths = as.numeric(abs(read.table("${params.lengths}", h=F)[,1]))
    found = c(read.table("${result}", h=F)[,1])
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
                        VALUES = c(positive, length(lengths), "${params.ntest}", Perm_mean, Perm_sd, z, z_P) )

    final_p = write.table(results, "pvalue.txt")
    /$
}