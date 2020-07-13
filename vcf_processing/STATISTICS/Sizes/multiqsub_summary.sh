#!/bin/bash

#qsub -t 1-9 -tc 9 -N SIZE_FB SummariseSizeTargeted.sh FBsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz FB FB_BQSR
qsub -t 1-9 -tc 9 -hold_jid_ad SIZE_FB -N SIZES_FB SummariseSize.sh FBsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz FB FB_BQSR

#qsub -t 1-9 -tc 9 -N SIZE_GA SummariseSizeTargeted.sh GATKsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz GATK GATK_BQSR
qsub -t 1-9 -tc 9 -hold_jid_ad SIZE_GA -N SIZES_GA SummariseSize.sh GATKsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz GATK GATK_BQSR

#qsub -t 1-9 -tc 9 -N SIZE_CA SummariseSizeTargeted.sh VGsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz VG VG_CACTUS
qsub -t 1-9 -tc 9 -hold_jid_ad SIZE_CA -N SIZES_CA SummariseSize.sh VGsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz VG VG_CACTUS

#qsub -t 1-9 -tc 9 -N SIZE_LI SummariseSizeTargeted.sh LINsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz LINEAR VG_LINEAR
qsub -t 1-9 -tc 9 -hold_jid_ad SIZE_LI -N SIZES_LI SummariseSize.sh LINsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz LINEAR VG_LINEAR

#qsub -t 1-9 -tc 9 -N SIZE_DI SummariseSizeTargeted.sh DIVsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz DIVERSITY VG_DIVERSITY
qsub -t 1-9 -tc 9 -hold_jid_ad SIZE_DI -N SIZES_DI SummariseSize.sh DIVsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz DIVERSITY VG_DIVERSITY

