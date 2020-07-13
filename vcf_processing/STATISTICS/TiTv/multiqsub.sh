#!/bin/bash

qsub -t 1-9 -tc 9 -N TT_FB TiTvTargeted.sh FBsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz FB FB_BQSR
#qsub -t 1-9 -tc 9 -hold_jid_ad TT_FB -N TTS_FB TTarSummary.sh FBsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz FB FB_BQSR

qsub -t 1-9 -tc 9 -N TT_GA TiTvTargeted.sh GATKsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz GATK GATK_BQSR
#qsub -t 1-9 -tc 9 -hold_jid_ad TT_GA -N TTS_GA TTarSummary.sh GATKsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz GATK GATK_BQSR

qsub -t 1-9 -tc 9 -N TT_CA TiTvTargeted.sh VGsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz VG VG_CACTUS
#qsub -t 1-9 -tc 9 -hold_jid_ad TT_CA -N TTS_CA TTarSummary.sh VGsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz VG VG_CACTUS

qsub -t 1-9 -tc 9 -N TT_LI TiTvTargeted.sh LINsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz LINEAR VG_LINEAR
#qsub -t 1-9 -tc 9 -hold_jid_ad TT_LI -N TTS_LI TTarSummary.sh LINsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz LINEAR VG_LINEAR

qsub -t 1-9 -tc 9 -N TT_DI TiTvTargeted.sh DIVsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz DIVERSITY VG_DIVERSITY
#qsub -t 1-9 -tc 9 -hold_jid_ad TT_DI -N TTS_DI TTarSummary.sh DIVsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz DIVERSITY VG_DIVERSITY

