#!/bin/bash

qsub -t 1-9 -tc 9 -N AB_FB AllelicBalanceTargeted.sh FBsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz FB FB_BQSR
qsub -t 1-9 -tc 9 -hold_jid_ad AB_FB -N ABS_FB AllelicBalanceSummary.sh FBsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz FB FB_BQSR

qsub -t 1-9 -tc 9 -N AB_GA AllelicBalanceTargeted.sh GATKsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz GATK GATK_BQSR
qsub -t 1-9 -tc 9 -hold_jid_ad AB_GA -N ABS_GA AllelicBalanceSummary.sh GATKsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz GATK GATK_BQSR

qsub -t 1-9 -tc 9 -N AB_CA AllelicBalanceTargeted.sh VGsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz VG VG_CACTUS
qsub -t 1-9 -tc 9 -hold_jid_ad AB_CA -N ABS_CA AllelicBalanceSummary.sh VGsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz VG VG_CACTUS

qsub -t 1-9 -tc 9 -N AB_LI AllelicBalanceTargeted.sh LINsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz LINEAR VG_LINEAR
qsub -t 1-9 -tc 9 -hold_jid_ad AB_LI -N ABS_LI AllelicBalanceSummary.sh LINsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz LINEAR VG_LINEAR

qsub -t 1-9 -tc 9 -N AB_DI AllelicBalanceTargeted.sh DIVsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz DIVERSITY VG_DIVERSITY
qsub -t 1-9 -tc 9 -hold_jid_ad AB_DI -N ABS_DI AllelicBalanceSummary.sh DIVsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz DIVERSITY VG_DIVERSITY

