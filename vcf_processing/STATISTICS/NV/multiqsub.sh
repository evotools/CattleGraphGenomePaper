#!/bin/bash

qsub -t 1-9 -tc 9 -N NV_FB NVarTargeted.sh FBsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz FB FB_BQSR
qsub -t 1-9 -tc 9 -hold_jid_ad NV_FB -N NVS_FB NVarSummary.sh FBsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz FB FB_BQSR

qsub -t 1-9 -tc 9 -N NV_GA NVarTargeted.sh GATKsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz GATK GATK_BQSR
qsub -t 1-9 -tc 9 -hold_jid_ad NV_GA -N NVS_GA NVarSummary.sh GATKsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz GATK GATK_BQSR

qsub -t 1-9 -tc 9 -N NV_CA NVarTargeted.sh VGsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz VG VG_CACTUS
qsub -t 1-9 -tc 9 -hold_jid_ad NV_CA -N NVS_CA NVarSummary.sh VGsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz VG VG_CACTUS

qsub -t 1-9 -tc 9 -N NV_LI NVarTargeted.sh LINsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz LINEAR VG_LINEAR
qsub -t 1-9 -tc 9 -hold_jid_ad NV_LI -N NVS_LI NVarSummary.sh LINsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz LINEAR VG_LINEAR

qsub -t 1-9 -tc 9 -N NV_DI NVarTargeted.sh DIVsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz DIVERSITY VG_DIVERSITY
qsub -t 1-9 -tc 9 -hold_jid_ad NV_DI -N NVS_DI NVarSummary.sh DIVsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz DIVERSITY VG_DIVERSITY

