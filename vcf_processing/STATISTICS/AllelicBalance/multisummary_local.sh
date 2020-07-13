#$/bin/bash
./AllelicBalanceSummary_local.sh GATKsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz GATK GATK_BQSR
./AllelicBalanceSummary_local.sh VGsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz VG VG_CACTUS
./AllelicBalanceSummary_local.sh LINsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz LINEAR VG_LINEAR
./AllelicBalanceSummary_local.sh DIVsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz DIVERSITY VG_DIVERSITY
./AllelicBalanceSummary_local.sh FBsamples.txt ../../NewGraph_03112019/VCF/CACTUS_VARIANTS/CACTUSVAR.vcf.gz FB FB_BQSR
