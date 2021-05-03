#!/bin/bash
/exports/cmvm/eddie/eb/groups/prendergast_grp/Andrea/canu-1.8/Linux-amd64/bin/
canu \
    -p ankole_canu \
    -d ankole_canu \
    useGrid=remote \
    gridEngineThreadsOption="-R y -pe sharedmem THREADS" \
    gridEngineMemoryOption="-l h_vmem=MEMORY" \
    gridOptions="-S /bin/bash -P roslin_ctlgh" \
    genomeSize=2.7g \
    -pacbio-raw /exports/cmvm/eddie/eb/groups/prendergast_roslin/Andrea/AnkoleAssembly/MappingReads/Ankole_UG833.fasta.gz \
    gridEngineStageOption="-l h_fsize=50G" \
    stageDirectory="\${TMPDIR}/\${JOB_ID}-\${SGE_TASK_ID}" \
    corMhapSensitivity=high \
    corMhapFilterThreshold=0.0000000002 \
    corMhapOptions="--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50" \
    mhapMemory=96g \
    mhapBlockSize=500 \
    ovlMerDistinct=0.975 \
