#!/bin/bash
canu -p ndama -d ndama-data \
    useGrid=remote \
    gridEngineThreadsOption="-R y -pe sharedmem THREADS" \ 
    gridEngineMemoryOption="-l h_vmem=MEMORY" \
    gridOptions="-S /bin/bash" \
    genomeSize=2.7g \
    -pacbio-raw Ndama.subreads.fastq \
    gridEngineStageOption="-l h_fsize=50G" \
    stageDirectory="\${TMPDIR}/\${JOB_ID}-\${SGE_TASK_ID}" \
    corMinCoverage=0 \
    corMhapSensitivity=high \
    correctedErrorRate=0.105 \
    corMhapFilterThreshold=0.0000000002 \
    corMhapOptions="--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50" \
    mhapMemory=96g \
    mhapBlockSize=500 \
    ovlMerThreshold=500 \
    corOutCoverage=40
