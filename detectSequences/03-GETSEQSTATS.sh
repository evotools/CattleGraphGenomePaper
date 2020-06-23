#!/bin/bash
#Grid Engine options (lines prefixed with #$ or #!)
#$ -N SeqStats
#$ -cwd
#$ -l h_rt=47:59:59
#$ -pe sharedmem 8
#$ -R y
#$ -l h_vmem=8.0G
#$ -o ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ./LOGS/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -P roslin_ctlgh

genomes=`python -c "import sys; print(','.join(sys.argv[1:]))" $@`
genome=`python -c "import sys; print(sys.argv[int(sys.argv[-1])])" $@ $SGE_TASK_ID`
echo "All genomes: $genomes"
echo "Target genome: $genome"


if [ ! -e SEQUENCE_STATISTICS ]; then mkdir SEQUENCE_STATISTICS; fi

# Characterise CpG islands
echo "Getting CpG islands in ${genome}"
if [ ! -e SEQUENCE_STATISTICS/${genome} ]; then mkdir SEQUENCE_STATISTICS/${genome}; fi
./SCRIPTS/CpGdetector.sh -f SEQUENCES/${genome}/${genome}.specific.fasta -o SEQUENCE_STATISTICS/${genome} 

# Get Predicted ORF
echo "Getting ORF in ${genome}"
if [ ! -e SEQUENCE_STATISTICS/${genome}/ORF ]; then mkdir SEQUENCE_STATISTICS/${genome}/ORF; fi
./SCRIPTS/ORFfinder.sh SEQUENCES/${genome}/${genome}.specific.fasta SEQUENCE_STATISTICS/${genome}/ORF

# Get motifs enrichments
echo "Getting motifs enrichment in ${genome}"
./SCRIPTS/MotifDetection.sh -b $PWD/SEQUENCES/${genome}/${genome}.bed -g $PWD/SEQUENCES/${genome}/${genome}.fasta -o $PWD/SEQUENCE_STATISTICS/${genome} 

# Characterise repetitive elements
echo "Getting repetitive elements in ${genome}"
if [ ! -e SEQUENCE_STATISTICS/${genome}/repetitiveElements ]; then mkdir SEQUENCE_STATISTICS/${genome}/repetitiveElements; fi
# General count
python SCRIPTS/CountLowerCaseBases.py SEQUENCES/${genome}/${genome}.specific.fasta > SEQUENCE_STATISTICS/${genome}/repetitiveElements/maskedBasesBySeq.tsv
# Repeat masker
python -c "import sys;[sys.stdout.write(line.strip().upper() + '\n') if '>' not in line else sys.stdout.write('>{}\n'.format( '_'.join(line.strip().split(':')[-2:]) )) for line in open(sys.argv[1]) ]" SEQUENCES/${genome}/${genome}.specific.fasta > SEQUENCE_STATISTICS/${genome}/repetitiveElements/TMP.fasta
python ./SCRIPTS/rename.py SEQUENCE_STATISTICS/${genome}/repetitiveElements/TMP.fasta SEQUENCE_STATISTICS/${genome}/repetitiveElements/renamed.fasta && rm SEQUENCE_STATISTICS/${genome}/repetitiveElements/TMP.fasta
./SCRIPTS/RepeatMasker.sh -f $PWD/SEQUENCE_STATISTICS/${genome}/repetitiveElements/renamed.fasta -s cow -n 4 -o $PWD/SEQUENCE_STATISTICS/${genome}/repetitiveElements/ && rm SEQUENCE_STATISTICS/${genome}/repetitiveElements/TMP.fasta


