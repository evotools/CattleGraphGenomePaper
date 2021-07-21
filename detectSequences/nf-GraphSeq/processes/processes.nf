

process non_ref_nodes {
    label "large"
    tag "non_ref_nodes"
    publishDir "${params.outfolder}", mode: 'symlink'

    input:
    path pg

    output:
    path "non_ref_nodes.bed"

    script:
    """
    01A-NonRefNodes -i ${pg} -o non_ref_nodes -r ${params.reference}
    """
}


process ref_nodes {
    label "large"
    tag "non_ref_nodes"
    publishDir "${params.outfolder}/00-initial_nodes", mode: 'symlink'

    input:
    path pg

    output:
    path "ref_nodes.bed"

    script:
    """
    01B-RefNodes -i ${pg} -o ref_nodes -r ${params.reference}
    """
}


process add_support_vector {
    label "medium"
    tag "supp_vec"

    input:
    path non_ref_nodes

    output:
    path "non_ref_nodes.supvec.bed"

    script:
    """
    02-ClassifyNodes -i ${non_ref_nodes} -o non_ref_nodes.supvec
    """
}

process get_gaps {
    label "medium"
    tag "get_gaps"

    input:
    path pooled_genomes

    output:
    path "gaps.bed"

    script:
    """
    faToTwoBit ${pooled_genomes} ${pooled_genomes.simpleName}.2bit
    twoBitInfo -nBed ${pooled_genomes.simpleName}.2bit stdout | \
        awk -v var=${params.gap_flanks} 'BEGIN{OFS="\t"}; \$2-var < 0{print \$1,"0",\$3+var}; \$2-var >= 0{print \$1,\$2-var,\$3+var}' | \
        bedSort stdin gaps.bed
    """
}


process add_gap_info {
    label "small_long"
    tag "supp_vec"
    publishDir "${params.outfolder}/01-filter_nmers", mode: 'symlink'

    input:
    path with_support
    path gaps

    output:
    path "*.labeled.bed"
    path "*.pdf"
    path "*.log"

    script:
    """
    bname=`basename -s '.bed' $with_support`
    02-ClassifyNodes -i ${with_support} -o \$bname
    bedtools intersect -a \${bname}.noNmers.bed -b ${gaps} -v | awk 'BEGIN{OFS="\t"};{print \$0, "0"}' > \${bname}.noOverlaps.bed
    bedtools intersect -a \${bname}.noNmers.bed -b ${gaps} -u | awk 'BEGIN{OFS="\t"};{print \$0, "1"}' > \${bname}.Overlaps.bed
    cat \${bname}.noOverlaps.bed \${bname}.Overlaps.bed | bedSort stdin \${bname}.labeled.bed && \
        rm \${bname}.noOverlaps.bed \${bname}.Overlaps.bed
    04B-nodesStats \${bname}.labeled.bed upsetPlot.no_Nmers > upsetPlot.no_Nmers.log
    """
}

process frc_filter {
    label "small_long"
    tag "frc_filt"
    publishDir "${params.outfolder}/02-filtered_frc", mode: 'symlink'

    input:
    path intervals
    path frc

    output:
    path "*.no_frc.bed"
    path "upsetPlot.no_Nmers.no_frc*"

    script:
    """
    bname=`basename -s '.bed' $intervals`
    bedtools intersect -a $intervals -b $frc -v | awk 'BEGIN{OFS="\t"};{print \$0}' > \${bname}.no_frc.bed    
    04B-nodesStats \${bname}.no_frc.bed upsetPlot.no_Nmers.no_frc > upsetPlot.no_Nmers.no_frc.log
    """
}

process combine_regions {
    label "medium"
    tag "combine_reg"
    publishDir "${params.outfolder}/03-combined", mode: 'symlink'

    input:
    path intervals

    output:
    path "*.lengths.merged.bed"

    script:
    """
    oname=`basename -s '.bed' ${intervals}`

    # Example of header
    # Add region lengths
    awk 'BEGIN{OFS="\t"};{print \$0, \$3-\$2}' \$intervals > \${oname}.lengths.bed
    # Create output header
    echo | awk 'BEGIN{OFS="\t"};{print "SEQID","BPI","BPE","NODES","N_NODES","STRANDS","SEQS","N_CLOSE_TO_GAPS","NODES_LENGTH"}' > \${oname}.lengths.merged.bed
    # Perform region merging
    bedtools merge -c 4,4,5,6,8,9 -o collapse,count,collapse,collapse,sum,sum -d 5 -i \${oname}.lengths.bed >> \${oname}.lengths.merged.bed
    """
}

process label_regions {
    label "small"
    tag "label_reg"

    input:
    path intervals
    path scaffolds
    path contigs
    
    output:
    path "*.seqtype.bed"

    script:
    """
    bname=`basename -s '.bed' $intervals`
    06B-ClassifyRegions -i ${intervals} -o \${bname} -c ${contigs} -s ${scaffolds} 
    """
}

process get_repetitiveness {
    label "medium"
    tag "add_rept"
    publishDir "${params.outfolder}/04-labelled", mode: 'symlink'

    input:
    path intervals
    path sequences
    
    output:
    path "*.masked.bed"

    script:
    """
    bname=`basename -s '.bed' $intervals`

    # Extract upper-case sequences
    bedtools getfasta -fi ${sequences} -bed ${intervals} -tab 2> getfasta.err | \
        python -c "import sys; [sys.stdout.write(f'{line.strip().split()[0]}\t{sum([int(c.islower()) for c in line.strip().split()[1]])}\t{len(line.strip().split()[1])}\t{line.strip().split()[1]}\n') for line in sys.stdin]" > regions.txt
    # Add info about masked bases and sequence in the region
    07B-combine ${intervals} regions.txt > \${bname}.masked.bed
    """
}

process cleanup {
    label "medium"
    tag "cleanup"
    publishDir "${params.outfolder}/05-cleanup", mode: 'symlink'

    input:
    path intervals
    path repetitiveness

    output:
    path "*.long.novel.noTelomere.noFlankGaps.lowrep.candidate.bed"
    path "*.pdf"
    path "*.log"

    script:
    """
    name=`basename -s '.bed' $intervals`
    awk 'NR==1 {print}; NR>1 && \$11~"LONG"{print}' ${intervals} > \${name}.long.bed
    04B-nodesStats \${name}.long.bed upsetPlot.no_Nmers.no_frc.long > upsetPlot.no_Nmers.no_frc.long.log

    awk -v val=${params.novelty_cutoff} 'NR==1{print};NR>1 && \$9/\$10>val {print}' \${name}.long.bed > \${name}.long.novel.bed
    04B-nodesStats \${name}.long.novel.bed upsetPlot.no_Nmers.no_frc.long.novel > upsetPlot.no_Nmers.no_frc.long.novel.log

    awk 'NR==1{print};NR>1 && \$11!~"TELOMER"{print}' \${name}.long.novel.bed > \${name}.long.novel.noTelomere.bed
    04B-nodesStats \${name}.long.novel.noTelomere.bed upsetPlot.no_Nmers.no_frc.long.noTelomere.novel > upsetPlot.no_Nmers.no_frc.long.novel.noTelomere.log

    awk 'NR==1{print};NR>1 && \$11!~"FLANK"{print}' \${name}.long.novel.noTelomere.bed > \${name}.long.novel.noTelomere.noFlankGaps.bed
    04B-nodesStats \${name}.long.novel.noTelomere.noFlankGaps.bed upsetPlot.no_Nmers.no_frc.long.novel.noTelomere.noFlankGaps > upsetPlot.no_Nmers.no_frc.long.novel.noTelomere.noFlankGaps.log

    echo "" | awk 'BEGIN{OFS="\t"}; {print "SEQID", "BPI", "BPE", "NODES", "N_NODES", "STRANDS", "NODE_SEQUENCE", "N_CLOSE_TO_GAPS", "NODES_LENGTH", "REGION_SIZE", "CLASSIFICATION", "N_MASKED", "N_NT", "RATIO_MASKED", "ZSCORE", "PVAL", "SEQUENCE"}' > \${name}.long.novel.noTelomere.noFlankGaps.lowrep.candidate.bed
    09A-FilterRepetitive \${name}.long.novel.noTelomere.noFlankGaps.lowrep.candidate.bed \
        ${repetitiveness} \${name}.long.novel.noTelomere.noFlankGaps.lowrep.candidate.bed
    04B-nodesStats \${name}.long.novel.noTelomere.noFlankGaps.lowrep.candidate.bed upsetPlot.no_Nmers.no_frc.long.novel.noTelomere.noFlankGaps.lowrep > upsetPlot.no_Nmers.no_frc.long.novel.noTelomere.noFlankGaps.lowrep.log
    """
}

process bedToFasta {
    tag "bed2fa"
    label "small"

    input:
    path intervals
    
    output:
    tuple path("candidate.fa"), path("candidate.fa.fai")

    script:
    """
    python -c 'import sys; [sys.stdout.write( f">{line.strip().split()[0]}_{line.strip().split()[1]}-{line.strip().split()[2]}\n{line.strip().split()[-1]}\n" ) for line in open(sys.argv[1]) if "SEQID" not in line]' ${interval} > candidate.fa
    samtools faidx candidate.fa
    """
}

process selfalign {
    tag "selfalign"
    label "medium"

    input:
    tuple path(candidates), path(fai)

    output:
    path "alignments.blasttab"

    script:
    """
    minimap2 -x asm5 -t ${task.cpus} --cs=long ${candidate} ${candidate} | \
        paftools.js view -f maf - > alignments.maf
    maf_path=`which maf-convert`
    cp \${maf_path} ./maf_convert && 2to3 -n -o ${PWD} ./maf_convert && chmod a+x ./maf_convert 
    ./maf_convert blasttab alignments.maf > alignments.tmp
    08B-AddScoresToBlast6 alignments.maf alignments.tmp > alignments.blasttab
    """
}

process simplify {
    tag "simplify"
    label "medium"
    publishDir "${params.outfolder}/06-clumped", mode: 'symlink'


    input:
    tuple path(fasta), path(fai)
    path alignments

    output:
    path "candidate.clump.bed"
    path "candidate.clump.flank${params.flank}.bed"

    script:
    """
    08C-DetectDuplicateContigs ${alignments} ${fai} candidate.clump.txt
    09D-faiToBed candidate.clump.txt > candidate.clump.bed
    if [ ${params.flanks} ]; then 
        awk -v var=${params.flanks} '\$2-var <=0 {print \$1, "0", \$3 + var}; \$2-var >0 {print \$1, \$2-var, \$3 + var}' > candidate.clump.flank${params.flank}.bed
    fi 
    """
}

process getfasta {
    tag "getfasta"
    label "small"
    publishDir "${params.outfolder}/06-clumped", mode: 'symlink'


    input:
    path intervals
    path sequences

    output:
    path "*.fa"

    script:
    """
    name=`basename -s '.bed' $intervals`
    bedtools getfasta -fi ${sequences} -bed ${intervals} -fo \${name}.fa
    """
}

process getfasta_flanked {
    tag "getfasta"
    label "small"
    publishDir "${params.outfolder}/06-clumped", mode: 'symlink'


    input:
    path intervals
    path sequences

    output:
    path "*.flank.fa"
    path "*.flank.bed"

    script:
    """
    name=`basename -s '.bed' $intervals`
    11F-addFlanks ${intervals} ${params.flanks} > \${name}.flank.bed
    bedtools getfasta -fi ${sequences} -bed \${name}.flank.bed -fo \${name}.flank.fa
    """
}

process make_diamond_db {
    tag "makedb"
    label "large"

    input:
    path fasta

    output:
    path "mydb.dmnd"

    script:
    """
    diamond makedb -d mydb --in ${fasta}
    """

}

process blastx {
    tag "diamond_bx"
    label "large"
    publishDir "${params.outfolder}/07-predict/07A-blastx", mode: 'symlink'

    input:
    path intervals
    path protdb

    output:
    path "candidates.clump.bx.aligned.ev1e-10.id80.cov80.tab"

    script:
    """
    diamond blastx -d ${protdb} \
        -p ${tasks.cpus} -k 1 --more-sensitive \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore scovhsp \
        --evalue 1e-10 --id 80 -q ${intervals}.fa> candidates.clump.bx.aligned.ev1e-10.id80.tab
    wc -l candidates.clump.bx.aligned.ev1e-10.id80.tab

    # Filter low-covered genes
    awk '\$NF>80 {print}' candidates.clump.bx.aligned.ev1e-10.id80.tab > candidates.clump.bx.aligned.ev1e-10.id80.cov80.tab
    wc -l candidates.clump.bx.aligned.ev1e-10.id80.cov80.tab
    """
}

process abinitio {
    tag "augustus"
    label "large"
    publishDir "${params.outfolder}/07-predict/07B-augustus", mode: 'symlink'

    input:
    path fasta

    output:
    path "*abInitio.gff"
    path "*abInitio.complete.aa"
    path "*abInitio.complete.gff"
    path "*abInitio.incomplete.aa"
    path "*abInitio.incomplete.gff"

    script:
    """
    name=`basename -s '.fa' ${fasta}`
    augustus --species=human ${fasta} > \${name}.abInitio.gff
    11B-ExtractGenes -i \${name}.abInitio.gff -o \${name}.abInitio
    """
}

process abinitio_flank {
    tag "augustus_flank"
    label "large"
    publishDir "${params.outfolder}/07-predict/07B-augustus", mode: 'symlink'

    input:
    path fasta

    output:
    path "*abInitio.flank.gff"
    path "*abInitio.flank.complete.aa"
    path "*abInitio.flank.complete.gff"
    path "*abInitio.flank.incomplete.aa"
    path "*abInitio.flank.incomplete.gff"

    script:
    """
    name=`basename -s '.fa' ${fasta}`
    augustus --species=human ${fasta} > \${name}.abInitio.flank.gff
    11B-ExtractGenes -i \${name}.abInitio.flank.gff -o \${name}.abInitio.flank
    """
}

process filter_abinitio {
    tag "filter_abinit"
    label "large"
    publishDir "${params.outfolder}/07-predict/07B-augustus", mode: 'symlink'

    input:
    path fastaa
    path gff
    path protdb
    path proteins

    output:
    path "*.aligned.eval1e-10.id80.cov80.withNames.withSeq.tab"

    script:
    """
    name=`basename -s '.aa' $fastaa`
    # Align to protein db provided
    diamond blastp -d ${protdb} \
        -p 4 -k 1 --more-sensitive \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore scovhsp \
        --evalue 1e-10 --id 80 -q ${fastaa} > \${name}.aligned.eval1e-10.id80.tab
    wc -l \${name}.aligned.eval1e-10.id80.tab

    # Filter low-covered genes
    awk '\$NF > 80 {print}' \${name}.aligned.eval1e-10.id80.tab > \${name}.aligned.eval1e-10.id80.cov80.tab
    11C-protNamesFromHeader \${name}.aligned.eval1e-10.id80.cov80.tab ${proteins} > \${name}.aligned.eval1e-10.id80.cov80.withNames.tab
    11D-getSequenceNames \${name}.aligned.eval1e-10.id80.cov80.withNames.tab ${gff} > \${name}.aligned.eval1e-10.id80.cov80.withNames.withSeq.tab
    """
}

process filter_abinitio_flank {
    tag "filter_abinit_fl"
    label "large"
    publishDir "${params.outfolder}/07-predict/07C-flanked_augustus", mode: 'symlink'

    input:
    path fastaa
    path gff
    path bed
    path protdb
    path proteins

    output:
    path "*.aligned.eval1e-10.id80.cov80.withNames.withSeq.unflank.tab"

    script:
    """
    name=`basename -s '.aa' $fastaa`
    # Align to protein db provided
    diamond blastp -d ${protdb} \
        -p 4 -k 1 --more-sensitive \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore scovhsp \
        --evalue 1e-10 --id 80 -q ${fastaa} > \${name}.aligned.eval1e-10.id80.tab
    wc -l \${name}.aligned.eval1e-10.id80.tab

    # Filter low-covered genes
    awk '\$NF > 80 {print}' \${name}.aligned.eval1e-10.id80.tab > \${name}.aligned.eval1e-10.id80.cov80.tab
    11C-protNamesFromHeader \${name}.aligned.eval1e-10.id80.cov80.tab ${proteins} > \${name}.aligned.eval1e-10.id80.cov80.withNames.tab
    11D-getSequenceNames \${name}.aligned.eval1e-10.id80.cov80.withNames.tab ${gff} > \${name}.aligned.eval1e-10.id80.cov80.withNames.withSeq.tab
    11G-removeFlanks \${name}.aligned.eval1e-10.id80.cov80.withNames.withSeq.tab ${bed} > \${name}.aligned.eval1e-10.id80.cov80.withNames.withSeq.unflank.tab
    """
}


process consolidate {
    tag "consilidate"
    label "large"
    publishDir "${params.outfolder}/08-consolidate", mode: 'symlink'

    input:
    path blastx
    path augustus
    path augustus_flank



    script:
    """

    """
}