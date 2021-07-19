#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Log informations
log.info """\
Non-ref sequence   v 0.5a 
================================
PG                         : $params.pg
Reference genome           : $params.reference
Sequence pool              : $params.genome_pool
Contigs IDs                : $params.contigs
Scaffolds IDs              : $params.scaffolds
Autosomes' repetitiveness  : $params.repetitiveness
Proteins fasta             : $params.proteins
Flanking regions           : $params.flank
Gaps flanking regions      : $params.gap_flank
"""

// Evaluate input files
if (params.pg) { ch_pg = file(params.pg) } else { exit 1, 'Graph pg not specified!' }
if (params.genome_pool) { ch_pool = file(params.genome_pool) } else { exit 1, 'Pooled genomes not specified!' }
if (params.contigs) { ch_ctg = file(params.contigs) } else { exit 1, 'Contigs list not specified!' }
if (params.scaffolds) { ch_scaff = file(params.scaffolds) } else { exit 1, 'Scaffolds list not specified!' }
if (params.repetitiveness) { ch_rep = file(params.repetitiveness) } else { exit 1, 'Scaffolds list not specified!' }
if (params.proteins) { ch_proteins = file(params.proteins) } else { exit 1, 'Scaffolds list not specified!' }

// Run nextflow workflow
include {non_ref_nodes; ref_nodes} from './processes/processes.nf'
include {get_gaps; add_support_vector; add_gap_info} from './processes/processes.nf'
include {combine_regions; frc_filter; label_regions} from './processes/processes.nf'
include {get_repetitiveness; cleanup} from './processes/processes.nf'
include {selfalign; simplify; getfasta; getfasta_flanked} from './processes/processes.nf'
include {make_diamond_db; blastx; bedToFasta} from './processes/processes.nf'
include {abinitio; abinitio_flank} from './processes/processes.nf'
include {filter_abinitio; filter_abinitio_flank} from './processes/processes.nf'
include {consolidate} from './processes/processes.nf'


// Run workflow
workflow {
    // Make diamond database of proteins
    make_diamond_db( ch_proteins )

    // Initial nodes
    non_ref_nodes(ch_pg)
    ref_nodes(ch_pg)

    // Add support vector
    add_support_vector( non_ref_nodes.out )

    // Get gaps in genomes
    get_gaps( ch_pool )

    // Add gaps information
    add_gap_info( add_support_vector.out, get_gaps.out )

    // If specified, run Frc-based filtering
    if (params.frc){
        frc_ch = file(params.frc)
        frc_filter( add_gap_info.out, frc_ch )
        itv_ch = frc_filter.out[0]
    } else {
        itv_ch = add_gap_info.out[0]
    }

    // Combine regions
    combine_regions( itv_ch )

    // Label regions
    label_regions( combine_regions.out, ch_scaff, ch_ctg )

    // Annotate sequences' repetitiveness
    get_repetitiveness( label_regions.out, ch_pool )

    // Multistage filtering of the sequences
    cleanup( get_repetitiveness.out, ch_rep )

    // Generate candidate fasta
    bedToFasta( cleanup.out[0] )

    // Self-alignments to identify redundancies
    selfalign( bedToFasta.out )

    // Simplify redundancies
    simplify( bedToFasta.out, selfalign.out )

    // Make new fastas, plain and flanked
    getfasta( simplify.out[0], ch_pool )
    getfasta_flanked( simplify.out[1], ch_pool )

    // Run blastx for protein identification
    blastx( getfasta.out, make_diamond_db.out )

    // Run augustus for protein identification
    abinitio( getfasta.out[0] )
    filter_abinitio( abinitio.out[1], abinitio.out[2], make_diamond_db.out, ch_proteins )

    // Run augustus for protein identification on flanked sequences
    abinitio_flank( getfasta_flanked.out[0] )
    filter_abinitio_flank( abinitio_flank.out[1], abinitio_flank.out[2], getfasta_flanked.out[1], make_diamond_db.out, ch_proteins )

    // Consolidate results
    consolidate( filter_abinitio.out, filter_abinitio_flank.out, blastx.out )
}

