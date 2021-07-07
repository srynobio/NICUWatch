#!/usr/bin/env nextflow

/*
 * UCGD NeoSeq Pipeline
 * Author: Shawn Rynearson
 * Copyright (c) 2020, Utah Center for Genetic Discovery.
*/

// Initial Channels
Channel.fromPath("$params.crams/*.cram")
    .ifEmpty { error "Path to cram required --crams" }
    .filter( ~/.*.cram$/ )
    .set { cram_input_ch }

// Indiviual VCF files.
Channel.fromPath("$params.vcfs/*.vcf.gz")
    .ifEmpty { error "Path to indivudal vcfs required --vcfs" }
    .filter( ~/.*.vcf.gz$/ )
    .set { vcf_input_ch }

// Indiviual VCF files.
Channel.fromPath("$params.slivar/*.vcf.gz")
    .filter( ~/.*.vcf.gz$/ )
    .set { slivar_input_ch }

// --mosaic_project_id
if( !params.mosaic_project_id) { exit 1, "[--mosaic_project_id] required." }

// ------------------------------------------------------- //

process mosaic_cram_stats {
    tag { $cram }
    label 'mosaic'
    publishDir "${params.mosaicstats}", mode: "link", pattern: "*.json" 

    input:
    each cram from cram_input_ch

    output:
    file("*.mosaic.json") into cram_stats_out_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    mosaic bam -f !{params.reference} -r 38 -b !{params.regionBeds} "!{cram.baseName}" !{cram}
    '''
}

// ------------------------------------------------------- //

Channel.from cram_stats_out_ch 
    .collect()
    .set { cram_stats_collect_ch }

// ------------------------------------------------------- //

process mosaic_vcf_stats {
    tag { $vcf }
    label 'mosaic'
    publishDir "${params.mosaicstats}", mode: "link", pattern: "*.json"
    
    input:
    each vcf from vcf_input_ch

    output:
    file("*.mosaic.json") into vcf_stats_out_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    mosaic vcf -r 38 -b !{params.regionBeds} !{vcf}
    '''
}

// ------------------------------------------------------- //

Channel.from vcf_stats_out_ch 
    .collect()
    .set { vcf_stats_collect_ch }

// ------------------------------------------------------- //

process mosaic_gather_stats {
    tag { $vcf }
    label 'mosaic'
    publishDir "${params.mosaicstats}", mode: "link", pattern: "*.json.gz"
    
    input:
    file bams from cram_stats_collect_ch
    file vcf from vcf_stats_collect_ch

    output:
    file("gathered.mosaic.json.gz") into gather_stats_out_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    mosaic gather ./
    '''
}

// ------------------------------------------------------- //

process mosaic_push_stats {
    tag { $gathered }
    label 'mosaic'

    input:
    file gathered from gather_stats_out_ch

    output:
    val 'pushed' into push_stats_ch
    val 'pushed' into push_slivar_ch 

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    export MOSAIC_URL='**********************'
    export MOSAIC_USERNAME=*****************
    export MOSAIC_PASSWORD='****************'
    mosaic post !{params.mosaic_project_id} !{gathered}
    '''
}

// ------------------------------------------------------- //

process push_pedigree {
    label 'local'
    tag {"$params.project"}

    input:
    val 'pushed' from push_stats_ch

    output:
    file("${params.project}.ped") into ped_out_ch

    script:
    if (params.neotest == 'test')
    """
    NICUWatch make_ped_file -p ${params.project} --test
    NICUWatch mosaic -ma upload_pedigree --pedigree_file "${params.project}.ped" --mosaic_project_id ${params.mosaic_project_id} --test
    """
    else if (params.neotest == 'notest')
    """
    NICUWatch make_ped_file -p ${params.project}
    NICUWatch mosaic -ma upload_pedigree --pedigree_file "${params.project}.ped" --mosaic_project_id ${params.mosaic_project_id}
    """
}

// ------------------------------------------------------- //

process add_slivar_file {
    label 'local'
    tag {"$params.project"}

    input:
    file(slivar_vcf) from slivar_input_ch 
    val 'pushed' from push_slivar_ch 

    when:
    slivar_vcf

    script:
    if (params.neotest == 'test')
    """
    NICUWatch mosaic -ma upload_slivar --mosaic_project_id ${params.mosaic_project_id} --slivar_file ${slivar_vcf} --test
    """
    else if (params.neotest == 'notest')
    """
    NICUWatch mosaic -ma upload_slivar --mosaic_project_id ${params.mosaic_project_id} --slivar_file ${slivar_vcf}
    """
}

// ------------------------------------------------------- //
