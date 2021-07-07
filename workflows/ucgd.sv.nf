#!/usr/bin/env nextflow

/*
 * UCGD NeoSeq Pipeline
 * Author: Shawn Rynearson
 * Copyright (c) 2020, Utah Center for Genetic Discovery.
*/

// Collect CRAMs channel
Channel.fromPath("$params.crams/*.cram")
    .ifEmpty { error "Path to cram required --crams" }
    .filter( ~/.*.cram$/ )
    .into { cram_input_ch; genotype_in_ch; manta_collect_in }

// Collect CRAMs index channel
Channel.fromPath("$params.crams/*.crai")
    .ifEmpty { error "Path to cram required --crams" }
    .filter( ~/.*.crai$/ )
    .collect()
    .set { manta_collect_index_in }

//--project_path 
if( !params.project_path) { exit 1, "Path to processing directory required." }
println(params.project_path)

// check project name from the command line.
if( !params.project ) { exit 1, "Project name required, not given. (--project [string] )" }

// ------------------------------------------------------- //

process smoove_call {
    label 'smoove'
    tag { $cram }
    
    input:
    each cram from cram_input_ch
    
    output:
    file("*smoove.genotyped.vcf.gz") into smoove_out_ch 
    file("*smoove.genotyped.vcf.gz.csi") into smoove_index_ch 

    shell:
    """
    export TMPDIR="/scratch/local/\$USER/\$SLURM_JOB_ID"
    ln -s ${cram}.crai .
    smoove call --outdir ./ --name ${cram.baseName} -C ${task.exclude} --fasta "${params.reference}" -p 1 --genotype $cram 
    """
}

// ------------------------------------------------------- //

// Collectors.
smoove_out_ch 
    .collect()
    .set { genotype_set_ch }

smoove_index_ch
    .collect()
    .set { index_set_ch }

// ------------------------------------------------------- //

process smoove_merge {
    label 'smoove'
    tag {"${params.project}"}

    input:
    ///val vcf from genotype_set_ch
    file(vcf) from genotype_set_ch
    file(index) from index_set_ch
    
    output:
    file("${params.project}-merged*gz") into smoove_merge_out_ch 

    script:
    all_vcfs = vcf.join(' ')
    """
    export TMPDIR="/scratch/local/\$USER/\$SLURM_JOB_ID"
    smoove merge --name "${params.project}-merged" -f "${params.reference}" --outdir ./ ${all_vcfs}
    """
}

// ------------------------------------------------------- //

process smoove_genotype {
    label 'smoove'
    tag {"${params.project}"}

    input:
    file(merge) from smoove_merge_out_ch
    each cram from genotype_in_ch 
    
    output:
    set file("*smoove.genotyped.vcf.gz"), file("*smoove.genotyped.vcf.gz.csi") into smoove_genotype_out_ch

    shell:
    """
    export TMPDIR="/scratch/local/\$USER/\$SLURM_JOB_ID"
    ln -s ${cram}.crai .
    smoove genotype --duphold -x -p 2 --name ${cram.baseName} --outdir ./ --fasta "${params.reference}" --vcf $merge $cram
    """
}
 
// ------------------------------------------------------- //

// Collect all individual genotyped vcf
smoove_genotype_out_ch
    .collect()
    .set { genotype_collect_in }

// ------------------------------------------------------- //

process smoove_paste {
    label 'smoove'
    tag {"${params.project}"}

    input:
    file(vcfs) from genotype_collect_in
    
    output:
    file("*smoove.square.vcf.gz") into smoove_paste_out_ch

    shell:
    """
    export TMPDIR="/scratch/local/\$USER/\$SLURM_JOB_ID"
    smoove paste --name "${params.project}" *genotyped.vcf.gz
    """
}

// ------------------------------------------------------- //

process smoove_annotate {
    label 'smoove'
    tag {"${params.project}"}
    publishDir "${params.project_path}/VCF/Smoove", mode: 'link', pattern: "*.vcf.gz"
    publishDir "${params.project_path}/VCF/Smoove", mode: 'link', pattern: "*.vcf.gz.tbi"

    input:
    file(vcf) from smoove_paste_out_ch

    output:
    set file("${params.project}.smoove.square.anno.vcf.gz"), file("${params.project}.smoove.square.anno.vcf.gz.tbi") into smoove_out_files
    
    shell:
    """
    export TMPDIR="/scratch/local/\$USER/\$SLURM_JOB_ID"
    smoove annotate --gff ${params.gff3} $vcf | bgzip -c > "${params.project}.smoove.square.anno.vcf.gz"
    tabix -p vcf "${params.project}.smoove.square.anno.vcf.gz"
    """
}

// ------------------------------------------------------- //
// End of smoove
// Start of Manta
// ------------------------------------------------------- //

manta_collect_in
    .collect()
    .set { manta_collect_in }

// ------------------------------------------------------- //

process manta {
    label 'manta'
    tag { "manta" }
    publishDir "${params.manta}", mode: "link", pattern: "*vcf.gz"
    publishDir "${params.manta}", mode: "link", pattern: "*.vcf.gz.tbi"

    input:
    val crams from manta_collect_in
    val indexs from manta_collect_index_in

    output:
    set file("*_candidateSmallIndels.vcf.gz"), file("*_candidateSV.vcf.gz"), file("*_diploidSV.vcf.gz")
    set file("*_candidateSmallIndels.vcf.gz.tbi"), file("*_candidateSV.vcf.gz.tbi"), file("*_diploidSV.vcf.gz.tbi")

    script:
    input_files = crams.join(' --bam ')
    """
    export TMPDIR="/scratch/local/\$USER/\$SLURM_JOB_ID"

    configManta.py --bam ${input_files} --referenceFasta ${params.reference} --runDir ./
    \$PWD/runWorkflow.py

    ## Work post-process. 
    ln results/variants/* .

    ## rename vcf and indexes.
    ls *gz | perl -lane '\$sample_name = "${params.project}_\$_"; system("mv \$_ \$sample_name")'
    ls *gz.tbi | perl -lane '\$sample_name = "${params.project}_\$_"; system("mv \$_ \$sample_name")'

    ## make and move manta stats to reports directory.
    ln results/stats/* "${params.mantastats}/"
    """
}

// ------------------------------------------------------- //
