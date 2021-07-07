#!/usr/bin/env nextflow

/*
 * UCGD Pipeline
 * Author: Shawn Rynearson
 * Copyright (c) 2020, UStar Center for Genetic Discovery.
*/

// Log info.
log.info "================================="
log.info "     UCGD NICU MULTI FASTQ       "
log.info "================================="
log.info "Version: $params.version"
log.info ""
log.info "~~ Using the following settings ~~"
log.info "UCGD Individual : $params.accession"
log.info ""
log.info "================================="

// Example:
// nextflow run ucgd.nicu.nf --accession "123456TEST" --proband -c ucgd.master.config

// Is this sample a proband.
proband = false
if ( params.proband ) { proband = true }

// check accession from the command line.
if( !params.accession ) { exit 1, "ARUP accession not given." }

// Default ARUP NICU location.
processing = "$PROCESSING/ARUP_NICU/$params.accession/"

// fastq channel
Channel
    .fromPath("$processing/*fastq.gz")
    .collect()
    .into { check_bgzf_in_ch; make_manifest_file_ch }
    ///////.into { fastq_in_ch; check_bgzf_in_ch; make_manifest_file_ch }

// ------------------------------------- //

/* Needed DB interactions:
*/

// ------------------------------------- //

process check_bgzf {
    tag {"$params.accession"}
    label 'localterm'

    input:
    each fastq from check_bgzf_in_ch
 
    shell:
    '''
    check_bgzf.pl !{fastq}
    '''
}

// ------------------------------------- //

process create_manifest_file {
    tag {"$params.accession"}
    label 'localterm'

    input:
    file fastq from make_manifest_file_ch
    
    output:
    file ("${params.accession}.manifest.txt") into data_prepped_ch

    shell:
    fastqs = fastq.join(' ')
    """
    ls $fastqs > "${params.accession}.list"
    data_prep.pl --list "${params.accession}.list" --sample $params.accession > "${params.accession}.manifest.txt"
    """
}

// ------------------------------------- //

// Parse source manifest file.
data_prepped_ch
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.SAMPLE,row.FILE,row.PAIR_FILE,row.SIZE,row.PAIR_SIZE,row.FORMAT,row.COMPRESSION,row.PAIRING,row.PLATFORM,row.QFORMAT,row.RESTORE,row.RG_ID,row.LIBRARY,row.PU) }
    .into { multifq2bam_in_ch; multibam_in_ch }

// ------------------------------------- //

process multifastq2bam {
    tag { "${params.accession}" }
    label 'picard'

    input:
    val record from multifq2bam_in_ch

    output:
    set val(record), file("*.bam") into lossless_check_ch
    
    shell:
    """ 
    java -jar /opt/picard/picard.jar FastqToSam \\
        TMP_DIR=./ \\
        FASTQ=${record[1]} \\
        FASTQ2=${record[2]} \\
        OUTPUT="${params.accession}-${record[11]}.bam" \\
        SAMPLE_NAME=${record[0]} \\
        PLATFORM=${record[8]} \\
        READ_GROUP_NAME=${record[11]} \\
        LIBRARY_NAME=${record[12]} \\
        PLATFORM_UNIT=${record[13]} \\
    """
}

// ------------------------------------ //

process ubam_validate {
    tag { "${params.accession}" }
    label 'ucgdmods'

    input:
    set val(record), file(bam) from lossless_check_ch

    shell:
    """
    lossless_validator2.pl \\
    $bam \\
    ${record[1]} \\
    ${record[2]} \\
    -cpus \$SLURM_CPUS_ON_NODE \\
    -qual -id \\
    > "${params.accession}.lossless.result"

    cat "${params.accession}.lossless.result" | perl -lane 'if (\$_ =~ /FAILURE/) { die "PolishedUBAM not lossless."; exit(1) }'
    """
}

// ------------------------------------ //

// Need to collect and merge bams.
fq2bam_out_ch
    .groupTuple(by:0)
    .set { bam_tuple_ch }

// ------------------------------------ //

process mergeBams {
    tag { "${params.accession}" }
    label 'samtools'
    publishDir ${processing}, mode: 'link', pattern: "*.bam*"

    input:
    set val(sample_id), val(collect) from bam_tuple_ch

    output:
    file("${params.accession}.bam") into merge_complete_ch

    shell:
    bams = collect.join(' ')
    """
    samtools merge \\
    -n -c -p \\
    -@ 10 \\
    "${params.accession}.bam" \\
    $bams
    """
}    

// ------------------------------------ //

process bam2fastq {
    label 'ucgdmods'
    tag { "${params.accession}" }

    input:
    file(ubam) into merge_complete_ch

    output:
    set file("${params.accession}.align.txt") into bam2fastq_ch

    shell:
    '''
    bam2fastq2.pl !{ubam} -fix -id -z -align_file "!{params.accession}.align.txt" -c $SLURM_CPUS_ON_NODE
    '''
}

// ------------------------------------ //

// Split for fastq2bam and fastp
bam2fastq_ch.into{ fastq2bam_in_ch; parse_in_ch }

// ------------------------------------- //

process fastq2bam {
    label 'fqf'
    tag { "${params.accession}" }
 
    input:
    file(alignFile) from fastq2bam_in_ch
    
    output:
    set file("${params.accession}.bam"), file("${params.accession}.bam.bai") into polishedBam_ch

    shell:
    '''
    envclean.sh

    ibrun perl -S FastQforward.pl \\
    fastq2bam \\
    -align_file !{alignFile} \\
    -ref !{params.reference} \\
    -realign_indels \\
    -known_indels !{params.indel_mills} \\
    -sentieon \\
    -hyperthread \\
    -outfile "!{params.accession}.bam"

    envclean.sh
    '''
}

// ------------------------------------- //

process parseAlignFile {
    label 'local'
    tag { "${params.accession}" }

    input:
    file(alignFile) from parse_in_ch

    output:
    file("*.rg.txt") into align_out_ch

    shell:
    '''
    UCGDalign -af !{alignFile} 
    '''
}

// ------------------------------------- //

process fastp {
    label 'fastp'
    tag { "${params.accession}" }
 
    input:
    file(alignFiles) from align_out_ch.transpose()

    output:
    val 'ready' into fastp_out_ch

    shell:
    '''
    export FASTQLINE=$(cat !{alignFiles} | perl -F, -lane '$F[1] =~ s|\\s+| --in2=|; print $F[1]')

    fastp \\
    --thread 10 \\
    --in1=$FASTQLINE \\
    --json "!{params.accession}.!{alignFiles.baseName}.fastp.json"
    ''' 
}

// ------------------------------------- //

// split channel to create GVCFs and lossless check.
polishedBam_ch.into { makeGVCF_ch; run_alignstats_ch; validate_bam_ch; makeCRAM_ch }

// ------------------------------------- //

process polish_validate {
    label 'ucgdmods'
    tag { "${params.accession}" }

    input:
    set file(pbam), file(index) from validate_bam_ch

    output:
    val 'ready' into losslessVal_done_ch

    shell:
    '''
    lossless_validator2.pl \\
    "!{processing}/!{pbam}" \\
    !{pbam} \\
    -cpus $SLURM_CPUS_ON_NODE \\
    -qual -id \\
    > "!{params.accession}.lossless.result"

    cat "!{params.accession}.lossless.result" | perl -lane 'if ($_ =~ /FAILURE/) { die "PolishedBam not lossless."; exit(1) }'
    '''
}

// ------------------------------------- //

process bam2gvcf {
    label 'fqf'
    tag { "${params.accession}" }

    input:
    set file(polishedbam), file(bamIndex) from makeGVCF_ch
    
    output:
    file("${params.accession}.g.vcf.gz") into bam2gvcf_ch
    file("${params.accession}.g.vcf.gz.tbi") into bam2gvcf_index_ch

    script:
    if ( "$workflow.profile" == 'standard') 
    """
    envclean.sh

    ibrun perl -S FastQforward.pl \\
    bam2gvcf \\
    -ref ${params.reference} \\
    -i ${polishedbam} \\
    -recalibrate \\
    -known_snps ${params.dbsnp} \\
    -sentieon \\
    -hyperthread \\
    ${params.pcrLibrary} \\
    -annotation ${params.fqfAnnos} \\
    -outfile "${params.accession}.g.vcf.gz"

    ## index the new gvcf.    
    sentieon util vcfindex "${params.accession}.g.vcf.gz"

    envclean.sh
    """

    else if ( "$workflow.profile" == 'grch37') 
    """
    envclean.sh

    ibrun perl -S FastQforward.pl \\
    bam2gvcf \\
    -ref ${params.reference} \\
    -i ${polishedbam} \\
    -recalibrate \\
    -known_snps ${params.dbsnp} \\
    -sentieon \\
    -hyperthread \\
    -outfile "${params.accession}.g.vcf.gz"

    ## index the new gvcf.    
    sentieon util vcfindex "${params.accession}.g.vcf.gz"

    envclean.sh
    """
}

// ------------------------------------- //

process samtoolsCRAMer {
    label 'samtools'
    tag { "${params.accession}" }

    input:
    set file(pbam), file(bamIndex) from makeCRAM_ch

    shell:
    '''
    samtools view \\
    -T !{params.reference} \\
    !{pbam} \\
    -o "!{params.accession}.cram" \\
    -C \\
    -@ 10

    ## index the new cram.
    samtools index "!{params.accession}.cram" 
    '''
}

// ------------------------------------- //

process alignstats {
    label 'alignstats'
    tag { "${params.accession}" }

    input:
    set file(pbam), file(index) from run_alignstats_ch
   
    output:
    file("*alignstats.json")
    val 'ready' into alignstats_done_ch 
 
    shell:
    '''
    alignstats -P $SLURM_CPUS_ON_NODE -i !{pbam} -j bam -o "!{params.accession}.alignstats.json" -t !{params.regionBeds}
    '''
}

// ------------------------------------- //

