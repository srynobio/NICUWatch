#!/usr/bin/env nextflow

/*
 * UCGD NeoSeq Pipeline
 * Author: Shawn Rynearson
 * Copyright (c) 2020, Utah Center for Genetic Discovery.
*/
 
// Log info.
log.info "================================="
log.info "       UCGD NeoSeq Pipeline        "
log.info "================================="
log.info "Version: $params.version"
log.info "================================="

// --fastqs path to processing
if( !params.fastqs) { exit 1, "Path to Fastqs required." }
fastqs = Channel.fromPath("$params.fastqs/*.fastq.gz")
    .set{ check_bgzf_in_ch }  

// --project_path 
if( !params.project_path) { exit 1, "Path to processing directory required." }

// --project
if( !params.project) { exit 1, "Project name required." }

// Hard coded values for nicu processing.
params.regionBeds    = '/scratch/ucgd/lustre/common/data/Regions/GRCh38/UCGD.WGS.Region.bed'
params.backgrounds   = '/scratch/ucgd/lustre/common/data/1000G/GRCh38'
params.fqfAnnos      = 'all'
params.pcrLibrary    = 'nopcr'

// ------------------------------------- //

process reset_acl {
    label 'local'
    tag {"$params.project"}
    errorStrategy = 'ignore'

    output:
    val 'reset' into reset_done_ch

    shell:
    '''
    reset-acl -R !{params.irb_reset} 
    '''
}

// ------------------------------------- //

process fastq_bgzf_eof_check {
    label 'localterm'
    tag {"$params.project"}

    input:
    val 'reset' from reset_done_ch
    each fastq from check_bgzf_in_ch

    output:
    val 'complete' into fastq_checked_ch
 
    shell:
    '''
    fastq_bgzf_eof.pl !{fastq}
    check_bgzf.pl !{fastq} | perl -lane 'if (\$_ !~ /VALIDATION SUCCESS/) { die "BGZF format and expected size error."; exit(1) }'
    '''
}

// ------------------------------------- //

Channel.from fastq_checked_ch
    .toList()
    .set { checked_fastqs_ch }

// ------------------------------------- //

process run_data_prep {
    label 'localterm'
    tag {"$params.project"}

    input:
    val 'complete' from checked_fastqs_ch

    output:
    file("${params.project}.align.txt") into prepped_out_ch

    script:
    if (params.neotest == 'test')
    """
    NICUWatch create_prep_file -p ${params.project} --test
    arup_prep.pl --list prep_file.txt > ${params.project}.align.txt
    """
    else if (params.neotest == 'notest')
    """
    NICUWatch create_prep_file -p ${params.project}
    arup_prep.pl --list prep_file.txt > ${params.project}.align.txt
    """
}

// ------------------------------------- //

prepped_out_ch
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.SAMPLE,row.FILE,row.PAIR_FILE,row.SIZE,row.PAIR_SIZE,row.FORMAT,row.COMPRESSION,row.PAIRING,row.PLATFORM,row.QFORMAT,row.RESTORE,row.RG_ID,row.LIBRARY,row.PU) }
    .into { fastp_in_ch; fq2bam_in_ch }

// ------------------------------------- //

process fastp {
    label 'fastp'
    tag { "${sample[0]}" }
    publishDir "${params.project_path}/Reports/fastp", mode: 'link', pattern: "*.fastp.json"

    input:
    val sample from fastp_in_ch

    output:
    val 'complete' into fastp_done_ch
    file("*.fastp.json") into fastp_done

    shell:
    '''     
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    fastp \\
    --thread 10 \\
    --in1=!{sample[1]} \\
    --in2=!{sample[2]} \\
    --json !{sample[0]}-!{sample[11]}.fastp.json
    '''
}

// ------------------------------------- //
 
process fastq2bam {
    label 'fqf'
    tag { "${sample[0]}" }
 
    input:
    val sample from fq2bam_in_ch 
    
    output:
    set val("${sample[0]}"), val("${sample[1]}"), val("${sample[2]}"), file("*.bam"), file("*.bam.bai") into polishedBam_ch

    shell:
    '''
    export ALIGN="Files=!{sample[1]},!{sample[2]};Type=!{sample[7]};RG=@RG\\tID:!{sample[11]}\\tSM:!{sample[0]}\\tPL:!{sample[8]}\\tPU:!{sample[13]}\\tLB:!{sample[12]}"
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    envclean.sh

    ibrun perl -S FastQforward.pl \\
    fastq2bam \\
    -align $ALIGN \\
    -ref !{params.reference} \\
    -realign_indels \\
    -known_indels !{params.indel_mills} \\
    -sentieon \\
    -hyperthread \\
    -outfile !{sample[0]}.bam

    envclean.sh
    '''
}

// ------------------------------------- //

// split channel to create GVCFs, lossless check & QC.
polishedBam_ch.into { makeGVCF_ch; run_alignstats_ch; validate_bam_ch; makeCRAM_ch }

// ------------------------------------- //

process bam2gvcf {
    label 'fqf'
    tag { $sample }
    publishDir "${params.project_path}/VCF/GVCFs", mode: 'link', pattern: "*.g.vcf.gz"
    publishDir "${params.project_path}/VCF/GVCFs", mode: 'link', pattern: "*.g.vcf.gz.tbi"

    input:
    set val(sample), val(fq1), val(fq2), file(pbam), file(index) from makeGVCF_ch 

    output:
    file("${sample}.g.vcf.gz") into bam2gvcf_ch
    file("${sample}.g.vcf.gz.tbi") into bam2gvcf_index_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    envclean.sh

    ibrun perl -S FastQforward.pl \\
    bam2gvcf \\
    -ref !{params.reference} \\
    -i !{pbam} \\
    -recalibrate \\
    -known_snps !{params.dbsnp} \\
    -sentieon \\
    -hyperthread \\
    !{params.pcrLibrary} \\
    -annotation !{params.fqfAnnos} \\
    -outfile "!{sample}.g.vcf.gz"

    ## index the new gvcf.    
    sentieon util vcfindex "!{sample}.g.vcf.gz"
    envclean.sh
    '''
}

// ------------------------------------- //

process losslessValidate {
    label 'ucgdmods'
    errorStrategy = 'terminate'
    tag { $sample }

    input:
    set val(sample), val(fq1), val(fq2), file(pbam), file(index) from validate_bam_ch

    shell:
    '''
    lossless_validator2.pl \\
    !{pbam} \\
    !{fq1} \\
    !{fq2} \\
    -cpus $SLURM_CPUS_ON_NODE \\
    -qual -id \\
    > "!{sample}.lossless.result"

    cat "!{sample}.lossless.result" | perl -lane 'if ($_ =~ /FAILURE/) { die "PolishedBam not lossless."; exit(1) }'
    '''
} 

// -----------------------------------

process samtoolsCRAMer {
    label 'samtools'
    tag { $sample }
    publishDir "${params.project_path}/Data/PolishedCrams", mode: 'link', pattern: "*.cram"
    publishDir "${params.project_path}/Data/PolishedCrams", mode: 'link', pattern: "*.crai"

    input:
    set val(sample), val(fq1), val(fq2), file(pbam), file(index) from makeCRAM_ch

    output:
    val 'ready' into run_goleft_ch
    set val("${sample}"), file("${sample}.cram"), file("${sample}.cram.crai") into cram_out_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    samtools view \\
    -T !{params.reference} \\
    !{pbam} \\
    -o "!{sample}.cram" \\
    -C \\
    -@ $SLURM_CPUS_ON_NODE 

    ## index the new cram.
    samtools index "!{sample}.cram" 
    '''
}

// ------------------------------------- //

process alignstats {
    label 'alignstats'
    tag { "${params.accession}" }
    publishDir "${params.project_path}/Reports/alignstats", mode: 'link', pattern: "*.json"

    input:
    set val(sample), val(fq1), val(fq2), file(pbam), file(index) from run_alignstats_ch
   
    output:
    file("*alignstats.json") into align_done
 
    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    alignstats -P $SLURM_CPUS_ON_NODE -i !{pbam} -j bam -o "!{sample}.alignstats.json" -t !{params.regionBeds}
    '''
}

// ------------------------------------- //

// collect all GVCFs
Channel.from bam2gvcf_ch
    .toSortedList()
    .set { gvcfCollect_ch }

// collect all GVCFs
Channel.from bam2gvcf_index_ch
    .toSortedList()
    .set { gvcfCollect_index_ch }

// Collect standard region bed files GVCFtyper.
// Always run on full genome.
Channel.fromPath("$params.wgsBeds/*bed")
    .map { file -> tuple(file.baseName, file) }
    .set { bedRegions_in_ch } 

// Collect backgrounds
Channel.fromPath("$params.backgrounds/*vcf.gz")
    .unique()
    .toSortedList()
    .set { backgroundInput }

// ------------------------------------- //

// split the channel for typer and vep.
bedRegions_in_ch.into { typer_bed_in_ch; vep_bed_in_ch }

// ------------------------------------- //

process gvcfTyper {
    label 'sentieon'
    tag {"$params.project"}

    input:
    val gvcfs from gvcfCollect_ch
    val index from gvcfCollect_index_ch
    val bks from backgroundInput
    set val(bed), file(bedFile) from typer_bed_in_ch

    output:
    file ("${bed}.g.vcf.gz") into gvcf_typed_out_ch  

    script:
    backgrounds = bks.join(' -v ')
    inputGVCFs  = gvcfs.join(' -v ')
    """
    export TMPDIR="/scratch/local/\$USER/\$SLURM_JOB_ID"

    sentieon driver \\
    --thread_count \$SLURM_CPUS_ON_NODE \\
    -r ${params.reference} \\
    --interval ${bedFile} \\
    --algo GVCFtyper \\
    "${bed}.g.vcf.gz" \\
    ${params.gvcf_annos} \\
    -v ${backgrounds} \\
    -v ${inputGVCFs}
    """
}

// ------------------------------------ //

// Channel to collect chr files.
Channel.from gvcf_typed_out_ch  
    .unique()
    .toSortedList()
    .set { chrTyped_in_ch }

// ------------------------------------ //

process mergeGVCFs {
    label 'bcftools'
    tag {"$params.project"}

    input:
    val chrFiles from chrTyped_in_ch

    output:
    set file("${params.project}_merged.vcf.gz"), file("${params.project}_merged.vcf.gz.tbi") into merged_ch

    shell:
    inputVCFs = chrFiles.join(' ')
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    bcftools concat --threads 10 -O z !{inputVCFs} -o "!{params.project}_merged.vcf.gz"
    tabix -p vcf "!{params.project}_merged.vcf.gz"
    '''
}

// ------------------------------------ //

process varCalSnp {
    label 'sentieon'
    tag {"$params.project"}

    input:
    set file(merged_vcf), file(index) from merged_ch

    output:
    set file(merged_vcf), file(index), file("${params.project}_recal.tranches.snp"), file("${params.project}_recal.snp.vcf.gz"), file("${params.project}_recal.snp.vcf.gz.tbi") into VarCalSNP_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    sentieon driver \\
    --thread_count \$SLURM_CPUS_ON_NODE \\
    -r !{params.reference} \\
    --interval !{params.regionBeds} \\
    --algo VarCal \\
    --vcf !{merged_vcf} \\
    --var_type SNP \\
    --tranches_file "!{params.project}_recal.tranches.snp" \\
    --resource !{params.hapmap}  \\
    --resource_param !{params.hapmap_par} \\
    --resource !{params.omni} \\
    --resource_param !{params.omni_par} \\
    --resource !{params.snp_1G} \\
    --resource_param !{params.snp_1G_par} \\
    --resource !{params.dbsnp} \\
    --resource_param !{params.dbsnp_par} \\
    !{params.vqsr_annos} \\
    "!{params.project}_recal.snp.vcf.gz"
    '''
}

// ------------------------------------ //

process applyVarCalSNP {
    label 'sentieon'
    tag {"$params.project"}

    input:
    set file(merged_vcf), file(merged_index), file(snpTranchFile), file(recalVCF), file(index) from VarCalSNP_ch

    output:
    set file("${params.project}_applyRecal.snp.vcf.gz"), file("${params.project}_applyRecal.snp.vcf.gz.tbi") into applyedSNP_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    sentieon driver \\
    --thread_count $SLURM_CPUS_ON_NODE \\
    -r !{params.reference} \\
    --algo ApplyVarCal \\
    --vcf !{merged_vcf} \\
    --var_type SNP  \\
    --recal !{recalVCF} \\
    --tranches_file !{snpTranchFile} \\
    --sensitivity 99.9 \\
    "!{params.project}_applyRecal.snp.vcf.gz"
    '''
}

// ------------------------------------ //

process varCalIndel {
    label 'sentieon'
    tag {"$params.project"}

    input:
    set file(snp_recal_file), file(index) from applyedSNP_ch

    output:
    set file(snp_recal_file), file(index), file("${params.project}_recal.tranches.indel"), file("${params.project}_recal.indel.vcf.gz"), file("${params.project}_recal.indel.vcf.gz.tbi") into VarCalINDEL_ch
    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    sentieon driver \\
    --thread_count \$SLURM_CPUS_ON_NODE \\
    --interval !{params.regionBeds} \\
    -r !{params.reference} \\
    --algo VarCal \\
    --vcf !{snp_recal_file} \\
    --var_type INDEL \\
    --tranches_file "!{params.project}_recal.tranches.indel" \\
    --resource !{params.indel_mills} \\
    --resource_param !{params.mills_par} \\
    !{params.vqsr_annos} \\
    "!{params.project}_recal.indel.vcf.gz"
    '''
}

// ------------------------------------ //

process applyVarCalINDEL {
    label 'sentieon'
    tag {"$params.project"}

    input:
    set file(snp_recal_file), file(recal_index), file(indelTranch), file(recalVCF), file(index) from VarCalINDEL_ch

    output:
    set file("${params.project}_complete.vcf.gz"), file("${params.project}_complete.vcf.gz.tbi") into applyedINDEL_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    sentieon driver \\
    --thread_count $SLURM_CPUS_ON_NODE \\
    -r !{params.reference} \\
    --algo ApplyVarCal \\
    --vcf !{snp_recal_file} \\
    --var_type INDEL  \\
    --recal !{recalVCF} \\
    --tranches_file !{indelTranch} \\
    --sensitivity 99.9 \\
    "!{params.project}_complete.vcf.gz"
    '''
}

// ------------------------------------ //

// Split to run stats and epoch the file.
applyedINDEL_ch.into { runVEP_ch; runVCFstats_ch }

// ------------------------------------ //

process generateSampleFile {
    label 'localterm'
    tag {"$params.project"}

    output:
    file('samples.list') into sampleList_out_ch

    script:
    if (params.neotest == 'test')
    """
    NICUWatch samples_list -p ${params.project} --test
    """
    else if (params.neotest == 'notest')
    """
    NICUWatch samples_list -p ${params.project}
    """
}

// ------------------------------------ //

process finalStats {
    label 'bcftools'
    tag {"$params.project"}
    publishDir "${params.project_path}/Reports/stats", mode: 'link', pattern: "*.stats"

    input:
    file(sampleFile) from sampleList_out_ch
    set file(vcf), file(index) from runVCFstats_ch

    output:
    val 'complete' into finalStats_done_ch 
    file("*.stats") into stats_file

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    bcftools stats --threads 10 --samples-file !{sampleFile} !{vcf} > "!{params.project}_FinalVCF.stats"
    '''
}

// ------------------------------------ //

process updateBED {
    label 'local'
	tag {"$params.project"}

	input:
    set val(bed), file(bedFile) from vep_bed_in_ch

	output:
	file("${bed}.ubed") into bed_out_ch

	shell:
	'''
	perl -lane '@x = split//, $_; print "$F[0]:$F[1]-$F[2]"' !{bedFile} > "!{bed}.ubed"
	'''
}

// ------------------------------------ //

// Channel to collect chr files to split vep run.
Channel.from bed_out_ch
    .toList()
    .set { bedVep_ch }

// ------------------------------------ //

/// look at adding: https://github.com/konradjk/loftee

process runVEP {
    label 'vep'
	tag {"$params.project"}

	input:
    set file(vcf), file(tbi) from runVEP_ch
    each file(bed) from bedVep_ch

	output:
	file("${bed.baseName}.vep.vcf.gz") into vep_out_ch

	shell:
	'''
    export RANGE=$(cat !{bed})
    export PERL5LIB=!{params.vep_plugin}
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    tabix -h !{vcf} $RANGE | \\
    vep --cache --dir_cache !{params.vep_cache} \\
    --format vcf \\
    --fork 10 \\
    --output_file "!{bed.baseName}.vep.vcf.gz" \\
    --vcf --offline \\
    --plugin GO \\
    --plugin LoFtool \\
    --plugin SpliceRegion \\
    --assembly !{params.vep_assembly} \\
    --everything
	'''
}

// ------------------------------------ //

// Channel to collect chr files.
Channel.from vep_out_ch
    .unique()
    .toList()
    .set { vep_collect_out }

// ----------------------------------------------------- //

process mergeVEP {
    label 'bcftools'
	tag {"$params.project"}

	input:
	file(vcfs) from vep_collect_out.transpose()

	output:
	set file("*.vcf.gz"), file("*.vcf.gz.tbi") into merged_VEP_ch

    script:
    """
    export TMPDIR="/scratch/local/\$USER/\$SLURM_JOB_ID"

    # make file list.
    perl -le 'for(1..22, X, Y, M) { print "chr\$_.vep.vcf.gz" }' > files.txt

    ## epoch to mark final file.
    export EPOCH=\$(date +%s)

    bcftools concat \\
    --threads 10 \\
    --file-list files.txt \\
    -O z \\
    --output "${params.project}.vep.vcf.gz"

    tabix -p vcf "${params.project}.vep.vcf.gz"
    """
}

// ------------------------------------ //

process finalVCF {
    label 'local'
    tag {"$params.project"}
    publishDir "${params.project_path}/VCF/Complete/", mode: 'link', pattern: "*.gz"
    publishDir "${params.project_path}/VCF/Complete/", mode: 'link', pattern: "*.tbi"

    input:
    set file(vcf), file(index) from merged_VEP_ch

    output:
    set file("*.vcf.gz"), file("*.vcf.gz.tbi") into finalVCF_out_ch

    shell:
    '''
    export EPOCH=$(date +%s)

    ls "!{params.project}.vep.vcf.gz" | \\
    perl -lane '$orig = $_; $_ =~ s|vep.vcf.gz|Final_$ENV{EPOCH}.vcf.gz|; system "cp $orig $_"'

    ls "!{params.project}.vep.vcf.gz.tbi" | \\
    perl -lane '$orig = $_; $_ =~ s|vep.vcf.gz.tbi|Final_$ENV{EPOCH}.vcf.gz.tbi|; system "cp $orig $_"'
    '''
}

// ------------------------------------ //

process makePedFile {
    label 'local'
    tag {"$params.project"}
    publishDir "${params.project_path}/Reports/", mode: 'link', pattern: "*.ped"

    output:
    file("${params.project}.ped") into ped_out_ch

    script:
    if (params.neotest == 'test')
    """
    NICUWatch make_ped_file -p ${params.project} --test
    """
    else if (params.neotest == 'notest')
    """
    NICUWatch make_ped_file -p ${params.project}
    """
}

// ------------------------------------ //

process peddy {
    label 'peddy'
    tag {"$params.project"}
    publishDir "${params.project_path}/Reports/peddy/", mode: 'link', pattern: "*.csv"
    publishDir "${params.project_path}/Reports/peddy/", mode: 'link', pattern: "*.json"
    publishDir "${params.project_path}/Reports/peddy/", mode: 'link', pattern: "*.png"
    publishDir "${params.project_path}/Reports/peddy/", mode: 'link', pattern: "*.html"
    publishDir "${params.project_path}/Reports/peddy/", mode: 'link', pattern: "*.ped"

    input:
    set file(vcf), file(tbi) from finalVCF_out_ch
    file("${params.project}.ped") from ped_out_ch

    output:
    val 'complete' into peddy_done_ch 
    set file("*csv"), file("*json"), file("*png"), file("*html"), file("*ped") into peddy_done

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    peddy --plot -p 10 --sites hg38 !{vcf} "!{params.project}.ped"
    '''
}

// ------------------------------------ //

// Collect all crams before runing goleft.
run_goleft_ch
    .collect()
    .set { goleft_ready_ch }

// ------------------------------------ //

process goleftIndexCov {
    label 'goleftIndexCov'
    tag {"$params.project"}

    input:
    val 'ready' from goleft_ready_ch

    output:
    val 'complete' into goleft_done_ch
    
    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    goleft indexcov --directory !{params.goleft} -fai !{params.referencefai} --sex "chrX,chrY" !{params.polishedData}/*crai
    '''
}

// ------------------------------------ //

process multiqc {
    label 'multiqc'
    tag {"$params.project"}
    publishDir "${params.project_path}/Reports/", mode: 'link', pattern: "*.html"

    input:
    val 'complete' from fastp_done_ch
    val 'complete' from peddy_done_ch
    val 'complete' from goleft_done_ch
    val 'complete' from finalStats_done_ch

    output:
    file("${params.project}.multiqc.report.html") into multiqc_done

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    multiqc !{params.reports} --force --no-data-dir --filename "!{params.project}.multiqc.report"
    '''
}

// ------------------------------------ //

process rerun_reset_acl {
    tag {"$params.project"}
    label 'local'
    errorStrategy = 'ignore'

    shell:
    '''
    reset-acl -R !{params.irb_reset}
    '''
}

// ------------------------------------ //
