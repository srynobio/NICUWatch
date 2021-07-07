#!/usr/bin/env nextflow

/*
 * UCGD NeoSeq Pipeline
 * Author: Shawn Rynearson
 * Copyright (c) 2020, UStar Center for Genetic Discovery.
*/
 
// Log info.
log.info "================================="
log.info "       UCGD NeoSeq Post Process    "
log.info "================================="
log.info "Version: $params.version"
log.info "================================="

// --vcf $path to processing
if(!params.vcf) { exit 1, "Path to Final VCF file required." }
Channel.fromPath("$params.vcf/*.vcf.gz")
    .set{ vcf_ch }

// --ped $path to processing
if(!params.ped) { exit 1, "Path to project ped file required." }
Channel.fromPath("$params.ped/*.ped")
    .into{ ped_to_normalize_ch; ped_to_slivar }

// --project_path 
if(!params.project_path) { exit 1, "Path to processing directory required." }

// Hard coded values for nicuwatch processing.
params.regionBeds    = '/scratch/ucgd/lustre/common/data/Regions/GRCh38/UCGD.WGS.Region.bed'
params.backgrounds   = '/scratch/ucgd/lustre/common/data/1000G/GRCh38'
params.fqfAnnos      = 'all'
params.pcrLibrary    = 'nopcr'

// ------------------------------------- //

process normalize {
    tag {"$params.project"}
    label 'bcftools'
    publishDir "${params.project_path}/Analysis/", mode: 'link', pattern: "*.gz"
    publishDir "${params.project_path}/Analysis/", mode: 'link', pattern: "*.gz.tbi"

    input:
    file(vcf) from vcf_ch 
    file(ped) from ped_to_normalize_ch 
    
    output:
    file("sample_ped_info.txt") into normalized_ped_out_ch
    set file("${params.project}_norm-fam.vcf.gz"), file("${params.project}_norm-fam.vcf.gz.tbi") into norm_fam_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    grep -iv 'Kindred_ID' !{ped} | cut -f 2 | sort > sample_ped_info.txt
    bcftools norm -Ou -c w -m - -f !{params.reference} --threads \$SLURM_CPUS_ON_NODE -w 10000 !{vcf} | \\
    bcftools view -Oz -S sample_ped_info.txt --no-update -a -c 1 - -o "!{params.project}_norm-fam.vcf.gz"
    tabix -p vcf "!{params.project}_norm-fam.vcf.gz"
    '''
}

// ------------------------------------- //

process slivar_af_filter {
    label 'slivar'
    tag {"$params.project"}

    input:
    set file(vcf), file(tbi) from norm_fam_ch

    output:
    file("${params.project}_norm-rare-fam.vcf") into slivar_af_norm_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    slivar expr --vcf !{vcf} --pass-only -g !{params.gnomad311} --info "INFO.gnomad311_AF < 0.02" --out-vcf "!{params.project}_norm-rare-fam.vcf"
    '''
}

// ------------------------------------- //

process tabix_slivar {
    label 'bcftools'
    tag {"$params.project"}
    publishDir "${params.project_path}/Analysis/", mode: 'link', pattern: "*.gz"
    publishDir "${params.project_path}/Analysis/", mode: 'link', pattern: "*.gz.tbi"

    input:
    file(vcf) from slivar_af_norm_ch

    output:
    set file("${params.project}_norm-fam-rare.vcf.gz"), file("${params.project}_norm-fam-rare.vcf.gz.tbi") into normalized_out_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    bgzip -c !{vcf} > "!{params.project}_norm-fam-rare.vcf.gz"
    tabix -p vcf "!{params.project}_norm-fam-rare.vcf.gz"
    '''
}

// ------------------------------------- //

// split channel to all post processing and slivar annotations.
normalized_out_ch.into { normalized_data_ch; to_slivar_ch } 

// Split norm ped data.
normalized_ped_out_ch
    .splitText()
    .map{ it -> it.trim() }
    .set { split_text_ch } 

// ------------------------------------- //

process split_vcf {
    // split vcf to upload to fabric
    label 'bcftools'
    tag {"$params.project"}
    publishDir "${params.project_path}/Analysis/Fabric_VCFs/", mode: 'link', pattern: "*.gz"
    publishDir "${params.project_path}/Analysis/Fabric_VCFs/", mode: 'link', pattern: "*.gz.tbi"
    publishDir "${params.project_path}/Analysis/Fabric_VCFs/", mode: 'link', pattern: "*.md5"

    input:
    each sample from split_text_ch
    set file(vcf), file(index) from normalized_data_ch

    output:
    set file("${sample}_rare-allpass.vcf.gz"), file("${sample}_rare-allpass.vcf.gz.tbi"), file("${sample}.fabric.md5") into split_finshed_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    bcftools view -s !{sample} !{vcf} | \\
    bcftools annotate -Oz -x FILTER - -o "!{sample}_rare-allpass.vcf.gz"
    tabix -p vcf "!{sample}_rare-allpass.vcf.gz"

    md5sum "!{sample}_rare-allpass.vcf.gz" | perl -lane 'print "md5:$F[0]"' > "!{sample}.fabric.md5"
    '''
}

// ------------------------------------- //

process get_kindred_id_ped {
    tag {"$params.project"}

    input:
    file(ped) from ped_to_slivar

    output:
    file("*_slivar.ped") into slivar_ped_out_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"
    grep -iv 'Kindred_ID' !{ped} | cut -f 1-6 | sort > "!{params.project}_slivar.ped"
    '''
}

// ------------------------------------------------------- //

// Slivar processing.
// https://github.com/UCGD/NeoSeq/blob/master/Slivar/SOP.md

// split the slivar ped
slivar_ped_out_ch.into { slivar_ped_one; slivar_ped_two; slivar_ped_three; slivar_ped_four }

// ------------------------------------------------------- //

process slivar_denovo_first {
    // 1a from SOP
    label 'slivar'
    tag {"$params.project"}

    input:
    file ped from slivar_ped_one
    set file(vcf), file(index) from to_slivar_ch

    output:
    file("*.slivar-first.vcf") into slivar_denovo_first_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    slivar expr \\
    --vcf !{vcf} \\
    --ped !{ped} \\
    --pass-only \\
    -g !{params.gnomad_v2} \\
    -g !{params.gnomad_2x} \\
    -g !{params.gnomad311} \\
    --info "INFO.impactful && INFO.gnomad_popmax_af_controls < 0.02 && INFO.gnomad_popmax_af < 0.02 && \\
        INFO.gnomad311_AF_controls < 0.02 && INFO.gnomad311_AF < 0.02 && INFO.gnomad311_AF_nfe < 0.02 && INFO.gnomad311_AF_amr < 0.02" \
    --js !{params.slivar_js} \\
    --trio "denovo:denovo(kid, mom, dad) && INFO.gnomad_popmax_af_controls < 0.001 && INFO.gnomad_popmax_af < 0.001 && \\
        INFO.gnomad311_AF_controls < 0.001 && INFO.gnomad311_AF < 0.001 && INFO.gnomad311_AF_nfe < 0.001 && INFO.gnomad311_AF_amr < 0.001 && \\
        (kid.sex != 'male' || variant.CHROM != 'chrX' || variant.POS <= 2781479 || variant.POS >= 155701383)" \\
    --trio "x_denovo:x_denovo(kid, mom, dad) && INFO.gnomad_popmax_af_controls < 0.001 && INFO.gnomad_popmax_af < 0.001 && kid.hom_alt && \\
        INFO.gnomad311_AF_controls < 0.001 && INFO.gnomad311_AF < 0.001 && INFO.gnomad311_AF_nfe < 0.001 && INFO.gnomad311_AF_amr < 0.001 && \\
        kid.sex == 'male' && variant.CHROM == 'chrX' && variant.POS > 2781479 && variant.POS < 155701383" \\
    --trio "hom:hom(kid, mom, dad) && INFO.gnomad_nhomalt_controls < 2 && INFO.gnomad_nhomalt < 3 && \\
        INFO.gnomad311_nhomalt_controls < 2 && INFO.gnomad311_nhomalt < 3 && \\
        (kid.sex != 'male' || variant.CHROM != 'chrX' || variant.POS <= 2781479 || variant.POS >= 155701383)" \\
    --trio "x_hemi:x_hemi(kid, mom, dad) && INFO.gnomad_popmax_af_controls < 0.001 && INFO.gnomad_popmax_af < 0.001 && kid.hom_alt && \\
        INFO.gnomad311_AF_controls < 0.001 && INFO.gnomad311_AF < 0.001 && INFO.gnomad311_AF_nfe < 0.001 && INFO.gnomad311_AF_amr < 0.001 && \\
        INFO.gnomad2b38_AC_male_controls < 2 && INFO.gnomad2b38_AC_male < 3 && INFO.gnomad311_AC_controls_XY < 2 && INFO.gnomad311_AC_XY < 3 && \\
        INFO.gnomad_nhomalt_controls < 2 && INFO.gnomad_nhomalt < 3 && INFO.gnomad311_nhomalt_controls < 2 && INFO.gnomad311_nhomalt < 3 && \\
        kid.sex == 'male' && variant.CHROM == 'chrX' && variant.POS > 2781479 && variant.POS < 155701383" \\
    --trio "comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt_controls < 2 && INFO.gnomad_nhomalt < 3 && \\
        INFO.gnomad311_nhomalt_controls < 2 && INFO.gnomad311_nhomalt < 3" \\
    -o "!{params.project}.slivar-first.vcf"
    '''
}

// ------------------------------------------------------- //

// Split the first denovo file.
slivar_denovo_first_ch.into { second_pass_in; first_filter_in } 

// ------------------------------------------------------- //

process slivar_denovo_second {
    // 1b from SOP
    label 'slivar'
    tag {"$params.project"}

    input:
    file ped from slivar_ped_two
    file vcf from second_pass_in

    output:
    file("*.slivar-denovo-second.vcf") into denovo_second_out

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    slivar expr \\
    --vcf !{vcf} \\
    --ped !{ped} \\
    --pass-only \\
    --js !{params.slivar_js} \\
    --trio "denovo:denovo(kid, mom, dad) && (kid.sex != 'male' || variant.CHROM != 'chrX' || variant.POS <= 2781479 || variant.POS >= 155701383)" \\
    --trio "x_denovo:x_denovo(kid, mom, dad) && kid.sex == 'male' && variant.CHROM == 'chrX' && variant.POS > 2781479 && variant.POS < 155701383" \\
    --trio "hom:hom(kid, mom, dad) && (kid.sex != 'male' || variant.CHROM != 'chrX' || variant.POS <= 2781479 || variant.POS >= 155701383)" \\
    --trio "x_hemi:x_hemi(kid, mom, dad) && kid.sex == 'male' && variant.CHROM == 'chrX' && variant.POS > 2781479 && variant.POS < 155701383" \\
    -o "!{params.project}.slivar-denovo-second.vcf"
    '''
}

// ------------------------------------------------------- //

process slivar_filter_first {
    // 2 from SOP
    label 'slivar'
    tag {"$params.project"}

    input:
    file vcf from first_filter_in
    file ped from slivar_ped_three

    output:
    file('*.slivar-comphets.vcf') into second_pass_out_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    slivar compound-hets \\
    --vcf !{vcf} \\
    --ped !{ped } \\
    --sample-field denovo \\
    --sample-field x_denovo \\
    --sample-field comphet_side \\
    -o "!{params.project}.slivar-comphets.vcf"
    '''
}

// ------------------------------------------------------- //

process slivar_bcftools {
    // 3 from SOP
    tag {"$params.project"}
    label 'bcftools'
    publishDir "${params.project_path}/Analysis/Slivar/", mode: 'link', pattern: "*.slivar-all.vcf.gz*"

    input:
    file(vcf1b) from denovo_second_out
    file (vcf2) from second_pass_out_ch

    output:
    set file("*.slivar-all.vcf.gz"), file("*.slivar-all.vcf.gz.tbi") into bcftools_out_ch

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    ## from 1b 
    bgzip !{vcf1b} 
    bcftools index -t "!{vcf1b}.gz"

    ## from 2
    bgzip !{vcf2} 
    bcftools index -t "!{vcf2}.gz"
    
    bcftools concat -Ou -a -D "!{vcf1b}.gz" "!{vcf2}.gz" | bcftools sort -Oz - -o "!{params.project}.slivar-all.vcf.gz"
    bcftools index -t "!{params.project}.slivar-all.vcf.gz"
    '''
}

// ------------------------------------------------------- //

process slivar_tsv {    
    // 4 from SOP
    label 'slivar'
    tag {"$params.project"}
    publishDir "${params.project_path}/Analysis/Slivar/", mode: 'link', pattern: "*.tsv"

    input:
    file(ped) from slivar_ped_four
    set file(vcf), file(tbi) from bcftools_out_ch

    output:
    file("*.tsv") into slivar_tsv_out

    shell:
    '''
    export TMPDIR="/scratch/local/$USER/$SLURM_JOB_ID"

    slivar tsv \\
    --ped !{ped} \\
    -s denovo \\
    -s x_denovo \\
    -s hom \\
    -s x_hemi \\
    -s slivar_comphet \\
    -i gnomad_popmax_af -i gnomad311_AF -i gnomad_nhomalt -i gnomad311_nhomalt -i gnomad2b38_AC_male -i gnomad311_AC_XY \\
    -c CSQ \\
    -g !{params.slivar_pli} \\
    -g !{params.slivar_clinvar} \\
    !{vcf} \\
    -o "!{params.project}.slivar-all.tsv"
    '''
}

// ------------------------------------------------------- //
// End of Slivar section.
// ------------------------------------------------------- //

process push_to_fabric {
    label 'localterm'
    tag {"$params.project"}

    input:
    set file(vcf), file(index), file(md5) from split_finshed_ch

    output:
    val 'pushed' into fabric_push_ch

    script:
    if (params.neotest == 'test')
    """
    export NAME=\$(ls ${vcf} | cut -d '_' -f1)
    export CHECKSUM=\$(cat ${md5})
    NICUWatch fabric -fa upload_genome -p ${params.project} --sample \$NAME --checksum_file \$CHECKSUM --vcf_file ${vcf} --test
    """
    else if (params.neotest == 'notest')
    """
    export NAME=\$(ls ${vcf} | cut -d '_' -f1)
    export CHECKSUM=\$(cat ${md5})
    NICUWatch fabric -fa upload_genome -p ${params.project} --sample \$NAME --checksum_file \$CHECKSUM --vcf_file ${vcf}
    """
}

// ------------------------------------- //

fabric_push_ch
    .collect()
    .set { all_files_pushed_ch }

// ------------------------------------- //

process create_case {
    label 'localterm'
    tag {"$params.project"}

    input:
    val 'all_pushed' from all_files_pushed_ch 

    script:
    if (params.neotest == 'test')
    """
    NICUWatch fabric -fa create_case --hpo_file ${params.hpo_file} -p ${params.project} --test
    """
    else if (params.neotest == 'notest')
    """
    NICUWatch fabric -fa create_case --hpo_file ${params.hpo_file} -p ${params.project}
    """
}

// ------------------------------------- //
