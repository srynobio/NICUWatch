#!/usr/bin/env nextflow

/*
 * UCGD NeoSeq Pipeline
 * Author: Shawn Rynearson
 * Copyright (c) 2020, Utah Center for Genetic Discovery.
*/

// hard coded manifest used.
manifest_input = Channel.fromPath("manifest/NeoSeq_Test_Manifest.xlsx").ifEmpty { error "Template NeoSeq_Test_Manifest.xlsx not found." }

// Set needed irb paths for testing 
test_irb_path  = '/scratch/ucgd/lustre/UCGD_Datahub/IRBs/IRB_00000000/'
test_nicu_drop = '/scratch/ucgd/lustre/UCGD_Datahub/IRBs/IRB_00000000/Staging'

// ------------------------------------- //

process upload_test_manifest {
    executor 'local'
    module 'ucgd_modules/dev'

    input:
    file manifest from manifest_input

    output:
    val 'copied' into upload_done_ch

    shell:
    '''
    rclone \\
    --config /uufs/chpc.utah.edu/common/HIPAA/ucgd-pepipeline/.config/rclone/rclone.conf \\
    copyto \\
    !{manifest} \\
    ubox:UCGD_Team/NeoSeq/NeoSeq_manifest_final/!{manifest} \\
    -L
    '''
}

// ------------------------------------- //

process manifest_check {
    executor 'local'
    module 'ucgd_modules/dev'

    input:
    val 'copied' from upload_done_ch

    output:
    val 'ready' into manifest_checked_ch

    shell:
    '''
    NICUWatch manifest_check --irb_path !{test_irb_path} --nicu_drop !{test_nicu_drop} --test
    '''
}

// ------------------------------------- //

process drop_check {
    executor 'local'
    module 'ucgd_modules/dev'

    input:
    val 'ready' from manifest_checked_ch

    output:
    val 'checked' into drop_checked_ch 

    shell:
    '''
    NICUWatch drop_check --irb_path !{test_irb_path} --nicu_drop !{test_nicu_drop} --test 
    '''
}

// ------------------------------------- //

process process_check {
    executor 'local'
    module 'ucgd_modules/dev'

    input:
    val 'checked' from drop_checked_ch 

    output:
    val 'checked' into process_checked_ch 

    shell:
    '''
    NICUWatch process_check --irb_path !{test_irb_path} --nicu_drop !{test_nicu_drop} --test 
    ''' 
}

// ------------------------------------- //

process run_pipeline_all {
    executor 'local'
    module 'ucgd_modules/dev'

    input:
    val 'checked' from process_checked_ch 

    shell:
    '''
    NICUWatch run_pipeline --pipeline_step ALL --irb_path !{test_irb_path} --nicu_drop !{test_nicu_drop} --test
    '''
}

// ------------------------------------- //

println('Delete test information from:')
println('ubox: https://uofu.app.box.com/folder/79662845561')
println('Mosaic: https://mosaic.chpc.utah.edu/#/projects')
println('Fabric: https://app.fabricgenomics.com/w/2286/p/#/projects')
println('UCGD Database sample and project name: T101-101-NEO-Cyberdyne')
