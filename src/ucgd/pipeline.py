#!/usr/bin/env python3

import os
import sys
import re
import subprocess
from multiprocessing import Process, Manager

class ProjectError(Exception):
    """
    Raised when processing space does not exist.
    """
    pass

class PipelineUCGD:

    def __init__(self, args):

        ## put hardcoded path in object.
        self.nicu_drop = args.nicu_drop
        self.processing_space = args.irb_path
        self.rclone_config = args.rclone_config

        ## ubox manifest drop location.
        self.ubox_drop = "ubox:UCGD_Team/NeoSeq/NeoSeq_manifest_final/"
        self.nicu_project_ubox = "ubox:UCGD_Team/NeoSeq/NeoSeq_projects/"

        nf_scripts = {
                'var' : 'ucgd.nicu.nf',
                'sv' : 'ucgd.sv.nf',
                'process_vcf' : 'ucgd.vcf.preprocess.nf',
                'mosaic' : 'ucgd.mosaic.nf',
                'config' : 'ucgd.master.config',
        }
        self.nf_scripts = nf_scripts

    ## ---------------------------------------- ##
    
    def __get_project_root(self, project, nicu):
        path = self.processing_space + project

        try:
            if not os.path.exists(path):
                raise ProjectError('Project Path {} does not exist'.format(path))
        except ProjectError as e:
            nicu.log.error(e)
                
        for root, dirs, files in os.walk(path):
            for f in files:
                if re.match(r'^ucgd.*nf', f):
                    return root

    ## ---------------------------------------- ##
   
    def reset_acl(self, nicu):
        acl_cmd = 'reset-acl -R {}'.format(self.processing_space)
        return_state  = self.__cmd_worker(acl_cmd)

        if return_state.returncode == 0:
            nicu.log.info('reset-acl could be ran on directory: {}'.format(self.processing_space))
        else:
            nicu.log.info('reset-acl ran on directory: {}'.format(self.processing_space))
    
    ## ---------------------------------------- ##
    
    def md5_check(self, accession, project, nicu):
        """
        Will collect md5 files from drop location and run
        md5sum --check on files. Return dict of [file:exitcode]
        """
        acc_dir = self.nicu_drop + accession
        os.chdir(acc_dir)
   
        md5_files = []
        for root, dirs, files in os.walk(os.getcwd()):
            for i in files:
                if 'md5' in i:
                    md5_files.append(i)

        if not md5_files:
            nicu.sns.missing_checksum(accession)
            nicu.log.error('No accession {} md5sum files found for project {}'.format(accession, project))

        ## multiprocessing manager.
        manager = Manager()
        return_dict = manager.dict()

        jobs = []
        for md5 in md5_files:
            p = Process(target=self.__md5_worker, args=(md5, return_dict))
            jobs.append(p)
            p.start()
      
        ## wait for finish.
        for proc in jobs:
            proc.join()

        return return_dict, md5_files

    ## ---------------------------------------- ##
    
    def __md5_worker(self, md5, return_dict):
        check_run = subprocess.run(["md5sum", "--check", md5], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return_dict[md5] = check_run.returncode

    ## ---------------------------------------- ##

    def launch_projects(self, project, args, nicu, mosaic_id):
        """
        Will launch each step of the NeoSeq pipelie base on request.
        """
        ## get base directory of project.
        project_root = self.__get_project_root(project, nicu)
        ## check project ready
        try:
            if project_root == None:
                raise ProjectError
        except ProjectError as e:
                nicu.sns.project_issue('Error occured trying to validate processing path for project {}'.format(project))
                nicu.log.error("Error occured trying to validate processing path for project {}".format(project))

        ## build the remaining location paths.
        project_setup = self.processing_space + project + '/Project_Setup'
        vcf_dir = project_root + '/VCF/Complete/'
        ped_dir = project_root + '/Reports/'
        cram_dir = project_root + '/Data/PolishedCrams/'
        fabric_dir = project_root + '/Analysis/Fabric_VCFs/'
        slivar_dir = project_root + '/Analysis/Slivar/'
        irb_reset = self.processing_space + project
        log_file  = '{}/{}_{}.log'.format(project_root, project, args.pipeline_step)
        trace_log = '{}/UCGD_Trace_Report_{}.txt'.format(project_root, args.pipeline_step)

        ## if running test.
        test_run = 'notest'
        if args.test:
            test_run = 'test'

        ## VAR command.
        var_template = 'nextflow -log {} run {}/{} -c {}/{} --fastqs {} --project {} -work-dir {}/work -with-trace {} --project_path {} --irb_reset {} --neotest {} -resume'
        var_command = var_template.format(log_file, project_root, self.nf_scripts['var'], project_root, self.nf_scripts['config'], project_setup, project, project_root, trace_log, project_root, irb_reset, test_run)

        ## SV command.
        sv_template = 'nextflow -log {} run {}/{} -c {}/{} --project {} -work-dir {}/work -with-trace {} --project_path {} --crams {} --neotest {} -resume'
        sv_command = sv_template.format(log_file, project_root, self.nf_scripts['sv'], project_root, self.nf_scripts['config'], project, project_root, trace_log, project_root, cram_dir, test_run)

        ## Preprocess vcf to fabric upload.
        post_template = 'nextflow -log {} run {}/{} -c {}/{} --project {} -work-dir {}/work -with-trace {} --project_path {} --vcf {} --ped {} --neotest {} -resume'
        post_command = post_template.format(log_file, project_root, self.nf_scripts['process_vcf'], project_root, self.nf_scripts['config'], project, project_root, trace_log, project_root, vcf_dir, ped_dir, test_run)

        ## mosaic json and upload steps.
        mosaic_template = 'nextflow -log {} run {}/{} -c {}/{} --project {} -work-dir {}/work -with-trace {} --project_path {} --vcfs {} --crams {} --slivar {} --mosaic_project_id {} --neotest {} -resume'
        mosaic_command = mosaic_template.format(log_file, project_root, self.nf_scripts['mosaic'], project_root, self.nf_scripts['config'], project, project_root, trace_log, project_root, fabric_dir, cram_dir, slivar_dir, mosaic_id, test_run)

        commands = {
                'VAR' : var_command,
                'SV'  : sv_command,
                'POST' : post_command,
                'MOSAIC'  : mosaic_command,
                'ALL' : [var_command, post_command, sv_command, mosaic_command],
        }

        ## Run all command or as requested.
        if args.pipeline_step == 'ALL':
            nicu.ucgd_db.update_project_status(project, 'processing')
            nicu.log.info('Processing for project {} started.'.format(project))
            nicu.sns.processing_started(project)
            nicu.sns.time_tracker(project, 'NeoSeq pipeline ALL processing started')
            for cmd in commands['ALL']:
                nicu.log.info('Processing command: {}'.format(cmd))
                return_state  = self.__cmd_worker(cmd)

                if return_state.returncode == 0:
                    ## send sns for each step.
                    for key, value in commands.items():
                        if value == cmd:
                            nicu.sns.data_alert(project, key) 
                            nicu.sns.time_tracker(project, key)
                else:
                    return return_state, self.processing_space
        else:
            nicu.ucgd_db.update_project_status(project, 'processing')
            nicu.log.info('Processing for project {} started. With command: {}'.format(project, commands[args.pipeline_step]))
            nicu.sns.processing_started(project)
            nicu.sns.time_tracker(project, args.pipeline_step)
            return_state  = self.__cmd_worker(commands[args.pipeline_step])
            return return_state, self.processing_space

    ## ---------------------------------------- ##

    def __cmd_worker(self, cmd):
        run = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return run

    ## ---------------------------------------- ##
    
    def get_ubox_projects(self, nicu):
        """
        Collect all ubox nicu projects and returns one per request.
        """
        xlsx_check = subprocess.run(["rclone", "--config", self.rclone_config, "ls", self.ubox_drop], stdout=subprocess.PIPE, check=True)
        cmd_stdout = xlsx_check.stdout
        files = cmd_stdout.decode().split()

        return_xlsx = []
        if not files:
            return False
        else:
            for manifest in files:
                if not manifest.endswith('xlsx'):
                    continue
                dest_path = '.'
                drop_file = self.ubox_drop + manifest
                xlsx_get = subprocess.run(["rclone", "--config", self.rclone_config, "copy", drop_file, dest_path], stdout=subprocess.PIPE, check=True)
                return_xlsx.append(manifest) 
                if xlsx_get.returncode != 0:
                    nicu.log.error('Could not create UCGD NeoSeq ubox folder.')
                    return False
                return_xlsx.append(manifest) 

        return return_xlsx.pop()
    
    ## ---------------------------------------- ##
    
    def build_ubox_project(self, project_name, manifest, nicu):
        """
        1. Move xlsx files from ubox 'final' location.
        2. Build ubox project directory 
        3. Move xlsx file to new project directory.
        4. Delete local copy.
        """
        new_project_path = self.nicu_project_ubox + project_name + "/" + manifest
        project_build = subprocess.run(["rclone", "--config", self.rclone_config, "moveto", manifest, new_project_path], stdout=subprocess.PIPE, check=True)

        if project_build.returncode != 0:
            nicu.log.error('Issue creating ubox project.')
            os.exit(1)
        else:
            old_xlsx = self.ubox_drop + manifest
            delete_ubox_xlsx = subprocess.run(["rclone", "--config", self.rclone_config, "deletefile", old_xlsx], stdout=subprocess.PIPE, check=True)
            if project_build.returncode != 0:
                nicu.sns.project_issue('NeoSeq manifest likely has file lock. Log into ubox and manually delete as soon as possible.')
                nicu.log.info('NeoSeq manifest likely has file lock. Log into ubox and manually delete as soon as possible.')

    ## ---------------------------------------- ##

