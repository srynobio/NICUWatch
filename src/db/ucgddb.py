#!/usr/bin/env python3

import os
import sys
import re
import datetime
import shutil
import configparser
import subprocess
import pandas as pd
import sqlalchemy as db
from sqlalchemy.sql import table, column, select, update
from collections import defaultdict

class PedigreeError(Exception):
    pass

class ManifestError(Exception):
    pass

class ProjectError(Exception):
    pass

class ProjectCreationError(Exception):
    pass

class DatabaseUCGD:

    def __init__(self, args):
        ## Open found configure file.
        config = configparser.ConfigParser()
        config.read(args.config)

        # database info
        dbname   = config['ucgddb']['dbname']
        host     = config['ucgddb']['host']
        port     = config['ucgddb']['port']
        username = config['ucgddb']['username']
        password = config['ucgddb']['password']

        url = 'postgresql://{}:{}@{}:{}/{}'
        url = url.format(username, password, host, port, dbname)
        conn = db.create_engine(url, client_encoding='utf8', pool_recycle=3600)
        meta = db.MetaData(bind=conn)
        meta.reflect()

        ## put hardcoded paths in object.
        self.nicu_drop = args.nicu_drop
        self.processing_space = args.irb_path

        ## add path to pipeline.
        apps = os.environ['APPS'] + '/'

        if args.test:
            self.workflows = apps + 'NICUWatch/dev/workflows/'
        else:
            self.workflows = apps + 'NICUWatch/production/workflows/'
    
        ## Required fields and print.
        required_fields = {
            "sample_id",
            "arup_accession",
            "pi_last_name",
            "kindred_id",
            "paternal_id",
            "maternal_id",
            "sex",
            "affection_status",
            "phenotype_description",
            "incidental_consent", 
            "hpo_terms",
        }
        self.required_fields = required_fields

        self.conn = conn
        self.meta = meta

    ## ---------------------------------- ##
 
    def build_processing_space(self, nicu):
        """
        Builds needed IRB processing space directories.
        """
        project_name = self.project_name

        ## get the assembly from database.
        projects = self.meta.tables['projects']
        statement = select([projects.c.assembly]).where(projects.c.project == project_name)
        result = self.conn.execute(statement)
        assembly = [x[0] for x in result]

        path_to = self.processing_space + project_name + '/UCGD/' + assembly.pop() + '/'
        process_space_path = self.processing_space + project_name + '/Project_Setup/'

        ## where new directories are added.
        project_dirs = {
            process_space_path, 
            path_to + "Analysis",
            path_to + "Analysis/Fabric_VCFs",
            path_to + "Analysis/Rufus",
            path_to + "Analysis/Slivar",
            path_to + "Analysis/Manta",
            path_to + "Analysis/Smoove",
            path_to + "Analysis/vIQ", 
            path_to + "Analysis/Phenotype_and_Genes",
            path_to + "Data/PolishedCrams",
            path_to + "Reports/fastp",
            path_to + "Reports/goleft",
            path_to + "Reports/peddy",
            path_to + "Reports/RunLogs",
            path_to + "Reports/alignstats",
            path_to + "Reports/stats",
            path_to + "Reports/mosaic",
            path_to + "Reports/manta",
            path_to + "VCF/GVCFs",
            path_to + "VCF/Complete",
            path_to + "VCF/Manta",
            path_to + "VCF/Smoove",
        }

        ## add dirs to project             
        for path in project_dirs:
            if not os.path.exists(path):
                os.makedirs(path)

        ## add nf scripts to project
        for root, dirs, files in os.walk(self.workflows):
            for nf in files:
                path_file = root + nf
                project_file = path_to + nf
                shutil.copyfile(path_file, project_file)

        ## Collect project hpo terms.
        obo_ids = nicu.ucgd_db.get_project_hpo_terms(project_name)
        discovered, undiscovered = nicu.h_api.get_obo_term(obo_ids)

        ## Create term to id file.
        term_to_id_file  = os.path.join(path_to + "Analysis/Phenotype_and_Genes/obo_term_to_id.txt")
        id_to_gene_file  = os.path.join(path_to + "Analysis/Phenotype_and_Genes/obo_id_to_genes.txt")
        hpo_unfound_file = os.path.join(path_to + "Analysis/Phenotype_and_Genes/undiscovered_obo_terms.txt")

        id_file    = open(term_to_id_file, 'w')
        gene_file  = open(id_to_gene_file, 'w')
        unfound_id = open(hpo_unfound_file, 'w')

        for hpo_term, hpo_id in discovered.items():
            id_file.write('{}\t{}\n'.format(hpo_term, hpo_id))
            gene_list = nicu.h_api.get_genes_from_id(hpo_id)
            gene_file.write('{}\t{}\t{}\n'.format(hpo_term, hpo_id, ','.join(gene_list)))
        id_file.close()
        gene_file.close()

        ## Create unfound term file.
        if undiscovered:
            for keys in undiscovered:
                unfound = 'Undiscovered Term: {}\n'.format(keys)
                unfound_id.write(unfound)
        unfound_id.close()

    ## ---------------------------------- ##
    
    def add_fabric_project_id(self, project, project_id):
        projects = self.meta.tables['projects']
        id_stmt = update(projects).where(projects.c.project == project).values(fabric_id = project_id)
        result = self.conn.execute(id_stmt)
        
        if not result.rowcount != 0:
            nicu.log.info('Could not add fabric_project_id to database {}.'.format(project_name))
    
    ## ---------------------------------- ##

    def add_fabric_genome_id(self, genome_id, sample_name, project_name, nicu):
        samples = self.meta.tables['samples']
        sample_stmt = update(samples).where(samples.c.sample_id == sample_name).values(fabric_genome_id = genome_id )
        result = self.conn.execute(sample_stmt)

        if not result.rowcount != 0:
            nicu.log.info('Could not add fabric_genome_id to project {}.'.format(project_name))
    
    ## ---------------------------------- ##
   
    def add_mosaic_project_id(self, mosaic_project_id, project, nicu):
        projects = self.meta.tables['projects']
        id_stmt = update(projects).where(projects.c.project == project).values(mosaic_id = mosaic_project_id)
        result = self.conn.execute(id_stmt)

        if result.rowcount == 0:
            nicu.log.info('Mosaic project_id: {} could not be updated for project {}'.format(mosaic_project_id, project))
        else:
            nicu.log.info('Mosaic project_id: {} updated for project {}'.format(mosaic_project_id, project))

    ## ---------------------------------- ##
    
    def get_nicu_drop(self):
        return self.nicu_drop
    
    ## ---------------------------------- ##

    def get_project_name(self):
        return self.project_name
    
    ## ---------------------------------- ##
    
    def get_kindred_id(self, project):
        samples = self.meta.tables['samples']
        status_stmt = select([samples.c.kindred_id]).where(samples.c.project == project)
        result = self.conn.execute(status_stmt)

        return result.first()[0]
    
    ## ---------------------------------- ##

    def get_fabric_project_id(self, project):
        projects = self.meta.tables['projects']
        status_stmt = select([projects.c.fabric_id]).where(projects.c.project == project)
        result = self.conn.execute(status_stmt)

        return result.first()[0]

    ## ---------------------------------- ##

    def get_fabric_sample_data(self, project):
        """
        Will return dict structured as fabric report API expects for upload.
        """
        samples = self.meta.tables['samples']
        status_stmt = select([samples.c.sample_id, samples.c.affection_status, samples.c.sex, samples.c.maternal_id, samples.c.paternal_id, samples.c.fabric_genome_id]).where(samples.c.project == project)
        result = self.conn.execute(status_stmt)

        samples = {}
        family_count = 1
        for sample_data in result:
            ## update sample sex.
            if sample_data[2].lower() == 'male':
                sex = 'm'
            elif sample_data[2].lower() == 'female':
                sex = 'f'
            ## build dict based on collected values.
            if (sample_data[1] == 'affected' and sample_data[3] != 0 and sample_data[4] != 0):
                samples.update({"proband" : { 'genome_id' : int(sample_data[5]), 'sex' : sex }})
                continue
            if (sex == 'f'):
                member = 'family_{}'.format(family_count)
                affected_status = True if sample_data[1] == 'affected' else False
                samples.update({ member: { 'genome_id' : int(sample_data[5]), 'sex' : sex, 'affected' : affected_status, 'relationship' : 'mother' }})
                family_count += 1
                continue
            if (sex == 'm'):
                member = 'family_{}'.format(family_count)
                affected_status = True if sample_data[1] == 'affected' else False
                samples.update({ member: { 'genome_id' : int(sample_data[5]), 'sex' : sex, 'affected' : affected_status, 'relationship' : 'father' }})
                family_count += 1
                continue

        return samples

    ## ---------------------------------- ##

    def get_project_hpo_terms(self, project):
        """
        Will return a dict of all project table information.
        """
        samples = self.meta.tables['samples']
        status_stmt= select([samples.c.hpo_terms]).where(samples.c.project == project)
        result = self.conn.execute(status_stmt)

        obo_terms = set()
        for hpo_result in result:
            for terms in hpo_result:
                term_list = terms.split(',')
                obos = [i.lower() for i in term_list]
                [obo_terms.add(name) for name in obos]

        return obo_terms

    ## ---------------------------------- ##

    def get_all_samples(self, project, depth='ids'):
        """
        Will return all samples for a given project. 
        depth = [id|all] 
        id: [sample_id], all: [{sample_id:mosaid_sample_id}]
        """
        samples = self.meta.tables['samples']
        status_stmt = select([samples.c.sample_id, samples.c.mosaic_sample_id]).where(samples.c.project == project)
        result = self.conn.execute(status_stmt)
 
        samples = []
        for data in result:
            if depth == 'ids':
                    samples.append(data[0])
            elif depth == 'all':
                samples.append({data[0]:data[1]})

        return samples

    # ---------------------------------- ##

    def get_ready_sample_list(self, project):

        samples = self.meta.tables['samples']
        status_stmt = select([samples.c.sample_id, samples.c.sample_status]).where(samples.c.project == project)
        result = self.conn.execute(status_stmt)
       
        samples = []
        for value in result:
            if value[1] == 'ready':
                samples.append(value[0])
        try:
            if len(samples) == 0:
                raise ValueError
        except ValueError:
            print('ERROR: Could not collect samples for project ', project)
            sys.exit(1)

        ready_list = open('samples.list', 'w')
        for peep in samples:
            who_peep = '{}\n'
            ready_list.write(who_peep.format(peep))
    
    ## ---------------------------------- ##
   
    def get_project_from_accession(self, accession, nicu):
        """
        Will take a assession number and return project associated
        """
        samples = self.meta.tables['samples']
        status_stmt = select([samples.c.project]).where(samples.c.sc_uuid == accession)
        result = self.conn.execute(status_stmt)
        project = [x[0] for x in result]

        try:
            if len(project) == 0:
                raise ProjectCreationError
        except ProjectCreationError:
            nicu.log.info('No project found for accession {}. Use --manifest_check option first.'.format(accession))
            return False
    
        return project.pop()

    ## ---------------------------------- ##
    
    def get_project_id_from_name(self, project):
        """
        Take a given project name and return mosaic_id stored in ucgddb.
        """
        projects = self.meta.tables['projects']
        id_stmt = select([projects.c.mosaic_id]).where(projects.c.project == project)
        result = self.conn.execute(id_stmt)
        mosaic_meta = [x[0] for x in result]
   
        return mosaic_meta[0]
    
    ## ---------------------------------- ##
    
    def get_sample_id_from_accession(self, accession, nicu):
        """
        Will take a assession number and return samples_id associated.
        """
        samples = self.meta.tables['samples']
        status_stmt = select([samples.c.sample_id]).where(samples.c.sc_uuid == accession)
        result = self.conn.execute(status_stmt)
        sample = [x[0] for x in result]

        try:
            if len(sample) == 0:
                raise ProjectCreationError
        except ProjectCreationError:
            nicu.sns.project_issue('No sample found for accession {}. Use --manifest_check option first.'.format(accession))
            nicu.log.error('No sample found for accession {}. Use --manifest_check option first.'.format(accession))
    
        return sample.pop()
    
    ## ---------------------------------- ##

    def get_accession_from_sample_id(self, sample):
        samples = self.meta.tables['samples']
        status_stmt = select([samples.c.sc_uuid]).where(samples.c.sample_id == sample)
        result = self.conn.execute(status_stmt)
        accession = [x[0] for x in result]

        return accession[0]
    
    ## ---------------------------------- ##
    
    def get_iterator_value(self):
        """
        Gets current iterator from db and increases it.
        """
        inter = self.meta.tables['project_iterator']

        status_stmt = select([inter.c.iterator])
        result = self.conn.execute(status_stmt)

        digit = 0
        plus_digit = 0
        for num in result:
            digit = num[0]
            plus_digit = digit + 1

        project_stmt = update(inter).values(iterator=plus_digit)
        self.conn.execute(project_stmt)

        return digit

    ## ---------------------------------- ##
    
    def get_assembly(self, project):
        """
        Will return projects build assembly or reference.
        """
        projects = self.meta.tables['projects']
        id_stmt = select([projects.c.assembly]).where(projects.c.project == project)
        result = self.conn.execute(id_stmt)
        reference = [x[0] for x in result]

        return reference[0]
    
    ## ---------------------------------- ##

    def get_mosaic_project_id(self, project):
        """
        From project name will return mosaic project id.
        """
        projects = self.meta.tables['projects']
        stmt = select([projects.c.mosaic_id]).where(projects.c.project == project)
        result = self.conn.execute(stmt)

        mosaic_result = [x[0] for x in result]
        mosaic_id = mosaic_result.pop()

        return mosaic_id
    
    ## ---------------------------------- ##
   
    def get_individual_sample(self, project, sample):
        """
        Will return tuple of all info for a sample in ucgddb.
        """
        samples = self.meta.tables['samples']
        stmt = select([samples]).where(samples.c.sample_id == sample).where(samples.c.project == project)
        result = self.conn.execute(stmt)

        for send in result:
            return send
    
    ## ---------------------------------- ##
    
    def get_family_type(self, project):
        """
        Will return the collected family type of project.
        """
        projects = self.meta.tables['projects']
        stmt = select([projects.c.family_type]).where(projects.c.project == project)
        result = self.conn.execute(stmt)
        
        family_result = [x[0] for x in result]
        family_type = family_result.pop()
        
        return family_type

    ## ---------------------------------- ##

    def make_ped_file(self, project):
        samples = self.meta.tables['samples']
        status_stmt = select([samples.c.kindred_id, samples.c.sample_id, samples.c.paternal_id, samples.c.maternal_id, samples.c.sex, samples.c.affection_status, samples.c.project]).where(samples.c.project == project)
        execute = self.conn.execute(status_stmt)
  
        ped = project + '.ped'
        ped_file = open(ped, 'w')
        ped_header = "#Kindred_ID\tSample_ID\tPaternal_ID\tMaternal_ID\tSex\tAffection_Status\tProject\n"
        ped_file.write(ped_header)

        for retn in execute:
            if retn[4] == 'male':
                sex = 1
            elif retn[4] == 'female':
                sex = 2
            else:
                sex = 0
    
            if retn[5] == 'affected':
                aff_status = 2
            elif retn[5] == 'unaffected':
                aff_status = 1
            elif retn[5] == 'missing':
                aff_status = 0

            pedigree = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'
            ped_file.write(pedigree.format(retn[0], retn[1], retn[2], retn[3], sex, aff_status, project))
    
    ## ---------------------------------- ##
    
    def make_prep_file(self, project):
        samples = self.meta.tables['samples']
        status_stmt = select([samples.c.sc_uuid, samples.c.sample_id]).where(samples.c.project == project)
        execute = self.conn.execute(status_stmt)

        ## TODO: think of a better way then this.
        all_found_files = []
        project_processing = self.processing_space + project

        for root, dirs, files in os.walk(project_processing):
            for f in files:
                abs_file = root + '/' + f
                if os.path.islink(abs_file):
                    continue
                if re.match(r'.*.fastq.gz$', f):
                    all_found_files.append(abs_file)

        ## create map of accession to sample_id.
        acc_map = {}
        for result in execute:
            try:
                if result[0] == None:
                    raise ProjectError
            except ProjectError:
                print("Requested project not found.")
                sys.exit(1)
            r = re.compile(result[0])
            wanted_list = list(filter(r.search, all_found_files))
            if wanted_list:
                acc_map[result[1]] = wanted_list
   
        prep_file = open('prep_file.txt', 'w')
        for sample, files in acc_map.items():
            for x in files:
                line = '{}\t{}\n'
                prep_file.write(line.format(sample, x))
        prep_file.close()

    ## ---------------------------------- ##
   
    def check_project_status(self):
        """
        Will collect new projects, and check status of each sample.
        if all samples are 'queue' status will update project and samples table
        and return pipeline request. 
        If samples are 'awaiting_data' for longer then 3 days sns message will be sent.
        """
        ## get project elapsed time.
        project_time = ''
        time_now = datetime.datetime.now()

        ## set up projects and samples tables
        projects = self.meta.tables['projects']
        samples = self.meta.tables['samples']

        ## collect all new projects from projects table.
        status_stmt = select([projects.c.project, projects.c.date]).where(projects.c.scope_work == 'NEO').where(projects.c.status == 'new_project')
        results = self.conn.execute(status_stmt)

        ## collect meta data for each sample per project returned.
        project_state = defaultdict(list)
        for n_proj in results:
            project_time = n_proj[1]
            elapsed_time = time_now - project_time

            sample_stmt = select([samples.c.sc_uuid, samples.c.sample_status]).where(samples.c.project == n_proj[0])
            sample_return = self.conn.execute(sample_stmt)

            for result in sample_return:
                meta_d = { 
                    'accession' : result[0],
                    'status' : result[1],
                }
                project_state[n_proj[0]].append(meta_d)

        ## work with collected meta data.
        sample_count, wait_count, queue_count, project_name, wait_samples = 0, 0, 0, '', []
        for project, project_data in project_state.items():
            project_name = project
            for data in project_data:
                sample_count += 1
                if data['status'] == 'queue':
                    queue_count += 1
                elif data['status'] == 'awaiting_data':
                    wait_count += 1
                    wait_samples.append(data['accession'])
        
        ## if all samples are ready, update database and return pipeline request.
        if queue_count == sample_count and queue_count > 0 and sample_count > 0:
            ## update project 
            project_stmt = update(projects).where(projects.c.project == project_name).values(status='ready')
            self.conn.execute(project_stmt)
            
            ## update samples.
            sample_stmt = update(samples).where(samples.c.project == project_name).values(sample_status='ready')
            self.conn.execute(sample_stmt)
   
        if wait_count > 0 and elapsed_time.days > 3:
            samples = ', '.join(wait_samples)
            return True, samples
        else:
            return False, ''
             
    ## ---------------------------------- ##
    
    def check_ready_projects(self):
        """
        Will check project db for NEO ready projects
        verifies sample status and count. returns list of ready projects.
        """
        ## set up projects and samples tables
        projects = self.meta.tables['projects']
        samples = self.meta.tables['samples']

        ## collect all new projects from projects table.
        project_stmt = select([projects.c.project]).where(projects.c.scope_work == 'NEO').where(projects.c.status == 'ready')
        results = self.conn.execute(project_stmt)

        ## Collect count of samples and status.
        projects_state = defaultdict(list)
        for project in results:
            ## update samples.
            sample_stmt = select([samples.c.sample_status]).where(samples.c.project == project[0])
            sample_results = self.conn.execute(sample_stmt)

            sample_count = 0
            ready_count = 0
            for status in sample_results:
                sample_count += 1
                if 'ready' in status[0]:
                    ready_count +=1

            projects_state[project[0]].append(sample_count)
            projects_state[project[0]].append(ready_count)

        ## check per-project.        
        ready_projects = []
        for project, counts in projects_state.items():
            if counts[0] == counts[1]:
                ready_projects.append(project)
        
        return ready_projects
    
    ## ---------------------------------- ##

    def update_ucgddb(self, xlsx_file, nicu, args):
        """
        Will run manifest validation and upload project and sample 
        tables in the UCGDDB database.
        """
        ## make dataframe of excel sheet.
        prep_manifest = pd.read_excel(xlsx_file)
        prep_manifest.columns = prep_manifest.columns.str.lower()

        ## Check for all fields of file & create project name.
        self.__check_manifest_requires(prep_manifest, nicu)
        self.__check_pedigree(prep_manifest, nicu)
        self.__build_nicu_project_name(prep_manifest, args)

        ## Update project and samples table
        self.__update_db_projects(prep_manifest, nicu)
        self.__update_db_samples(prep_manifest, nicu)
        
    ## ---------------------------------- ##
    
    def update_project_status(self, project, new_status):
        projects = self.meta.tables['projects']
        stmt = update(projects).where(projects.c.project == project).values(status = new_status)
        self.conn.execute(stmt)
    
    ## ---------------------------------- ##
    
    def update_sample_status(self, accession, new_status):
        samples = self.meta.tables['samples']
        stmt = update(samples).where(samples.c.sc_uuid == accession).values(sample_status= new_status)
        self.conn.execute(stmt)

    ## ---------------------------------- ##

    def __check_manifest_requires(self, manifest, nicu):
        manifest.fillna(0, inplace=True)
        failed = False
        for x, content in manifest.iterrows():
            for req in self.required_fields:
                try:
                    ## check that all required value are present other than
                    ## paternal and maternal requirements.
                    if not ('ternal' in req):
                        if content.get(req) == 0:
                            message = 'Required field {} requires value, none given.'
                            raise ManifestError(message.format(req))
                except ManifestError as e:
                    failed = True
                    nicu.log.info(e)
                    nicu.sns.manifest_issue(e)
                    sys.exit(1)
                try:
                    if (content.get(req) == 0 and content.get('Sample_ID')):
                        error_m = "Needed value {} for sample {} not given."
                        raise ManifestError(error_m.format(req, content.get('Sample_ID')))
                except ManifestError as e:
                    failed = True
                    nicu.log.error(e)
                    nicu.sns.manifest_issue(e)
                    sys.exit(1)

        if failed == True:
            sys.exit(1)

    ## ---------------------------------- ##

    def __check_pedigree(self, manifest, nicu):
        ## update any missing values with 0.
        manifest.fillna(0, inplace=True)

        for index, data in manifest.iterrows():
            if not data.get('kindred_id'):
                nicu.sns.manifest_issue('kindred_id missing.')
                nicu.log.error('kindred_id missing.')

            if not data.get('affection_status'):
                nicu.sns.manifest_issue('affection_status missing.')
                nicu.log.error('affection_status missing.')

            if not data.get('sex'):
                nicu.sns.manifest_issue('sample sex info missing.')
                nicu.log.error('sample sex info missing.')

            if data.get('maternal_id') == data.get('sample_id'):
                nicu.sns.manifest_issue('maternal_id missing or incorrect.')
                nicu.log.error('maternal_id missing or incorrect.')

            if data.get('paternal_id') == data.get('sample_id'):
                nicu.sns.manifest_issue('paternal_id missing.')
                nicu.log.error('paternal_id missing.')

    ## ---------------------------------- ##

    def __build_nicu_project_name(self, manifest, args):

        if args.test:
            self.project_name = 'T101-101-NEO-Cyberdyne'
            return
        else:
            interator = self.get_iterator_value()
            ucgd_number = 'A' + str(interator)

            kindred_id = manifest.get('kindred_id')[1]
            phenotype_list = manifest.get('hpo_terms')[1].split(',')
            hpo_term = phenotype_list[0].replace(' ', '-')

            name_format = "{}-{}-NEO-{}"
            project_name = name_format.format(ucgd_number, kindred_id, hpo_term.capitalize())
            self.project_name = project_name

    ## ---------------------------------- ##

    def __update_db_projects(self, manifest, nicu):
        now = datetime.datetime.now()
        projects = self.meta.tables['projects']

        entry = {
            'project': self.project_name,
            'gnomex_analysis_group': self.project_name,
            'gnomex_year': now.year,
            'data_namespace': 'UCGD',
            'description': None,
            'project_manager_last_name': 'Boyden',
            'project_manager_first_name': 'Steven',
            'start_date': now,
            'primary_data_path': '/scratch/ucgd/lustre/UCGD_Datahub/IRBs/IRB_00125940',
            'scheduled_start': now,
            'vc_framework': 'Sentieon',
            'gnomex_data_path': '/scratch/ucgd/lustre/UCGD_Datahub/Repository/AnalysisData',
            'date': now,
            'status': 'new_project',
            'scope_work': 'NEO',
            'seq_design': 'rWGS',
            'sequence_center': 'ARUP',
            'box_url': 'https://uofu.app.box.com/folder/79662845561',
            'assembly': manifest.get('reference')[1],
            'irb_institution' : manifest.get('irb_institution')[1],
            'irb_number' : manifest.get('irb_number')[1],
            'notes': manifest.get('notes')[1],
            'pi_first_name': manifest.get('pi_first_name')[1],
            'pi_last_name': manifest.get('pi_last_name')[1],
            'phenotype': manifest.get('phenotype_description')[1],
            'family_type': manifest.get('family_type')[1],
        }
        try:
            self.conn.execute(projects.insert(), entry)
        except:
            print(self.conn.error)
            nicu.log.error('Project {} updating issue or project exists.'.format(self.project_name))
            sys.exit(1)


    ## ---------------------------------- ##

    def __update_db_samples(self, manifest, nicu):
        samples = self.meta.tables['samples']

        for x, content in manifest.iterrows():
            ## check for and remove '-' in arup accession id
            arup_accession = content.get('arup_accession')
            if isinstance(arup_accession, str):
                if arup_accession.find('-'):
                    arup_accession = re.sub(r'-', '', arup_accession)
            entry = {
                'sample_id' : content.get('sample_id'),
                'sc_uuid' : arup_accession,
                'project' : self.project_name,
                'pi_first_name' : content.get('pi_first_name'),
                'pi_last_name' : content.get('pi_last_name'),
                'irb_institution' : content.get('irb_institution'),
                'irb_number' : content.get('irb_number'),
                'kindred_id' : content.get('kindred_id'),
                'paternal_id' : content.get('paternal_id'),
                'maternal_id' : content.get('maternal_id'),
                'sex' : content.get('sex'),
                'affection_status' : content.get('affection_status'),
                'phenotype_description' : content.get('phenotype_description'),
                'hpo_terms' : content.get('hpo_terms'),
                'birth_year' : content.get('birth_year'),
                'assessment_year' : content.get('assessment_year'),
                'ancestry' : content.get('ancestry'),
                'race' : content.get('ethnic_group'),
                'consanguinity' : content.get('consanguinity'),
                'tissue_type' : content.get('tissue_type'),
                'tissue_condition' : content.get('tissue_condition'),
                'molecule_type' : content.get('molecule_type'),
                'sequence_center' : content.get('sequence_center'),
                'seq_design' : content.get('seq_design'),
                'sample_status' : 'awaiting_data',
                'capture' : content.get('capture'),
                'library_kit' : content.get('library_kit'),
                'library_pcr' : content.get('library_pcr'),
                'instrument' : content.get('instrument'),
                'target_depth' : content.get('target_depth'),
                'incidental_consent' : content.get('incidental_consent'),
                'notes' : content.get('notes'),
            }
            try:
                self.conn.execute(samples.insert(), entry)
            except:
                print(self.conn.error)
                nicu.log.error('Sample for project {} updating issue or sample[s] exists'.format(self.project_name))
                sys.exit(1)

    ## ---------------------------------- ##
