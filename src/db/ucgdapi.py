#!/usr/bin/env python3

import os
import sys
import sqlalchemy as db
import configparser
from sqlalchemy.sql import table, column, select

class ProcessSpaceError(Exception):
    pass

class PortalUCGD:

    def __init__(self, args):

        ## Open found configure file.
        config = configparser.ConfigParser()
        config.read(args.config)

        # database info
        dbname   = config['webportal']['dbname']
        host     = config['webportal']['host']
        port     = config['webportal']['port']
        username = config['webportal']['username']
        password = config['webportal']['password']

        url = 'postgresql://{}:{}@{}:{}/{}'
        url = url.format(username, password, host, port, dbname)
        conn = db.create_engine(url, client_encoding='utf8')
        meta = db.MetaData(bind=conn, reflect=True)

        ## put hardcoded path in object.
        self.nicu_drop = args.nicu_drop
        self.processing_space = args.irb_path

        self.conn = conn
        self.meta = meta

    ## ------------------------------------ ##

    def get_current_accession_ids(self):

        transfer = self.meta.tables['portal_api_datatransfer']
        statement = select([transfer.c.src_id])
        result = self.conn.execute(statement)
        
        src_ids = []
        for i in result.fetchall():
            src_ids.append(i[0])
        
        return src_ids

    ## ------------------------------------ ##
    
    def move_drop_data(self, accession, project, nicu):
        """Collect and move fastq files to processing location."""

        transfer_id = 'ARUP-' + accession
        transfer_file = self.meta.tables['portal_api_datatransferfile']
        statement = select([transfer_file.c.filename]).where(transfer_file.c.datatransfer_id == transfer_id)
        result = self.conn.execute(statement)
        drop_location = nicu.ucgd_db.get_nicu_drop()

        ## project processing space.
        project_setup = self.processing_space + project + '/Project_Setup/'
        try:
            if not (os.path.exists(project_setup)):
                raise ProcessSpaceError
        except ProcessSpaceError as e:
            nicu.sns.project_issue('Could not move accession data, project {} Project_Setup directory not found, or data received before manifest added.'.format(project))
            nicu.log.error('Could not move accession data, project {} Project_Setup directory not found, or data received before manifest added.'.format(project))

        for files in result:
            acc_dir = drop_location + accession
            drop_file = drop_location + files[0]
            drop_md5 = drop_location + files[0] + '.md5'
            file_meta = files[0].split('/')
            new_file = project_setup + file_meta[-1] 
            new_md5  = project_setup + file_meta[-1] + '.md5'

            if not (os.path.exists(new_file)):
                os.rename(drop_file, new_file)
                os.rename(drop_md5, new_md5)
                if not os.listdir(acc_dir):
                    os.rmdir(acc_dir)
                nicu.log.info('File: {} moved into processing space for project: {}'.format(drop_file, project))
            
    ## ------------------------------------ ##

