#!/usr/bin/env python3

import os
import sys
import requests
from requests.auth import HTTPBasicAuth
import json
import configparser

class fabricAPI:

    def __init__(self, args):
        ## open config file.
        config = configparser.ConfigParser()
        config.read(args.config)

        if args.test:
        # fabric dev api
            fabric_login = config['fabric-dev-nicu-api']['login']
            fabric_pass  = config['fabric-dev-nicu-api']['password']
            fabric_url = '*******************'
            fabric_panel_id = '***************'
        else:
            # fabric production api
            fabric_login = config['fabric-production-nicu-api']['login']
            fabric_pass  = config['fabric-production-nicu-api']['password']
            fabric_url = '******************'
            fabric_panel_id = '*****************'

        ## create request auth
        auth = HTTPBasicAuth(fabric_login, fabric_pass)

        self.fabric_url = fabric_url
        self.fabric_panel_id = fabric_panel_id
        self.auth = auth
    
    ## ------------------------------------- ##

    def create_project(self, project_name, nicu):
        """Will create a new UCGD project via fabric api."""
        path = "{}/projects/"
        url = path.format(self.fabric_url)

        description = 'UCGD NeoSeq project {}'.format(project_name)
        project_payload = {
            'project_name': project_name,
            'description' :  description,
            'share_role'  : 'CONTRIBUTOR',
        }
        result = requests.post(url, data=project_payload, auth=self.auth)

        if not result.ok:
            nicu.sns.fabric_issue('Could not create fabric project.'.format(nicu.args.project))
            nicu.log.info('Could not create fabric project {}.'.format(nicu.args.project))
            nicu.log.info('Fabric create project error response: {}.'.format(result.text))

        nicu.log.info('Fabric API create_project response: {}'.format(result.text))

        ## Add fabric project id to database.
        json_response = result.json()
        project_id = json_response.get('id')
        nicu.ucgd_db.add_fabric_project_id(project_name, project_id)

    ## ------------------------------------- ##

    def get_all_projects(self):
        """Will collect all fabric id and projects and return them as a 
        dict -> {id:project}"""
        path = "{}/projects/"
        url = path.format(self.fabric_url)
        result = requests.get(url, auth=self.auth)

        result = result.json()

        id_project = {}
        for found in result['objects']:
            id_project[found['id']] = found['project_name']
        return id_project
    
    ## ------------------------------------- ##

    def add_genome(self, fabric_project_id, sex, vcf_file, sample_name, project_name, checksum, nicu):
        """Will add a vcf genome file to an existing fabric project and update ucgddb with genomic id."""
        path = "{}/projects/{}/genomes?genome_label={}&genome_sex={}&assembly_version=b38&checksum={}"
        url = path.format(self.fabric_url, fabric_project_id, sample_name, sex, checksum)

        with open(vcf_file, 'rb') as file_handle:
            # Post request and return id of newly uploaded genome
            result = requests.put(url, auth=self.auth, data=file_handle)
            if not result.ok:
                nicu.sns.fabric_issue('Could not load Slivar VCF sample {} to fabric projct.'.format(sample_name))
                nicu.log.info('Could not load Slivar VCF sample {} to fabric project.'.format(sample_name))
                nicu.log.info('Fabric add_genome error response: {}.'.format(result.text))

            ## get genome id to add to ucgddb. 
            original_response = result.json()
            genome_id = original_response.get('genome_id')
            sample_name = original_response.get('genome_label')
            nicu.ucgd_db.add_fabric_genome_id(genome_id, sample_name, project_name, nicu)

        nicu.log.info('Fabric API add_genome response: {}'.format(result.text))

    ## ------------------------------------- ##
   
    def delete_fabric_project(self, fabric_id):
        """Given a fabric project id will delete project from Opal"""
        path = "{}/projects/{}"
        url = path.format(self.fabric_url, fabric_id)
        result = requests.delete(url, auth=self.auth)
    
    ## ------------------------------------- ##

    def create_case(self, report_json, args, nicu):
        """Take created json report and posts to fabric via api."""

        path = "{}/reports/"
        url = path.format(self.fabric_url)
        result = requests.post(url, auth=self.auth, data=report_json)

        if not result.ok:
            nicu.sns.fabric_issue('Could not create case for project {}.'.format(args.project))
            nicu.log.info('Could not create case for project {}.'.format(args.project))
            nicu.log.info('Fabric create case failed due to {}.'.format(result.text))
            nicu.log.info('Upload case issue: {}'.format(result.text))

        ## record payload and response.
        nicu.log.info('Fabric create case payload: {}'.format(report_json))
        nicu.log.info('Fabric API create_case response: {}'.format(result.text))

    ## ------------------------------------- ##
