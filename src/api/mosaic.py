#!/usr/bin/env python3

import os
import sys
import requests
import json
import configparser
from requests.auth import HTTPBasicAuth
            
class ProjectCreationError(Exception):
    pass

class RolesAddError(Exception):
    pass

class ConnectionError(Exception):
    pass

class mosaicAPI:

    def __init__(self, args):
        ## open config file.
        config = configparser.ConfigParser()
        config.read(args.config)

        # database info
        mosaic_username = config['ucgd-mosaic-cli']['MOSAIC_USERNAME']
        mosaic_password = config['ucgd-mosaic-cli']['MOSAIC_PASSWORD']

        ## create request auth
        mosaic_bearer = config['ucgd-mosaic']['authorization']

        mosaic_auth = {
            'url' :'**************************',
            'authorization': mosaic_bearer,
            'content-type': 'application/json',
        }

        ## set default project mosaic attribures
        ## Mosaic Attribute id : Attribute default value
        default_attributes = {
            123 : 'Undiagnosed', #diagnostic_status
            122 : 'In-progress', #processing_status
            121 : 'Trio',
        }            

        ## set user and roles.
        ## [mosaic id : role type]
        if args.test:
            default_users = {
                # Removed
            }
        }
        else:
            default_users = {
                # Removed
            }
        self.default_users = default_users
        self.mosaic_auth = mosaic_auth
        self.default_attributes = default_attributes

    ## ------------------------------------------------------------------- ##

    def get_all_projects(self, args, nicu):
        """Will collect all current mosaic projects and return json"""
        path = "{}/projects"
        url = path.format(self.mosaic_auth['url'])

        headers = {
            'Authorization' : self.mosaic_auth['authorization'],
            'Content-Type': 'application/json',
        }
        params = {
            'limit' : 100
        }
        result = requests.get(url, headers=headers, params=params)

        if not result.ok:
            nicu.log.info('Could not get all mosaic projects {}.'.format(args.mosaic_project_id))
            nicu.log.error('Mosaic get all projects failed due to: {}.'.format(result.text))
   
        nicu.log.info('Mosaic API get_all_projects response: {}'.format(result.text))

        return result.json()

    ## ------------------------------------------------------------------- ##
    
    def create_project(self, project, nicu):
        path = "{}/projects"
        url = path.format(self.mosaic_auth['url'])
        description = 'UCGD NeoSeq project {}'.format(project)
        reference = nicu.ucgd_db.get_assembly(project)

        header = {
            'Authorization' : self.mosaic_auth['authorization'],
        }
        info = {
            'name' : project,
            'description' : description,
            'reference' : reference, 
        }
        result = requests.post(url, headers=header, data=info) 

        ## check if project created.
        try:
            if result.status_code != 200:
                raise ProjectCreationError
        except ProjectCreationError as e:
            nicu.log.error('Could not create Mosaic project {} {}'.format(project, e))

        nicu.log.info('Mosaic API create_project response: {}'.format(result.text))

        ## get project mosaic id
        project_info = result.json().get('project_settings')

        ## Add mosaic id to ucgddb project table.
        nicu.ucgd_db.add_mosaic_project_id(project_info['project_id'], project, nicu)

        ## add correct roles for users.
        for i, r_type in self.default_users.items():
            self.add_role(i, r_type, project_info['project_id'], nicu)
    
        ## add attributes defaults.
        self.add_default_project_attributes(project_info['project_id'], nicu)
        self.add_project_attributes_to_dashboard(project_info['project_id'], nicu)

    ## ------------------------------------------------------------------- ##
    
    def add_default_project_attributes(self, mosaic_project_id, nicu):
        """
        Will add selected default attributes to the project
        """
        for id, value in self.default_attributes.items():
            path = "{}/projects/{}/attributes/import"
            url = path.format(self.mosaic_auth['url'], mosaic_project_id)

            headers = {
                'Authorization' : self.mosaic_auth['authorization'],
            }
            data = {
                'attribute_id' : id,
                'value' : value,
            }
            result = requests.post(url, data=data, headers=headers)

            if not result.ok:
                nicu.log.info('Could not add default mosaic project attributes. ERROR: {}'.format(result.text))

            nicu.log.info('Mosaic API add_default_project_attributes response: {}'.format(result.text))

    ## ------------------------------------------------------------------- ##

    def add_project_attributes_to_dashboard(self, mosaic_project_id, nicu):
        """ 
        Will pin default attributes to the dashboard of the project,
        """
        for id, value in self.default_attributes.items():
            path = "{}/projects/{}/dashboard"
            url = path.format(self.mosaic_auth['url'], mosaic_project_id)

            headers = {
                'Authorization' : self.mosaic_auth['authorization'],
            }
            data = {
                'type' : 'project_attribute',
                'is_active' : 'true',
                'attribute_id' : id,
            }
            result = requests.post(url, data=data, headers=headers)

            if not result.ok:
                nicu.log.info('Could not add default mosaic project attributes. ERROR: {}'.format(result.text))

            nicu.log.info('Mosaic API add_project_attributes_to_dashboard response: {}'.format(result.text))

    ## ------------------------------------------------------------------- ##

    def add_role(self, user, role_type, project_id, nicu):
        path = "{}/projects/{}/roles"
        url = path.format(self.mosaic_auth['url'], project_id) 
        
        headers = {
            'Authorization' : self.mosaic_auth['authorization'],
        }
        data = {
                'user_id' : user,
                'role_type_id' : role_type,
        }
        params = {
                'project_id' : project_id,
        }
        result = requests.post(url, headers=headers, data=data, params=params)

        ## check if project roles added.
        ## non fatal as can be added later.
        try:
            if result.status_code != 200:
                raise RolesAddError
            ProjectCreationError
        except RolesAddError as e:
            nicu.log.info('Could not add role for user {} {}'.format(user, e))
        
        nicu.log.info('Mosaic API add_role response: {}'.format(result.text))


    ## ------------------------------------------------------------------- ##

    def post_slivar_file(self, nicu, args):
        """Updates mosaic with provided slivar file."""
        path = "{}/projects/{}/variants/upload"
        url = path.format(self.mosaic_auth['url'], args.mosaic_project_id)

        files = {'file' : open(args.slivar_file, 'rb') }
        headers = {
            'Authorization' : self.mosaic_auth['authorization'],
        }
        result = requests.post(url, headers=headers, files = files)

        if not result.ok:
            nicu.log.info('Could not add slivar file to mosaic project {}.'.format(args.mosaic_project_id))
            nicu.log.info('Mosaic upload of slivar file due to: {}.'.format(result.text))

        nicu.log.info('Mosaic API post_slivar_file response: {}'.format(result.text))

    ## ------------------------------------------------------------------- ##

    def post_pedigree(self, nicu, args):
        """
        Will post pedigree data to mosaic project.
        """
        path = "{}/projects/{}/pedigree"
        url = path.format(self.mosaic_auth['url'], args.mosaic_project_id)

        headers = {
            'Authorization' : self.mosaic_auth['authorization'],
        }
        params = {
            'project_id' : args.mosaic_project_id,
        }            
        files = {'file': (args.pedigree_file, open(args.pedigree_file, 'r'),'text/csv')}

        result = requests.post(url, headers=headers, params=params, files=files)

        if not result.ok:
            nicu.log.info('Could not add pedigree file to mosaic project {}.'.format(args.mosaic_project_id))
            nicu.log.info('Mosaic pedigree upload error due to: {}.'.format(result.text))

        nicu.log.info('Mosaic API post_pedigree response: {}'.format(result.text))

    ## ------------------------------------------------------------------- ##
