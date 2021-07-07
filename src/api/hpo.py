#!/usr/bin/env python3

import os
import sys
import re
import requests
from requests.auth import HTTPBasicAuth
import json
import argparse
import urllib3

class OBOApi:

    def __init__(self):
        obo_url = 'https://hpo.jax.org'
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

        ## create request auth, hpo api does not require username and password.
        auth = HTTPBasicAuth('','')

        self.obo_url = obo_url
        self.auth = auth


    ## ------------------------------------- ##

    def api_up_check(self):
        """General term to check if api is up."""
        path = '{}/api/hpo/term/HP:0000006'
        url = path.format(self.obo_url)
        result = requests.get(url, auth=self.auth)

        if not result.ok:
            print('The HPO API appears to be down.')

    ## ------------------------------------- ##

    def get_obo_term(self, obo_terms):
        """Will search hpo api for like terms and return dict of found terms:id and dict of undiscovered terms."""
        discovered_terms = {}
        undiscovered_terms = {}
        for obo in obo_terms:
            lc_obo = obo.lower().replace(' ', '_')
            path = '{}/api/hpo/search/?q={}'
            url = path.format(self.obo_url, lc_obo) 
            result = requests.get(url, auth=self.auth, verify=False)

            result_json = result.json()

            if not result_json['terms']:
                undiscovered_terms[obo] = None

            for found in result_json['terms']:
                lc_name = found['name'].lower().replace(' ', '_')
                if re.fullmatch(lc_obo, lc_name):
                    discovered_terms[lc_name] = found['id']

        return discovered_terms, undiscovered_terms

    ## ------------------------------------- ##

    def get_genes_from_id(self, obo_id):
        path = '{}/api/hpo/term/{}/genes?max=-1&offset=1'
        url = path.format(self.obo_url, obo_id)
        result = requests.get(url, auth=self.auth, verify=False)

        result = result.json()

        gene_list = []
        for collect in result['genes']:
            for key, value in collect.items():
                if key == 'entrezGeneSymbol':
                    gene_list.append(value)

        return gene_list

    ## ------------------------------------- ##
