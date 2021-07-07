#!/usr/bin/env python3

import os
import sys
import configparser
import networkx
import obonet

class OBO:

    def __init__(self, args):
        ## open config file.
        config = configparser.ConfigParser()
        config.read(args.config)

        if args.hpo_file:
            self.hpo_file = args.hpo_file
        else:
            self.hpo_file = None

    ## ------------------------------------- ##

    def term_map_parse(self):
        """Will open and parse obo file to name -> id map"""
        try:
            if self.hpo_file == None:
                raise ValueError
        except ValueError:
            print('hpo.obo file required to create fabric trio case.')
            sys.exit(1)

        graph = obonet.read_obo(self.hpo_file)
        name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}
        return name_to_id
    
    ## ------------------------------------- ##

    def term_to_id(self, project_hpo, fabric_project_id, kindred_id, family_type, nicu):
        id_map = self.term_map_parse()

        term_mapped = {}
        for name, hpo_term in id_map.items():
            fitted_name = name.lower().replace(' ', '_')
            term_mapped[fitted_name] = hpo_term

        ## map project term to parsed hpo file,
        to_hpo_ids = []
        no_match = []
        for map_set in project_hpo:
            if map_set in term_mapped:
                to_hpo_ids.append(term_mapped[map_set])
            else:
                no_match.append(map_set)

        ## send sns of non-match terms
        if len(no_match) < 0:
            nicu.sns.ontology_issue('Matching hpo term could not be found for: {}'.format(', '.join(no_match)))

        ## create template dict entry.
        ## collect family type from db.
        fabric_terms = ', '.join(to_hpo_ids)
        report_type = {
            "report_type": family_type,
            "accession_id": kindred_id,
            "project_id": fabric_project_id,
            "panel_id": nicu.f_api.fabric_panel_id,
            "phevor_terms": fabric_terms,
        }
        return report_type
        
    ## ------------------------------------- ##
