#!/usr/bin/env python3

import requests
from requests.auth import HTTPBasicAuth
import json


class UCGDPortal:

    def __init__(self, args):
        self.portal_url = '**********************'

    ## ---------------------------------- ##

    def portal_health(self, nicu):
        """Will collect all current mosaic projects and return json"""
        path = '{}'
        url = path.format(self.portal_url)

        try:
            response = requests.get(url)
        except:
            nicu.sns.portal_issue('The UCGD Neoseq API portal appears to be down.')
            nicu.log.info('The UCGD Neoseq API portal appears to be down.')
    
    ## ---------------------------------- ##
