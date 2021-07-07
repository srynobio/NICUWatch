#!/usr/bin/env python3

import os
import sys
import logging
import datetime

class LogUCGD:

    def __init__(self, args):
        self.log_file = args.log_file
    
    ## ---------------------------------- ##

    def get_current_date(self):
        return datetime.datetime.now()

    ## ---------------------------------- ##
    
    def info(self, message):
        logging.basicConfig(filename=self.log_file, level=logging.INFO)
        date = self.get_current_date()
        mesg  = '{} : {}'
        logging.info(mesg.format(date, message))
    
    ## ---------------------------------- ##

    def error(self, message):
        logging.basicConfig(filename=self.log_file, level=logging.ERROR)
        date = self.get_current_date()
        mesg  = '{} : {}'
        logging.info(mesg.format(date, message))
        sys.exit(1)
    
    ## ---------------------------------- ##

