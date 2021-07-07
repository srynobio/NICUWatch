#!/usr/bin/env python3

import os
import sys
import boto3
import argparse
import hashlib
import datetime
import json
import configparser
import datetime

class SNS:

    def __init__(self, args):
        ## Open found configure file.
        config = configparser.ConfigParser()
        config.read(args.config)

        # database info
        aws_access_key_id = config['aws']['aws_access_key_id']
        aws_secret_access_key = config['aws']['aws_secret_access_key']
        region_name = config['aws']['region_name']

        # Daily Message recorder.
        today = datetime.date.today()
        self.record_file = args.message_tracker_path + str(today) + '_message_tracker.txt'
       
        ## create boto3 SNS client.
        client = boto3.client(
                "sns",
                aws_access_key_id = aws_access_key_id,
                aws_secret_access_key = aws_secret_access_key,
                region_name = region_name 
        )
        self.client = client 

        ## add if testing.
        self.test = True if args.test else False

    ## ------------------------------ ##
   
    def message_sent(self, message):
        """Returns True if message has already been sent."""
        m = hashlib.sha256()
        m.update(message.encode('utf-8'))
        message_hash = m.hexdigest()

        ## open or create record file.        
        if os.path.exists(self.record_file):
            ledger = open(self.record_file, 'r+') 
        else: 
            ledger = open(self.record_file, 'w+') 
        
        for f in ledger:
            f_clean = f.replace('\n', '')
            if f_clean == message_hash:
                return True

        print(message_hash, file=ledger)
        ledger.close()
        return False

    ## ------------------------------ ##

    def neoseq_test(self, message):
        topic = self.client.create_topic(Name='NeoSeq_Test')
        topic_arn = topic['TopicArn']

        alert = """
UCGD NeoSeq Test Message: {}
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)
    
    ## ------------------------------ ##
    
    def portal_issue(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='Portal_Issue')
        topic_arn = topic['TopicArn']

        alert = """
UCGD Portal Alert: {}
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)
 
    ## ------------------------------ ##

    def neoseq_general(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='NeoSeq_General')
        topic_arn = topic['TopicArn']

        alert = """
UCGD NeoSeq General Message: {}
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)
    
    ## ------------------------------ ##
    
    def portal_issue(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='Portal_Issue')
        topic_arn = topic['TopicArn']

        alert = """
UCGD Portal Alert: {}
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)
    
    ## ------------------------------ ##

    def time_tracker(self, project, step):
        if self.test == True:
            return
        topic = self.client.create_topic(Name='Time_Tracker')
        topic_arn = topic['TopicArn']

        ## get current time.
        current_time = datetime.datetime.now()

        alert = """
Time tracker message for project: {}
The following step has completed: {} at {}
        """
        tweet = alert.format(project, step, current_time)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)

    ## ------------------------------ ##
    
    def project_issue(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='Project_Issue')
        topic_arn = topic['TopicArn']

        alert = """
Project issue discovered: {}
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)

    ## ------------------------------ ##

    def manifest_issue(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='Manifest_Issue')
        topic_arn = topic['TopicArn']

        alert = """
Manifest Issue discovered: {}
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)

    ## ------------------------------ ##
    
    def new_project_discovered(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='New_Project')
        topic_arn = topic['TopicArn']

        alert = """
Project {} discovered. Processing will begin when data arrives ~48 hours.
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)

    ## ------------------------------ ##

    def rouge_data(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='Rouge_Data')
        topic_arn = topic['TopicArn']

        alert = """
ARUP drop location discovered unknown or non-message data: {}
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)

    ## ------------------------------ ##

    def late_sample(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='Late_Sample')
        topic_arn = topic['TopicArn']

        alert = """
Awaiting (> 72 hours) accession data {}.'''
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)
 
    ## ------------------------------ ##

    def checksum_issue(self, accession, md5_file):
        if self.test == True:
            self.neoseq_test(md5_file)
            return
        topic = self.client.create_topic(Name='Check_Sum_Issue')
        topic_arn = topic['TopicArn']

        alert = """
Fastq MD5 validation failed for samples {}
        """
        tweet = alert.format(md5_file)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)

    ## ------------------------------ ##
   
    def missing_checksum(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='Check_Sum_Issue')
        topic_arn = topic['TopicArn']

        alert = """
Missing MD5 checksum file for message sample set: {}
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)

    ## ------------------------------ ##

    def processing_started(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='Processing_Started')
        topic_arn = topic['TopicArn']

        alert = """
Processing for message {} has started.
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)

    ## ------------------------------ ##
    
    def data_alert(self, project, step):
        if self.test == True:
            return
        topic = self.client.create_topic(Name='Data_Alert')
        topic_arn = topic['TopicArn']

        alert = ''
        if step == 'VAR':
            alert = """
BAM and Final VCF called for project {}
            """
        elif step == 'SV':
            alert = """
Smoove VCF file complete for project {}
            """
        elif step == 'POST':
            alert = """
Post VCF process and upload to Fabric complete for project {}
            """
        elif step == 'FASTQ':
            alert = """
Validated fastq for project {} ready.
            """
        elif step == 'MOSAIC':
            alert = """
Mosaic project {} ready.
            """
        tweet = alert.format(project)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)
    
    ## ------------------------------ ##

    def processing_completed(self, project, root):
        if self.test == True:
            self.neoseq_test('test processing comlete.')
            return
        topic = self.client.create_topic(Name='Processing_Complete')
        topic_arn = topic['TopicArn']

        project_path = root + project
        alert = """
Project {} has finalized processing.
Project path: {}
        """
        tweet = alert.format(project, project_path)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)
    
    ## ------------------------------ ##

    def fabric_issue(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='Fabric_Issue')
        topic_arn = topic['TopicArn']

        alert = """
Issue with Fabric api.
{}
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)

    ## ------------------------------ ##

    def ontology_issue(self, message):
        if self.test == True:
            self.neoseq_test(message) 
            return
        topic = self.client.create_topic(Name='Ontology_Issue')
        topic_arn = topic['TopicArn']

        alert = """
Missing Ontology term from HPO.
{}
        """
        tweet = alert.format(message)
        if self.message_sent(tweet) == False:
            self.client.publish(Message=tweet, TopicArn=topic_arn)

    ## ------------------------------ ##

