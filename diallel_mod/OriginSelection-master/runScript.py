# Generate diallel table 
# Yujing Cao, July 22, 2021

import pandas as pd
import numpy as np
import requests
import s3fs
import MidParentValue as MP
import rulefilter as rf
import utils
import master_table as mt
import collections
import re
import os
import warnings
warnings.filterwarnings("ignore")
import dask
from dask import delayed
import dask.dataframe as dd
from functools import reduce
import itertools
import time
import json
import hvac
import boto3
import argparse
import logging
from botocore.exceptions import ClientError

def parse_argument():
    parser = argparse.ArgumentParser(description='parse arguments for materialized view.')
    parser.add_argument('--gcafile', dest='gcafile', help='path to a gca master table.', required = True)
    parser.add_argument('--bucket', dest='bucket', help='S3 bucket where gca master table is saved.', required = True)
    parser.add_argument('--scafile', dest='scafile', help='path to a sca master table. This is optional input.')
    parser.add_argument('--keyinbred', dest='keyinbred', help='path to key inbred lists provide by breeders. The inbred list should be saved as csv file. This is optional input.')
    parser.add_argument('--het', dest='het', nargs = "+", help='heterotic group for filtering. This is optional input.')
    parser.add_argument('--config', dest='config', help='path to a configuration file where filtering criteria is stored as a dictionary in a txt file.')
    parser.add_argument('--local_path', dest='local_path', help='path to save output locally.', required = True)
    parser.add_argument('--s3_path', dest='s3_path', help='path to write output to s3')
    
    args = parser.parse_args()
    
    return args

def main(args):
    
    gcafilename = args.gcafile
    bucket = args.bucket
    
    if args.scafile is not None:
        scafielname = args.scafile
    
    if args.scafile is not None:
        scafielname = args.scafile
        
    if args.keyinbred is not None:
        path_keyinbred = args.keyinbred
        key_inbred = pd.read_csv(path_keyinbred)
        
    if args.het is not None:
        het = args.het
        print(het)
#           het = args.het
        
        
    if args.config is not None:
        path_config = args.config
        with open(path_config) as f:
            DEFAULTS = json.loads(f.read())
    
    ######## Load Parental Table ##########
    print('\n========= Loading Parental Table =========')
    parentalTable = utils.read_s3_file_generic(filename = 'Cucumber/Cucumber_Parentals.csv', bucket = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data')
    parentalTable['Pedigree_1'] = parentalTable['FemalePedigree'] + '__' + parentalTable['MalePedigree']
    parentalTable['Pedigree_2'] = parentalTable['MalePedigree'] + '__' + parentalTable['FemalePedigree']

    ######### Load master file ############
    print('\n========= Loading GCA Master Table =========')
    if 'ec2' in os.environ['HOME']:
        gca_dat = utils.read_as_dask_1(gcafilename, bucket)
        
        # add index
        gca_dat =  gca_dat.compute().set_index(pd.Index(range(gca_dat.compute().shape[0])))
    else:
        gca_dat = utils.read_as_dask(gcafilename, bucket)
        gca_dat =  gca_dat.compute().set_index(pd.Index(range(gca_dat.compute().shape[0])))
#     print(gca_dat.head())
        
    ######### Filter GCA table by key inbred list if it is available ########
    
    if args.keyinbred is not None:
        print('\n========= Filtering GCA Master Table By Key Inbred Lines =========')
        gca_dat_filtered = gca_dat[gca_dat['PEDIGREE_NAME'].isin(key_inbred['Pedigree'])]
    else: 
        gca_dat_filtered = gca_dat
    
    ######### filter GCA table by heterotic group if it is available #########
    if args.het is not None:
        print('\n========= Filtering GCA Master Table By Heterotic Group =========')
        gca_dat_filtered = rf.filter_by_HETGP(gca_dat_filtered, het=het)
#         gca_dat_filtered = rf.filter_by_HETGP(gca_dat_filtered, het=['678', '83', '89', 'A', 'SLM804'])
#     print(gca_dat_filtered.head())
    
    ######### Extract traits
    
    print('\n========= Extract Focal Traits =========')
    
#     print(gca_dat_filtered.head())
    traits = [col for col in gca_dat_filtered.columns if 'predicted.value_' in col]
#     print(traits)
    traits_1 = []
    for tt in traits: 
        t_1 = re.sub('predicted.value_', '', tt)
        traits_1.append(t_1)
#     print(traits_1)
    
    ######### Calculate mid parent value 
    print('\n========= Calculating Mid-parent Values for All the Focal Traits =========')
    diallel_table = MP.gather_mid_value(dd.from_pandas(gca_dat_filtered, npartitions=10), traits_1)
    
    
    ############### Get reference data in one table
    
    print('\n========= Getting Reference Data in One Table =========')
    reference_cols = [col for col in gca_dat_filtered.columns if not any(tr in col for tr in traits_1)]
    reference_table = gca_dat_filtered[reference_cols]
    if 'PEDIGREE' in reference_table:
        reference_table = reference_table.drop(columns=['PEDIGREE'], axis=1)
    reference_table = reference_table.drop_duplicates(subset='PEDIGREE_NAME', keep='last') # there are duplicates

    ############### Add reference table to diallel table
    
    print('\n========= Adding Reference Table to Diallel Table =========')
    start = time.time()
    print(f'{"Master Shape BEFORE adding reference data:":20}{diallel_table.shape}')
    diallel_table = utils.add_reference_dat(diallel_table, reference_table)
    print(f'{"Master Shape AFTER adding reference data:":20}{diallel_table.shape}')
          
    diallel_table['clean_origin'] = diallel_table['P1_M.GERMPLASM.ORIGIN'] + '__' + diallel_table['P2_M.GERMPLASM.ORIGIN']
    diallel_table['clean_pedigree'] = diallel_table['P1_PEDIGREE'] + '__' + diallel_table['P2_PEDIGREE']
          
          
    ######### Remove combinations if parental origin = parental origin
    print('\n ======== Eliminate combinations with the same P1/P2 origin ======')
          
    print(f'{"Master Shape BEFORE eliminating the same origin:":20}{diallel_table.shape}')
    dialle_table = rf.remove_same_origin(diallel_table)
    print(f'{"Master Shape AFTER eliminating the same origin:":20}{diallel_table.shape}')
          
    ############### Configure file ###################
    # Given configuration file to create subset diallel table for origin selection, hybrid prediction, and Pre-plant HMU
          
    
    if args.config is not None:
          print('\n========= Filtering Diallel Table Based On A Given Configuration File =========')
          diallel_table = utils.create_column_combo(diallel_table, DEFAULTS)
          DataObject = collections.namedtuple('DataObject', 'data')
          in_source = {}
          in_source['diallel'] = DataObject(data = diallel_table)
          in_source_filtered = rf.apply_rules(in_source, DEFAULTS['Filters'])
          
          ####################### Write results to a csv file
          print('\n========= Writing Results to a Local folder =========')
          final_df = in_source_filtered['diallel'].data
        
          final_df['existing_status'] = [True if (x in parentalTable['Pedigree_1']) or (x in parentalTable['Pedigree_2']) else False for x in final_df['clean_pedigree']]
          final_df.to_csv(args.local_path)
    else: 
          print('\n========= Writing Results to a Local folder =========')
          diallel_table['existing_status'] = [True if (x in parentalTable['Pedigree_1']) or (x in parentalTable['Pedigree_2']) else False for x in diallel_table['clean_pedigree']]
          diallel_table.to_csv(args.local_path)
    

        

    ###################### Upload file to S3
    
    if args.s3_path is not None:
        print('\n========= Writing Results to S3 =========')
        utils.upload_file(file_name = args.local_path, bucket= bucket, object_name = args.s3_path)
          
if __name__ == "__main__":
    args = parse_argument()
    main(args)
