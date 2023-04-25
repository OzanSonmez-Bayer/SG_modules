# Generic help functions 
import pandas as pd
import boto3
import numpy as np
import json
import requests
import hvac
import dask
from dask.distributed import Client
from dask import delayed
import dask.dataframe as dd
from functools import reduce
import itertools
import os
import sys
import argparse

def read_s3_data(filename, bucket):
    with open('/mnt/vaultCredentials.json') as f:
        vaultCredentials = json.loads(f.read())
    
    client_vault = hvac.Client("https://vault.agro.services")
    client_vault.auth_approle(vaultCredentials['AppRoleCredsJson'][0]['role_id'],vaultCredentials['AppRoleCredsJson'][0]['secret_id'] )
    d = client_vault.read(vaultCredentials['VaultSecretPath'])
    appdata = d["data"]
    
    # how to store credentials
#     bucket = "veg-apd-sdi-predictiveanalytcs-prod-alwaysonblup"
    bucekt = bucket
    client = boto3.client(
        's3',
        aws_access_key_id=appdata['AWS_ACCESS_KEY_ID'],
        aws_secret_access_key=appdata['AWS_SECRET_ACCESS_KEY']
    )

    obj = client.get_object(Bucket= bucket, Key=filename)
    df = pd.read_csv(obj['Body']) 
    
    if 'Unnamed: 0' in df:
        df = df.rename(columns = {'Unnamed: 0':'PEDIGREE'})
    if 'Pedigree' in df:
        df = df.rename(columns = {'Pedigree':'PEDIGREE'})
    df.to_csv(os.path.join(os.getcwd(), 'inputfile.csv'), index=False)
    
    return df

def read_as_dask(filename, bucket):
    
    read_s3_data(filename, bucket)
    
    dat = dd.read_csv(os.path.join(os.getcwd(), 'inputfile.csv'))
    
    cols_1 = [col for col in dat.columns if '_GCV' in col]
    cols_2 = [col for col in dat.columns if '_TLM' in col]
    cols = cols_1 + cols_2 + ['M.GERMPLASM.CROSSNAME', 'M.GERMPLASM.X_ID']
    
#     del dat
    
    dtype = {cols_name: 'object' for cols_name in cols} 
    dat_dask = dd.read_csv(os.path.join(os.getcwd(), 'inputfile.csv'), 
                      dtype = dtype,
                    blocksize = 1024*1024)
    return dat_dask

def read_s3_data_1(filename,bucket, output=False):
    
    """
    read_s3_data_1 is used when the model is running on aws
    """
    
    bucket = "veg-apd-sdi-predictiveanalytcs-prod-alwaysonblup"
    client = boto3.client('s3')

    obj = client.get_object(Bucket= bucket, Key=filename)
    df = pd.read_csv(obj['Body']) 
    if 'Unnamed: 0' in df:
        df = df.rename(columns = {'Unnamed: 0':'PEDIGREE'})
    if 'Pedigree' in df:
        df = df.rename(columns = {'Pedigree':'PEDIGREE'})
    df.to_csv(os.path.join(os.getcwd(),'inputfile.csv'), index=False)
    
    if output==True: 
        return df

def read_as_dask_1(filename,bucket, output=False):
    
    """
    read_as_dask_1 is used when the model is running on aws
    """
    
    read_s3_data_1(filename, output)
    
    dat = dd.read_csv(os.path.join(os.getcwd(),'inputfile.csv'))
    
    cols_1 = [col for col in dat.columns if '_GCV' in col]
    cols_2 = [col for col in dat.columns if '_TLM' in col]
    cols = cols_1 + cols_2 + ['M.GERMPLASM.CROSSNAME', 'M.GERMPLASM.X_ID']
    
    del dat
    
    dtype = {cols_name: 'object' for cols_name in cols} 
    dat_dask = dd.read_csv(os.path.join(os.getcwd(),'inputfile.csv'), 
                      dtype = dtype,
                    blocksize = 1024*1024)
    return dat_dask


def read_geno_s3_data(filename):
    """
    filename: path to genotypic data in s3
    """
    with open('/mnt/vaultCredentials.json') as f:
        vaultCredentials = json.loads(f.read())
    
    client_vault = hvac.Client("https://vault.agro.services")
    client_vault.auth_approle(vaultCredentials['AppRoleCredsJson'][0]['role_id'],vaultCredentials['AppRoleCredsJson'][0]['secret_id'] )
    d = client_vault.read(vaultCredentials['VaultSecretPath'])
    appdata = d["data"]
    
    # how to store credentials
    bucket = "veg-apd-sdi-predictiveanalytics-prod-geno-data"
#     bucket = bucket
    client = boto3.client(
        's3',
        aws_access_key_id=appdata['AWS_ACCESS_KEY_ID'],
        aws_secret_access_key=appdata['AWS_SECRET_ACCESS_KEY']
    )

    obj = client.get_object(Bucket= bucket, Key=filename)
    df = pd.read_csv(obj['Body'])     
    return df


def read_s3_file_generic(filename, bucket):
    with open('/mnt/vaultCredentials.json') as f:
        vaultCredentials = json.loads(f.read())
    
    client_vault = hvac.Client("https://vault.agro.services")
    client_vault.auth_approle(vaultCredentials['AppRoleCredsJson'][0]['role_id'],vaultCredentials['AppRoleCredsJson'][0]['secret_id'] )
    d = client_vault.read(vaultCredentials['VaultSecretPath'])
    appdata = d["data"]
    
    # how to store credentials
#     bucket = "veg-apd-sdi-predictiveanalytcs-prod-alwaysonblup"
    bucekt = bucket
    client = boto3.client(
        's3',
        aws_access_key_id=appdata['AWS_ACCESS_KEY_ID'],
        aws_secret_access_key=appdata['AWS_SECRET_ACCESS_KEY']
    )

    obj = client.get_object(Bucket= bucket, Key=filename)
    df = pd.read_csv(obj['Body']) 
    return df


def upload_file(file_name, bucket, object_name=None):
    """Upload a file to an S3 bucket

    :param file_name: File to upload
    :param bucket: Bucket to upload to
    :param object_name: S3 object name. If not specified then file_name is used
    :return: True if file was uploaded, else False
    """

    # If S3 object_name was not specified, use file_name
    if object_name is None:
        object_name = file_name
    
    with open('/mnt/vaultCredentials.json') as f:
        vaultCredentials = json.loads(f.read())
    client_vault = hvac.Client("https://vault.agro.services")
    client_vault.auth_approle(vaultCredentials['AppRoleCredsJson'][0]['role_id'],vaultCredentials['AppRoleCredsJson'][0]['secret_id'] )
    d = client_vault.read(vaultCredentials['VaultSecretPath'])
    appdata = d["data"]


#     bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace"
    client = boto3.client(
        's3',
        aws_access_key_id=appdata['AWS_ACCESS_KEY_ID'],
        aws_secret_access_key=appdata['AWS_SECRET_ACCESS_KEY']
    )

    # Upload the file
    try:
        response = client.upload_file(file_name, bucket, object_name)
    except ClientError as e:
        logging.error(e)
        return False
    return True



def add_reference_dat(diallel_table, reference_dat):
    """
    diallel_table: panda data frame. diallel table
    reference_dat:reference data from filtered master file 
    """
    col_names = reference_dat.columns
    rename_dict_P1 = {cc : 'P1_' + cc for cc in col_names}
    rename_dict_P2 = {cc : 'P2_' + cc for cc in col_names}
    
    diallel_table_dd = dd.from_pandas(diallel_table, npartitions=10)
    join_table = delayed()(pd.merge)(diallel_table_dd, dd.from_pandas(reference_dat.rename(columns=rename_dict_P1), npartitions=3), 
                                     left_on='P1_PEDIGREE', right_on='P1_PEDIGREE_NAME', how='left')
    join_table = delayed()(pd.merge)(join_table, dd.from_pandas(reference_dat.rename(columns=rename_dict_P2), npartitions=3),
                                    left_on='P2_PEDIGREE', right_on='P2_PEDIGREE_NAME', how='left')
    
#     join_table = pd.merge(left=diallel_table, right=reference_dat.rename(columns=rename_dict_P1), left_on='P1_PEDIGREE', right_on='P1_PEDIGREE_NAME', how='left')
#     join_table = pd.merge(left=join_table, right=reference_dat.rename(columns=rename_dict_P2), left_on='P2_PEDIGREE', right_on='P2_PEDIGREE_NAME', how='left')
    
    return join_table.compute()

def split_diallel_table(diallel_table): 
    diallel_table['Goal'] = 'Origin Selection'
    for row_indx in range(1, diallel_table.shape[0]):
        if diallel_table.loc[row_indx, 'P1_HETGP_GCV'] != diallel_table.loc[row_indx, 'P2_HETGP_GCV']:
            diallel_table.loc[row_indx, 'Goal'] = 'Hybrid Prediction'
            
    devX_df = diallel_table[diallel_table['Goal'] == 'Origin Selection']
    hmu_df = diallel_table[diallel_table['Goal'] == 'Hybrid Prediction']
            
    return devX_df, hmu_df

# Find the columns of interest and create a combo column
def create_column_combo(dat, defaults): 
    """
    input:
        dat: unfiltered diallel table
        defaults: configuration file with filtering criteria
    output: diallel table with new columns
    """
    cols = [defaults['Filters']['diallel']['rules'][i]['trait'] for i in range(len(defaults['Filters']['diallel']['rules']))]
    
    for cc in cols: 
        print(cc)
        p1_cc = 'P1_' + cc
        p2_cc = 'P2_' + cc
        if (p1_cc in dat.columns) & (p2_cc in dat.columns): 
            dat[cc] = dat[p1_cc] + '_' + dat[p2_cc]
        else: 
            print("Column doesn't exsit, please check column name")
    return dat
    
