# Data headers are consistent with H2H data
# Yujing Cao

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
import time
import os
import sys
import argparse
import utils

def read_s3_data(filename):
    with open('/mnt/vaultCredentials.json') as f:
        vaultCredentials = json.loads(f.read())
    
    client_vault = hvac.Client("https://vault.agro.services")
    client_vault.auth_approle(vaultCredentials['AppRoleCredsJson'][0]['role_id'],vaultCredentials['AppRoleCredsJson'][0]['secret_id'] )
    d = client_vault.read(vaultCredentials['VaultSecretPath'])
    appdata = d["data"]
    
    # how to store credentials
    bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace"
    client = boto3.client(
        's3',
        aws_access_key_id=appdata['AWS_ACCESS_KEY_ID'],
        aws_secret_access_key=appdata['AWS_SECRET_ACCESS_KEY']
    )

    obj = client.get_object(Bucket= bucket, Key=filename)
    df = pd.read_csv(obj['Body']) 
    df = df.rename(columns = {'Unnamed: 0':'PEDIGREE'})
    df.to_csv(os.path.join(os.getcwd(), 'inputfile.csv'), index=False)
    
    return df

def read_as_dask(filename):
    
    read_s3_data(filename)
    
    dat = dd.read_csv(os.path.join(os.getcwd(), 'inputfile.csv'))
    
#     cols_1 = [col for col in dat.columns if '_GCV' in col]
#     cols_2 = [col for col in dat.columns if '_TLM' in col]
#     cols = cols_1 + cols_2 + ['M.GERMPLASM.CROSSNAME', 'M.GERMPLASM.X_ID']
    
# #     del dat
    
#     dtype = {cols_name: 'object' for cols_name in cols} 
    dat_dask = dd.read_csv(os.path.join(os.getcwd(), 'inputfile.csv'), 
#                       dtype = dtype,
                    blocksize = 1024*1024)
    return dat_dask


def parent_mean(dat, traitname):
#     dat = dat.rename(columns = {'M.GERMPLASM.PEDIGREE': 'PEDIGREE'})
#     trtcol = traitname+'_predicted.value'
    trtcol = 'predicted.value_' + traitname
    dat_trt = dat[['PEDIGREE_NAME', trtcol]]
    dat_trt.dropna(inplace=True)
    dat_trt.reset_index(drop=True, inplace=True)
    ID_combine = np.array(list(itertools.combinations_with_replacement(dat_trt.index, 2)))
    p1_labels = dat_trt[['PEDIGREE_NAME']].iloc[ID_combine[:,0]].values
    p2_labels = dat_trt[['PEDIGREE_NAME']].iloc[ID_combine[:,1]].values
    trtcol_vals = dat_trt[[trtcol]].values
    trt_col_avgs = (trtcol_vals[ID_combine[:,0]] + trtcol_vals[ID_combine[:,1]]) / 2
    dat_merge = pd.DataFrame(
        {
            "P1_PEDIGREE": p1_labels.flatten(),
            "P2_PEDIGREE": p2_labels.flatten(),
            'P1_' + traitname: trtcol_vals[ID_combine[:,0]].flatten(),
            'P2_' + traitname: trtcol_vals[ID_combine[:,1]].flatten(),
            f"Parent_Mean_{traitname}": trt_col_avgs.flatten()
        }
    )
    return dat_merge

def mid_reliability(dat, traitname):
    trtcol = traitname + '_reliability'
    xcol = traitname + '_n'
    dat_trt = dat[['PEDIGREE_NAME', trtcol, xcol]]
#     print(dat_trt)
    dat_trt.dropna(inplace=True)
    dat_trt.reset_index(drop=True, inplace=True)
#     print(dat_trt)
    ID_combine = np.array(list(itertools.combinations_with_replacement(dat_trt.index, 2)))
#     print(ID_combine.shape)
#     print(trtcol)
    trtcol_vals = dat_trt[[trtcol]].values
    p1_labels = dat_trt[['PEDIGREE_NAME']].iloc[ID_combine[:,0]].values
    p2_labels = dat_trt[['PEDIGREE_NAME']].iloc[ID_combine[:,1]].values
    trtcol_n = dat_trt[[xcol]].values
    trt_col_avgs = (trtcol_vals[ID_combine[:,0]] * trtcol_n[ID_combine[:, 0]] + trtcol_vals[ID_combine[:,1]] * trtcol_n[ID_combine[:,1]]) / \
    (trtcol_n[ID_combine[:, 0]] + trtcol_n[ID_combine[:, 1]])
    dat_merge = pd.DataFrame(
        {
            "P1_PEDIGREE": p1_labels.flatten(),
            "P2_PEDIGREE": p2_labels.flatten(),
            'P1_' + trtcol: trtcol_vals[ID_combine[:,0]].flatten(),
            'P2_' + trtcol: trtcol_vals[ID_combine[:,1]].flatten(),
            'P1_' + xcol: trtcol_n[ID_combine[:,0]].flatten(),
            'P2_' + xcol: trtcol_n[ID_combine[:,1]].flatten(),
            f"Mid_Reliability_{traitname}": trt_col_avgs.flatten()
        }
    )
    return dat_merge

def gather_mid_value(dat, trtlist):
    results = []
    for trait in trtlist:
        mid_blup = delayed()(parent_mean)(dat, trait)
#         mid_re = delayed()(mid_reliability)(dat, trait)
#         final=None
#         final = delayed()(pd.merge)(mid_blup, mid_re, on=['P1_PEDIGREE', 'P2_PEDIGREE'], how='inner')
#         final = delayed()(final.drop_duplicates)()
        final = delayed()(mid_blup.drop_duplicates)()
        results.append(final)
    df_merged = reduce(lambda  left,right: delayed()(pd.merge)(left, right, on=['P1_PEDIGREE', 'P2_PEDIGREE'], how='outer'), results)
    return df_merged.compute()

