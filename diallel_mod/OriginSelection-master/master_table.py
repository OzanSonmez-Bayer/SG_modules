# Functions to generate master table written by Xing Weng
# The original file is https://domino.science-at-scale.io/u/VPS/Always-On_BLUP/view/test/master_file_test.ipynb

import pandas as pd
import numpy as np

def get_vault_token(vault_role_id, vault_secret_id):
    url = 'https://vault.agro.services/v1/auth/approle/login'
    vault_cred = {'role_id': vault_role_id, 'secret_id': vault_secret_id}
    response = requests.post(url, json=vault_cred)
    assert response.status_code == 200, response.json()['errors']
    return response.json()['auth']['client_token']

def get_aws_cred(vault_path, token):
    url = f'https://vault.agro.services/v1/secret/{vault_path}'
    headers = {'X-Vault-Token':token}
    response = requests.get(url, headers=headers)
    return response.json()['data']

def combine_blup_id(blup_wide, id_table):
    print(f'{"Combine Id_Table":_^30}')
    print(f'{"Initial Size":20}{blup_wide.shape} {id_table.shape}')
    id_column = ['M.GERMPLASM.X_ID','M.GERMPLASM.PEDIGREE','M.GERMPLASM.CROSSNAME','M.GERMPLASM.ORIGIN',
                     'M.LINETYPE', 'M.GERMPLASM.LINECODE', 'M.GERMPLASM.CODEYEAR']
    id_table = id_table[id_column].drop_duplicates(keep='first')
    blup_id = blup_wide.merge(id_table, how='left', 
                              left_on='PEDIGREE_NAME', right_on=['M.GERMPLASM.PEDIGREE'])
    blup_id.drop(blup_id.columns[blup_id.isnull().values.all(axis=0)], axis=1, inplace=True)
    print(f'Removed {blup_id["M.GERMPLASM.X_ID"].isnull().sum()} null germplasm')
    blup_id.dropna(subset=['M.GERMPLASM.X_ID'], inplace=True)
    print(f'{"Master Shape:":20}{blup_id.shape}')

    return blup_id
          
          

def combine_blup_obs(master_file, gpc_df):
    print(f'{"Combine GPC_Table":_^30}')
    print(f'{"Initial Size":20}{master_file.shape} {gpc_df.shape}')
#     filtered_gpc = pd.DataFrame(master_file['M.GERMPLASM.X_ID']).merge(gpc_df, how='left',
#                     left_on=[('M.GERMPLASM.X_ID')], right_on=[('O.GCV.X_GPID')])[['O.OBSRVTN.OBSRVTN_REF_CD', 'O.GCV.X_VALUE', 'O.GCV.X_GPID']]
#     obs_ref_table = filtered_gpc.pivot_table(
#         values='O.GCV.X_VALUE', index='O.GCV.X_GPID', columns='O.OBSRVTN.OBSRVTN_REF_CD',
#         aggfunc=lambda x: '|'.join(str(v) for v in x))
#     trait_columns = obs_ref_table.columns
#     obs_ref_table.rename(columns = {f'{x}':f'{x}_GCV' for x in trait_columns}, inplace=True)
#     blup_obs = master_file.merge(obs_ref_table, how='left',
#                              left_on=['M.GERMPLASM.X_ID'], right_index=True)
#     blup_obs.drop(blup_obs.columns[blup_obs.isnull().values.all(axis=0)], axis=1, inplace=True)
#     print(f'{"Master Shape:":20}{blup_obs.shape}')
#     return blup_obs
          
    if 'Unnamed: 0' in gpc_df.columns:
        gpc_df = gpc_df.drop(['Unnamed: 0'],axis=1)
    columns = list([gpc_df.columns[0]])+list(gpc_df.columns[4:])
    gpc_df = gpc_df[columns]
    gpc_df.rename(columns = {f'{x}':f'{x}_GCV' for x in gpc_df.columns[1:]}, inplace=True)
    blup_obs = master_file.merge(gpc_df, how='left',
                                 left_on='PEDIGREE_NAME', right_on=['Pedigree'])
    blup_obs.drop(blup_obs.columns[blup_obs.isnull().values.all(axis=0)], axis=1, inplace=True)
    print(f'{"Master Shape:":20}{blup_obs.shape}')
    return blup_obs

def combine_blup_markerData(master_file, markerData):
    print(f'{"Combine Marker_Table":_^30}')
    print(f'{"Initial Size":20}{master_file.shape} {markerData.shape}')
    traitLocus_table = pd.pivot_table(markerData,values='ALLELE_VALUE',
        index=['GERMPLASM_ID'], columns=['TRAITLOCUSNAME'], aggfunc=lambda x: '|'.join(str(v) for v in x))
    trait_columns = traitLocus_table.columns
    traitLocus_table.rename(columns = {f'{x}':f'{x}_TLM' for x in trait_columns}, inplace=True)
    combined_df = master_file.merge(traitLocus_table, how='left',
                                    left_on=['M.GERMPLASM.X_ID'], right_index=True)
    combined_df = combined_df.drop(combined_df.columns[combined_df.isnull().values.all(axis=0)], axis=1)
    print(f'{"Master Shape:":20}{combined_df.shape}')
    return combined_df

def combine_blup_parental(master_file, parentalTable):
    """
    If FemaleCount > 1 Then FemalePedigree = NA
    If MaleCount > 1 Then MalePedigree = NA
    if FemalePedigree == MalePedigree then MalePedigree = NA
    if FemalePedigree == NA and MalePedigree != NA then MalePedigree = NA
    After all that FemalePedigree = P1
    MalePedigree = P2
    """
    print(f'{"Combine Parental_Table":_^30}')
    print(f'{"Initial Size":20}{master_file.shape} {parentalTable.shape}')

    parentalTable[['FemaleCount','MaleCount']] = parentalTable[['FemaleCount','MaleCount']].astype('float32')
    parentalTable.loc[parentalTable['FemaleCount'] > 1, 'FemalePedigree'] = np.NaN
    parentalTable.loc[parentalTable['MaleCount'] > 1, 'MalePedigree'] = np.NaN
    parentalTable.loc[parentalTable['FemalePedigree']==parentalTable['MalePedigree'], 'MalePedigree'] = np.NaN
    parentalTable.loc[parentalTable['FemalePedigree'].isna() & parentalTable['MalePedigree'].notna(), 'MalePedigree'] = np.NaN
    temp_df = parentalTable[['MalePedigree','FemalePedigree']].rename(columns={'FemalePedigree':'P1','MalePedigree':'P2'})
    temp_df.index = parentalTable['ChildGermID']
    merge_df = master_file.merge(right=temp_df, how='left', left_on=['M.GERMPLASM.X_ID'],
                                 right_index=True)
    merge_df = merge_df.drop(merge_df.columns[merge_df.isnull().values.all(axis=0)], axis=1)
    print(f'{"Master Shape:":20}{merge_df.shape}')
    return merge_df
          
          
def combine_blup_stage(master_file, stageData):
    print(f'{"Combine Stage_Table":_^30}')
    print(f'{"Initial Size":20}{master_file.shape} {stageData.shape}')
    stageData.dropna(inplace=True)
    stageData['stage'] = stageData['stage']+'_hybrid_stage_count'
    stage_count = stageData.pivot_table(
        values='stage_count', index='parent_pedigree_name', columns='stage')
    merge_df = master_file.merge(right=stage_count, how='left', left_on='M.GERMPLASM.PEDIGREE', right_index=True)
    merge_df.index = master_file.index
    merge_df = merge_df.drop(merge_df.columns[merge_df.isnull().values.all(axis=0)], axis=1)
    print(f'{"Master Shape:":20}{merge_df.shape}')
    return merge_df
