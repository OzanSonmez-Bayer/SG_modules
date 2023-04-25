# Functions to filter diallel table: Xing Weng and Yujing Cao

import operator
import pandas as pd
import numpy as np
import sys
import collections
 
op_map = {
'greater_than': operator.gt,
'less_than': operator.lt,
'equal_to': operator.eq,
'not_equal_to': operator.ne,
'range': 'between',
'in': 'contains'
}
out_op_map = {
'greater_than': '>',
'less_than': '<',
'equal_to': '==',
'not_equal_to': '!=',
'range': 'between',
'in':'in'
}

def apply_rules(dataObject_1, configurations):
    dataObject = dataObject_1.copy()
    if 'rules' in configurations.keys():
        conf = configurations['configurations']
    else:
        conf = configurations
    unfound_table = [i for i in conf.keys() if i not in dataObject.keys()]
    found_table = [i for i in conf.keys() if i in dataObject.keys()]
    if len(unfound_table) > 0:
        print(f'{" ".join(unfound_table)} rules skipped, {" ".join(unfound_table)} table not found')
    for i in found_table:
        filtered_parent = dataObject[i].data
        rules = []
        if isinstance(conf[i]['rules'], list):
            rules = conf[i]['rules']
        elif isinstance(conf[i]['rules'], dict):
            rules.append(conf[i]['rules'])
        print(f'{i} table: \n\t\tshape: {filtered_parent.shape}')
        print('-----------------------------------')
        unfound_column = [i['trait'] for i in rules if i['trait'] not in filtered_parent.columns]
        assert len(unfound_column) == 0, f'Column: {unfound_column} not found in {i} table'
        for j in rules:
            # operation for in_range
            assert j['operator'] in op_map.keys(), f'"{j["operator"]}" operation not found in op_map'
            if j['operator'] == 'range':
                values = [i.strip() for i in j['value'].split(',')]
                try:
                    values = [float(i) for i in values]
                    col_data = filtered_parent[j['trait']].astype('float')
                    rules_bool_result = col_data.between(values[0],values[1])
                    filtered_parent = filtered_parent[rules_bool_result]
                except:
                    e = sys.exc_info()[0]
                    print( "<p>Error: %s</p>" % e )
            elif j['operator'] == 'in':
                values = [i.strip() for i in j['value'].split(',')]
                col_data = filtered_parent[j['trait']]
                rules_bool_result = np.any([col_data.str.contains(i, case=False, na=False) for i in values], axis=0)
                filtered_parent = filtered_parent[rules_bool_result]
            elif j['operator'] == 'greater_than' or j['operator'] == 'less_than':
                col_data = filtered_parent[j['trait']].astype('float')
                rules_bool_result = op_map[j['operator']](col_data, float(j['value']))
                filtered_parent = filtered_parent[rules_bool_result]
            elif j['operator'] == 'equal_to' or j['operator'] == 'not_equal_to':
                col_data = filtered_parent[j['trait']]
                rules_bool_result = op_map[j['operator']](col_data, j['value'])
                filtered_parent = filtered_parent[rules_bool_result]
            print('{trait} {operation} {value} \n\t\tshape: {shape:}'.format(
                trait=j['trait'], operation=out_op_map[j['operator']],
                value=j['value'], shape=str(filtered_parent.shape))
            )
            filtered_parent.index = range(filtered_parent.shape[0])
        dataObject[i] = dataObject[i]._replace(data = filtered_parent)
    return dataObject
              
              
def filter_by_HETGP(master_file, het):
    """
    master_file: master file
    het: A list. Heterotic group to keep
    """
    dat = master_file[master_file['HETGP_GCV'].isin(het)]
    
    return dat
              

def filter_on_het_combination(diallel_master, het_comb):
    """
    diallel_master: diallel master table with mid parent value and reference data
    het_comb: hetertoc group combinations of interest
    """
    
    diallel_master['het_comb'] = diallel_master['P1_HETGP_GCV'] + '_' + diallel_master['P2_HETGP_GCV']
    diallel_master = diallel_master[diallel_master['het_comb'].isin(het_comb)]
    diallel_master = diallel_master.drop(['het_comb'],axis=1)
    
    return diallel_master


def filter_on_inbrdstg(diallel_master, inbred_stage_comb):
    """
    diallel_master: diallel master table with mid parent value and reference data
    inbred_stage_comb: inbred stage combination of interests
    """
    
    diallel_master['inbred_stage_comb'] = diallel_master['P1_SC_INBRDSTG_GCV'] + '_' + diallel_master['P2_SC_INBRDSTG_GCV']
    diallel_master = diallel_master[diallel_master['inbred_stage_comb'].isin(inbred_stage_comb)]
    diallel_master = diallel_master.drop(['inbred_stage_comb'],axis=1)
    
    return diallel_master
              

def remove_same_origin(diallel_master): 
    """
    diallel_master: diallel master table with mid parent value and reference data
    """
    diallel_master_1 = diallel_master[diallel_master['P1_M.GERMPLASM.ORIGIN']!=diallel_master['P2_M.GERMPLASM.ORIGIN']]
    
    return diallel_master_1
