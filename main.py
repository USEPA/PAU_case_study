# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Importing libraries
import argparse
import os
import xlrd
import re
from xlutils.copy import copy
from bayesian_network.bayesian_network import *
from fuzzy_analytical_hierarchy_process.fuzzy_inference import *
from chemical_flow_analysis.chemical_tracking import *
from ph import corrosiveness_score
import pandas as pd
import json

if __name__ == '__main__':

    parser = argparse.ArgumentParser(argument_default = argparse.SUPPRESS)

    parser.add_argument('-CAS', nargs = '+',
                        help = 'Enter the unhyphenated CAS Number(s) of the chemical(s) you want to analyze.',
                        type = str,
                        required = False,
                        default = None)

    parser.add_argument('-Y', '--Year', nargs = '+',
                        help = 'What TRI year do you want to retrieve?.',
                        type = str,
                        required = False,
                        default = None)

    args = parser.parse_args()

    CAS = args.CAS
    Years = args.Year
    dir_path = os.path.dirname(os.path.realpath(__file__)) # Current directory
    # Bayesian Network
    df_PAU, CAS_for_search, Options = building_bayesian_network_db(CAS, Years, dir_path)

    # Building models
    for chem_1, chem_2 in CAS_for_search.items():
        building_bayesian_network_model(dir_path, df_PAU, chem_1, chem_2)

    # Filling excel file which will be used to enter the information
    read_book = xlrd.open_workbook(dir_path + '/Inputs.xls', formatting_info = True)
    write_book = copy(read_book)
    write_specifications = write_book.get_sheet(3)
    initial_pos = {'Chemical price': 29, 'Waste flow': 39, 'PACE': 49, 'PAOC': 59}
    idx = 0
    for chem_1, chem_2 in CAS_for_search.items():
        write_specifications.write(28, idx + 1, chem_1)
        for key, values in Options.items():
            df_PAU_chem = df_PAU[df_PAU['CAS NUMBER'] == chem_2]
            states = df_PAU_chem[key].unique().tolist()
            list_position = list()
            pos = initial_pos[key]
            for val in values:
                number = re.search(r'(^\d{1,2}):*', val).group(1)
                if number in states:
                    write_specifications.write(pos, idx + 1, val)
                    pos = pos + 1
        idx = idx + 1
    write_book.save(dir_path + '/Inputs.xls')
    print('-'*120)
    print('Fill out the required information in the "Input.xls". Check the options in the sheet called "Specifications"')
    print()
    input('When you finish filling out the sheet "Input", please press enter')
    print('-'*120)
    print()
    # Calling the inputs
    df_inputs = pd.read_excel(dir_path + '/Inputs.xls',
                                sheet_name = 'Input',
                                header = None,
                                skiprows=[0, 1, 2, 15, 19, 20])
    df_inputs[df_inputs.notnull()] = df_inputs[df_inputs.notnull()].astype(str)
    df_inputs.set_index(0, inplace = True)
    df_inputs = df_inputs.T
    df_inputs.rename_axis(None, axis=1, inplace=True)
    df_inputs.reset_index(inplace = True, drop = True)
    Columns = set(df_inputs.columns.tolist()) - set(['CAS NUMBER', 'Stream #',
                                                     'Chemical flow [kg/yr]',
                                                     'Efficiency [%]',
                                                     'Concentration [%wt/wt]'])
    streams = df_inputs['Stream #'].unique().tolist()
    concerning_chemical_in_stream = {stream: df_inputs.loc[df_inputs['Stream #'] == stream, 'CAS NUMBER'].tolist() for stream in streams}
    del streams
    # Assesing the Bayesian Network under the input evidence
    for stream, chemicals in concerning_chemical_in_stream.items():
        for chem_1 in chemicals:
            chem_2 = CAS_for_search[chem_1]
            Input_dictionary_joint = dict()
            Input_dictionary_marginal = dict()
            for col in Columns:
                df_input_chem = df_inputs.loc[(df_inputs['CAS NUMBER'] == chem_1) & (df_inputs['Stream #'] == stream)]
                val = df_input_chem[col].values[0]
                Input_dictionary_joint.update({col: val})
                Input_dictionary_marginal.update({col: val})
            with open(f'{dir_path}/bayesian_network/models/BN_for_{chem_1}.json', 'r') as openfile:
                json_object = json.load(openfile)
                PAU_model = from_json(json_object)
            calculating_joint_probabilities(Input_dictionary_joint, chem_1, chem_2, stream, PAU_model, dir_path, df_PAU)
            calculating_marginal_probabilities(Input_dictionary_marginal, PAU_model, dir_path, chem_1, chem_2, stream, df_PAU)
    # Fuzzy analysis
    Prob_PAU_chemical_and_stream = pd.DataFrame()
    for stream, chemicals in concerning_chemical_in_stream.items():
        for chemical in chemicals:
            Path = f'/bayesian_network/probabilities/marginal/Marginal_probabilities_based_on_BN_for_{chemical}_in_stream_{stream}.csv'
            Prob_chem = pd.read_csv(dir_path + Path)
            Prob_chem = Prob_chem[Prob_chem['PAU-probability'] != 0.0]
            if Prob_chem.empty:
                print(f'Based on the information, there is not chance to isolate the chemical {chemical} in stream {stream} under the input conditions')
            else:
                Prob_chem['Stream'] = stream
                Prob_chem['Chemical'] = chemical
                Prob_PAU_chemical_and_stream = pd.concat([Prob_PAU_chemical_and_stream, Prob_chem],
                                                        axis = 0, sort = False, ignore_index = True)
    del Prob_chem

    ## Selecting PAUs
    Prob_PAU_chemical_and_stream = Prob_PAU_chemical_and_stream\
                                    .groupby('Stream', as_index = False)\
                                    .apply(lambda x: pairwise_comparison(x, objective = 'pau'))
    Prob_PAU_chemical_and_stream = Prob_PAU_chemical_and_stream[Prob_PAU_chemical_and_stream['Selected'] == 'Yes']
    Prob_PAU_chemical_and_stream.drop(columns = ['Selected'], inplace = True)
    Prob_PAU_chemical_and_stream.reset_index(inplace = True, drop = True)

    ## Sequence for PAUs
    df_input_for_chems =  pd.read_excel(dir_path + '/Inputs.xls',
                                sheet_name = 'Chemical',
                                header = None,
                                skiprows=[0, 11])
    df_input_for_chems.set_index(0, inplace = True)
    df_input_for_chems = df_input_for_chems.T
    df_input_for_chems.rename_axis(None, axis=1, inplace=True)
    df_input_for_chems.reset_index(inplace = True, drop = True)
    df_input_for_chems[['CAS NUMBER', 'Flammability', 'Instability']] =\
        df_input_for_chems[['CAS NUMBER', 'Flammability', 'Instability']].applymap(lambda x: str(int(x)))
    df_inputs = pd.merge(df_inputs, df_input_for_chems, on='CAS NUMBER', how='inner')
    del df_input_for_chems
    df_inputs['Corrosiveness'] = df_inputs.apply(lambda x: corrosiveness_score(
                                                x['pKa Acidic Apparent'],
                                                x['pKa Basic Apparent'],
                                                float(x['Concentration [%wt/wt]']),
                                                x['Density [g/cm3]'],
                                                x['Molecular Mass'])
                                                , axis=1)
    df_inputs = df_inputs[['Stream #', 'Chemical flow [kg/yr]', 'CAS NUMBER',
                           'Flammability', 'Instability', 'Corrosiveness',
                           'Efficiency [%]', 'As a byproduct',
                           'As a manufactured impurity', 'As a process impurity',
                           'Type of waste', 'Melting Point [°C]',
                           'Boiling Point [°C]', 'Water Solubility [mg/L]']]
    df_inputs.rename(columns = {'Stream #': 'Stream',
                                'Chemical flow [kg/yr]': 'Chemical flow',
                                'CAS NUMBER': 'Chemical',
                                'Water Solubility [mg/L]': 'Solubility',
                                'Boiling Point [°C]': 'Tb',
                                'Melting Point [°C]': 'Tf'}, inplace = True)
    df_inputs['Chemical flow'] = df_inputs['Chemical flow'].astype('float')
    df_inputs['Tf'] = df_inputs['Tf'].astype('float')
    df_inputs['Tb'] = df_inputs['Tb'].astype('float')
    df_inputs['Solubility'] = df_inputs['Solubility'].astype('float')
    df_inputs['Efficiency [%]'] = df_inputs['Efficiency [%]'].astype('float')
    df_inputs['Flammability'] = df_inputs['Flammability'].astype('int')
    df_inputs['Instability'] = df_inputs['Instability'].astype('int')
    df_inputs['Corrosiveness'] = df_inputs['Corrosiveness'].astype('int')
    Prob_PAU_chemical_and_stream = pd.merge(df_inputs, Prob_PAU_chemical_and_stream,
                                            on = ['Chemical', 'Stream'],
                                            how = 'inner')
    df_TRI_methods = pd.read_csv(dir_path + '/Methods_TRI.csv',
                                usecols = ['Code 2004 and prior', 'Objective',
                                           'Method 2004 and prior'])
    df_TRI_methods.rename(columns = {'Code 2004 and prior': 'PAU',
                                    'Method 2004 and prior': 'PAU name'}, inplace = True)
    Prob_PAU_chemical_and_stream = pd.merge(df_TRI_methods, Prob_PAU_chemical_and_stream,
                                            on = ['PAU'], how = 'inner')
    Result = Prob_PAU_chemical_and_stream.copy()
    Prob_PAU_chemical_and_stream = pd.DataFrame()
    for stream in Result['Stream'].unique():
        df = Result.loc[Result['Stream'] == stream]
        if df.shape[0] == 1:
            df['Position'] = 1
            Prob_PAU_chemical_and_stream = pd.concat([Prob_PAU_chemical_and_stream, df], axis = 0,
                                                    sort = False,
                                                    ignore_index = True)
        else:
            df = pairwise_comparison(df, objective = 'seq')
            Prob_PAU_chemical_and_stream = pd.concat([Prob_PAU_chemical_and_stream, df], axis = 0,
                                                    sort = False,
                                                    ignore_index = True)
    Prob_PAU_chemical_and_stream.sort_values(by = ['Stream', 'Position'],
                                             ascending = [True, True],
                                             inplace = True)
    Cols_for_saving = ['Stream', 'PAU', 'PAU name', 'Objective', 'Position', 'Chemical']
    Prob_PAU_chemical_and_stream[Cols_for_saving].to_csv(dir_path + '/fuzzy_analytical_hierarchy_process/PAU_selection_and_position_under_FAHP.csv',
                                                         sep = ',', index = False)
    # Chemical flow analysis
    Chemical_tracking = pd.DataFrame()
    PAUs = list()
    for stream in concerning_chemical_in_stream.keys():
        df_for_stream = Prob_PAU_chemical_and_stream.loc[Prob_PAU_chemical_and_stream['Stream'] == stream]
        paus_stream = '_'.join(set(df_for_stream['PAU'].tolist()))
        if paus_stream not in PAUs:
            drawing = True
            PAUs.append(paus_stream)
        else:
            drawing = False
        Chemical_tracking_aux = picture(df_for_stream, stream, dir_path, drawing, paus_stream)
        Chemical_tracking = pd.concat([Chemical_tracking, Chemical_tracking_aux], axis = 0,
                    sort = False,
                    ignore_index = True)
    # Removing rows without values (non information)
    Chemical_tracking = Chemical_tracking.loc[pd.notnull(Chemical_tracking[['Mean quantity [kg/yr]', 'Std [kg/yr]']]).all(axis=1)]
    # Removing rows with values 0
    Chemical_tracking = Chemical_tracking.loc[(Chemical_tracking[['Mean quantity [kg/yr]', 'Std [kg/yr]']] != 0.0).all(axis=1)]
    Chemical_tracking['CV'] = Chemical_tracking['Std [kg/yr]']/Chemical_tracking['Mean quantity [kg/yr]']
    cols = ['Stream', 'PAU', 'PAU name', 'Objective', 'Position', 'Type of stream', 'Phase',
            'Chemical', 'Mean quantity [kg/yr]', 'Std [kg/yr]', 'CV', 'Number of data']
    Chemical_tracking = pd.merge(df_TRI_methods, Chemical_tracking,
                                            on = ['PAU'], how = 'right')
    Chemical_tracking = Chemical_tracking[cols]
    Chemical_tracking.sort_values(by = ['Stream', 'Position', 'Chemical'],
                                        ascending = [True, True, True],
                                        inplace = True)
    Chemical_tracking.to_csv(dir_path + '/chemical_flow_analysis/Chemical_flow_tracking.csv',
                                        sep = ',', index = False)
