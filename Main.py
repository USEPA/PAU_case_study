# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Importing libraries
import argparse
import os
import xlrd, xlwt
from xlutils.copy import copy
from Bayesian_Network.Bayesian_Network import *
from Fuzzy_Analytical_Hierarchy_Process.Fuzzy_Inference import *
import pandas as pd

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
    N_Bins = 10
    dir_path = os.path.dirname(os.path.realpath(__file__)) # Current directory


    # Bayesian Network
    df_PCU, CAS_for_search, Options = building_bayesian_network_db(CAS, Years, N_Bins, dir_path)
    # Filling excel file which will be used to enter the information
    read_book = xlrd.open_workbook(dir_path + '/Inputs.xls', formatting_info = True)
    write_book = copy(read_book)
    write_sheet1 = write_book.get_sheet(1)
    initial_pos = {'Chemical price': 22, 'Waste flow': 32, 'PACE': 42, 'PAOC': 52}
    for key, values in Options.items():
        list_position = list()
        pos = initial_pos[key]
        for val in values:
            write_sheet1.write(pos, 1, val)
            pos = pos + 1
    write_book.save(dir_path + '/Inputs.xls')
    print('-'*120)
    print('Fill out the required information in the "Input.xls". Check the options in the sheet called "Specifications"')
    print()
    input('When you finish filling out the sheet "Input", please press enter')
    print('-'*120)
    print()
    df_inputs = pd.read_excel(dir_path + '/Inputs.xls',
                                sheet_name = 'Input',
                                header = None)
    df_inputs[df_inputs.notnull()] = df_inputs[df_inputs.notnull()].astype(str)
    df_inputs.set_index(0, inplace = True)
    df_inputs = df_inputs.T
    df_inputs.reset_index(inplace = True, drop = True)
    del df_inputs.index.name
    Columns = set(df_inputs.columns.tolist()) - set(['CAS NUMBER', 'Stream #', 'Chemical flow [kg/yr]',
                                                    'Flammability', 'Instability', 'Corrosivity'])
    streams = df_inputs['Stream #'].unique().tolist()
    concerning_chemical_in_stream = {stream: df_inputs.loc[df_inputs['Stream #'] == stream, 'CAS NUMBER'].tolist() for stream in streams}
    del streams
    for stream, chemicals in concerning_chemical_in_stream.items():
        for chem_1 in chemicals:
            chem_2 = CAS_for_search[chem_1]
            Input_dictionary_joint = dict()
            Input_dictionary_marginal = dict()
            for col in Columns:
                df_input_chem = df_inputs.loc[(df_inputs['CAS NUMBER'] == chem_1) & (df_inputs['Stream #'] == stream)]
                val = df_input_chem[col].values[0]
                Input_dictionary_joint.update({col: val})
                if col != 'Type of waste management':
                    Input_dictionary_marginal.update({col: val})
            PCU_model = Building_bayesian_network_model(dir_path, df_PCU, chem_1, chem_2)
            Calculating_joint_probabilities(Input_dictionary_joint, chem_1, chem_2, stream, PCU_model, dir_path, df_PCU)
            Calculating_marginal_probabilities(Input_dictionary_marginal, PCU_model, dir_path, chem_1, chem_2, stream, df_PCU)
    # Fuzzy analysis
    Prob_PCU_chemical_and_stream = pd.DataFrame()
    for stream, chemicals in concerning_chemical_in_stream.items():
        for chemical in chemicals:
            Path = '/Bayesian_Network/Probabilities/Marginal/Marginal_probabilities_based_on_BN_for_{}_in_stream_{}.csv'.format(chemical, stream)
            Prob_chem = pd.read_csv(dir_path + Path)
            Prob_chem = Prob_chem[Prob_chem['PCU-probability'] != 0.0]
            if Prob_chem.empty:
                print('Based on the information, there is not chance to isolate the chemical under the input conditions')
            else:
                Prob_chem['Stream'] = stream
                Prob_chem['Chemical'] = chemical
                Prob_PCU_chemical_and_stream = pd.concat([Prob_PCU_chemical_and_stream, Prob_chem],
                                                        axis = 0, sort = False, ignore_index = True)
    del Prob_chem
    ## Selecting PCUs
    Prob_PCU_chemical_and_stream = Prob_PCU_chemical_and_stream\
                                   .groupby('Stream', as_index = False)\
                                   .apply(lambda x: pairwise_comparison(x, objective = 'pcu'))
    Prob_PCU_chemical_and_stream = Prob_PCU_chemical_and_stream[Prob_PCU_chemical_and_stream['Selected'] == 'Yes']
    Prob_PCU_chemical_and_stream.drop(columns = ['Selected'], inplace = True)
    Prob_PCU_chemical_and_stream.reset_index(inplace = True, drop = True)
    ## Sequence for PCUs
    df_inputs = df_inputs[['Stream #', 'Chemical flow [kg/yr]', 'CAS NUMBER',
                           'Flammability', 'Instability', 'Corrosivity']]
    df_inputs.rename(columns = {'Stream #': 'Stream',
                                'Chemical flow [kg/yr]': 'Chemical flow',
                                'CAS NUMBER': 'Chemical'}, inplace = True)
    df_inputs['Chemical flow'] = df_inputs['Chemical flow'].astype('float')
    df_inputs['Flammability'] = df_inputs['Flammability'].astype('int')
    df_inputs['Instability'] = df_inputs['Instability'].astype('int')
    df_inputs['Corrosivity'] = df_inputs['Corrosivity'].astype('int')
    Prob_PCU_chemical_and_stream = pd.merge(df_inputs, Prob_PCU_chemical_and_stream,
                                            on = ['Chemical', 'Stream'],
                                            how = 'inner')
    df_TRI_methods = pd.read_csv(dir_path + '/Bayesian_Network/Methods_TRI.csv',
                                usecols = ['Code 2004 and prior', 'Objective'])
    df_TRI_methods.rename(columns = {'Code 2004 and prior': 'PCU'}, inplace = True)
    Prob_PCU_chemical_and_stream = pd.merge(df_TRI_methods, Prob_PCU_chemical_and_stream,
                                            on = ['PCU'], how = 'inner')
    Prob_PCU_chemical_and_stream =  Prob_PCU_chemical_and_stream\
                                    .groupby('Stream', as_index = False)\
                                    .apply(lambda x: pairwise_comparison(x, objective = 'seq'))
    Prob_PCU_chemical_and_stream.to_csv(dir_path + '/Fuzzy_Analytical_Hierarchy_Process/PCU_selection_and_position_under_FAHP.csv',
                                        sep = ',', index = False)
    # Chemical flow analysis
