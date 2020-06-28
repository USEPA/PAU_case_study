# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Importing libraries
import argparse
import os
import xlrd, xlwt
from xlutils.copy import copy
from Bayesian_Network.Bayesian_Network import *

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
    Columns = set(df_inputs.columns.tolist()) - set(['CAS NUMBER', 'Stream #', 'Chemical flow [kg/yr]'])
    for chem_1, chem_2 in CAS_for_search.items():
        Input_dictionary_joint = dict()
        Input_dictionary_marginal = dict()
        for col in Columns:
            df_input_chem = df_inputs[df_inputs['CAS NUMBER'] == chem_1]
            val = df_input_chem[col].values[0]
            Input_dictionary_joint.update({col: val})
            if col != 'Type of waste management':
                Input_dictionary_marginal.update({col: val})
        PCU_model = Building_bayesian_network_model(dir_path, df_PCU, chem_1, chem_2)
        Calculating_joint_probabilities(Input_dictionary_joint, chem_1, chem_2, PCU_model, dir_path, df_PCU)
        Calculating_marginal_probabilities(Input_dictionary_marginal, PCU_model, dir_path, chem_1, chem_2, df_PCU)
    # Fuzzy analysis
