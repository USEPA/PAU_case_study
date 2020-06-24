# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Importing libraries
import pandas as pd
pd.set_option('mode.chained_assignment', None)
import math
from pomegranate import *
import argparse
import os
import itertools
import pygraphviz as pgv

from interval import interval, inf, imath


def building_network(Probabilities, edges):
    PCU_model = BayesianNetwork('PCU Selection')
    # Adding states
    nodes =  dict()
    for Name, Prob in Probabilities.items():
        s = Node(Prob, name = Name)
        nodes.update({Name: s})
        PCU_model.add_node(s)
    # Adding edges
    for idx, row in edges.iterrows():
        PCU_model.add_edge(nodes[row['Source']], nodes[row['Target']])
    PCU_model.bake()
    return PCU_model


def drawing_network(df, dir_path):
    df_names = pd.read_csv(dir_path + '/Node_names.csv')
    names = {row['Name']:row['Node'] for idx, row in df_names.iterrows()}
    df = df.applymap(lambda x: names[x])
    Graph_bayesian = pgv.AGraph(strict = True, directed = True, ranksep = '3')
    values = [val for val in set([col for row in  df.values for col in row])]
    Graph_bayesian.add_nodes_from(values, color = '#99b898', shape = 'circle',
                                fontsize = '40', fontweight = 'bold',
                                fontfamily = 'arial', width = '3',
                                style = 'bold', penwidth = '5')
    for idx, row in df.iterrows():
        if row['Target'] == 'Node-17':
            color = '#1b6ca8'
        elif row['Target'] == 'Node-18':
            color = '#ff847c'
        elif row['Source'] == 'Node-17':
            color = '#e84a5f'
        else:
            color = 'black'
        Graph_bayesian.add_edge(*tuple(row.values) , color = color, alpha = '0.2',
                            arrowsize = '2', arrowType = 'open', style = 'bold')
    Graph_bayesian.draw(dir_path + '/Bayesian_Network_PCU.png' ,prog = 'dot')


def building_dataframe(dir_path, Years, values):
    cols_for_using = ['TRIFID', 'CAS NUMBER', 'PRODUCE THE CHEMICAL', 'IMPORT THE CHEMICAL',
    	               'SALE OR DISTRIBUTION OF THE CHEMICAL',
                       'AS A BYPRODUCT', 'AS A MANUFACTURED IMPURITY', 'USED AS A REACTANT',
                       'ADDED AS A FORMULATION COMPONENT', 'USED AS AN ARTICLE COMPONENT',
                       'REPACKAGING', 'AS A PROCESS IMPURITY', 'USED AS A CHEMICAL PROCESSING AID',
                       'USED AS A MANUFACTURING AID', 'ANCILLARY OR OTHER USE',
                       'WASTE STREAM CODE', 'RANGE INFLUENT CONCENTRATION',
                       'METHOD CODE - 2004 AND PRIOR', 'EFFICIENCY RANGE CODE',
                       'TYPE OF MANAGEMENT', 'PRIMARY NAICS CODE']
    df = pd.DataFrame()
    for year in Years:
        df_year = pd.read_csv(dir_path + '/Final_PCU_datasets/PCUs_DB_filled_{}.csv'.format(year),
                            usecols = cols_for_using, low_memory = False,
                            dtype = {'PRIMARY NAICS CODE': 'object'})
        df_year['REPORTING YEAR'] = str(year)
        df = pd.concat([df, df_year], axis = 0,
                        sort = False,
                        ignore_index = True)
        del df_year
    df.drop_duplicates(keep = 'first', inplace = True)
    df = df[~df['METHOD CODE - 2004 AND PRIOR'].str.contains('\+')]
    cols_for_change = {'WASTE STREAM CODE': 'Type of waste',
                    'RANGE INFLUENT CONCENTRATION': 'Concentration',
                    'METHOD CODE - 2004 AND PRIOR': 'PCU',
                    'EFFICIENCY RANGE CODE': 'Efficiency',
                    'TYPE OF MANAGEMENT': 'Type of waste management'}
    cols_for_change.update({col: col.capitalize() for col in set(cols_for_using) - set(cols_for_change.keys()) if not col in \
                            ['CAS NUMBER', 'PRIMARY NAICS CODE', 'TRIFID', 'REPORTING YEAR']})
    df.rename(columns = cols_for_change, inplace = True)
    for key, val in values.items():
        df[key] = 'No'
        df.loc[(df[val] == 'Yes').any(axis = 1), key] = 'Yes'
    df['Concentration'] = df['Concentration'].apply(lambda x: str(int(abs(x))))
    df = df.loc[df['Type of waste'] != 'X']
    df = df.applymap(lambda x: x.strip())
    return df


def buidling_probabilities(df, edges, Inputs_list):
    children = edges['Target'].unique().tolist()
    parents = {child: edges.loc[edges['Target'] == child, 'Source'].tolist() for child in children}
    without_parents = edges.loc[~edges['Source'].isin(children), 'Source'].unique().tolist()
    Probabilities = {root: DiscreteDistribution({key: value for key, value in df[root].value_counts(normalize = True).items()}) for root in without_parents}
    Conditional_probabilities = dict()
    for key, value in parents.items():
        contingency_table = pd.crosstab(df[key], [df[col] for col in value], normalize = 'columns').T
        if len(value) != 1:
            level_values = tuple(Inputs_list[val] for val in contingency_table.index.names)
            combinations = [list(com) for com in itertools.product(*level_values)]
        else:
            level_values = list(Inputs_list[contingency_table.index.name])
            combinations = [[com] for com in level_values]
        column_names = Inputs_list[key]
        list_of_lists = list()
        for combination in combinations:
            for column_name in column_names:
                list_aux = combination.copy()
                try:
                    if len(value) == 1:
                        frac = contingency_table.loc[combination, column_name].values[0]
                    else:
                        frac = frac = contingency_table.xs(tuple(combination), level = tuple(value))[column_name].values[0]
                except (KeyError, IndexError):
                    frac = 0
                list_aux  = list_aux + [column_name, frac]
                list_of_lists.append(list_aux)
        Probabilities.update({key: ConditionalProbabilityTable(list_of_lists, [Probabilities[val] for val in value])})
    return Probabilities


def Calculating_probabilities(Input_dictionary, chem, Model, dir_path):
    Input_dictionary.update({'Chemical': chem})
    PCUs = set([vals[2] for state in Model.states for vals in state.distribution.parameters[0] if state.name == 'PCU'])
    value_event = [Input_dictionary[state.name] for state in Model.states if state.name != 'PCU']
    values_to_assess = [value_event + [PCU] for PCU in PCUs]
    Probabilities = Model.probability(values_to_assess)
    df_values = {key: [value]*len(PCUs) for key, value in Input_dictionary.items()}
    df_values.update({'PCU': list(PCUs), 'Probability': Probabilities})
    df = pd.DataFrame(df_values)
    df.to_csv(dir_path + '/Probabilities_based_on_BN_for_{}.csv'.format(chem), sep = ',', index = False)


def Building_flows_dataset(dir_path, Years, nbins, df_PCU, CAS_for_search):
    df =  pd.DataFrame()
    for Year in Years:
        df_year = pd.read_csv(dir_path + '/Waste_flow/Waste_flow_to_PCUs_{}_10.csv'.format(Year),
                     usecols = ['METHOD CODE', 'MIDDLE WASTE FLOW', 'TRIFID',
                                'WASTE STREAM CODE', 'PRIMARY NAICS CODE',
                                'CAS NUMBER'],
                     low_memory = False, dtype = {'PRIMARY NAICS CODE': 'object'})
        df_year['REPORTING YEAR'] = str(Year)
        df_year['MIDDLE WASTE FLOW'] = df_year.groupby(['TRIFID', 'METHOD CODE',\
                                                        'WASTE STREAM CODE',
                                                        'PRIMARY NAICS CODE'],\
                                                        as_index = False)['MIDDLE WASTE FLOW']\
                                             .transform('mean')
        df = pd.concat([df, df_year], axis = 0,
                        sort = False,
                        ignore_index = True)
        del df_year
    df.drop_duplicates(keep = 'first', inplace = True)
    df.rename(columns = {'METHOD CODE': 'PCU',
                         'WASTE STREAM CODE': 'Type of waste'},
              inplace = True)
    df = pd.merge(df_PCU, df, on = ['TRIFID', 'CAS NUMBER',
                                         'PRIMARY NAICS CODE',
                                         'PCU', 'Type of waste',
                                         'REPORTING YEAR'],
                         how = 'left')
    df['MIDDLE WASTE FLOW'] = df['MIDDLE WASTE FLOW'].fillna(df.groupby(['PCU', 'PRIMARY NAICS CODE'])\
                                                            ['MIDDLE WASTE FLOW'].transform('mean'))
    df['MIDDLE WASTE FLOW'] = df['MIDDLE WASTE FLOW'].fillna(df.groupby(['PCU'])\
                                                            ['MIDDLE WASTE FLOW'].transform('mean'))
    df = df.loc[pd.notnull(df).all(axis = 1)]
    Max_value = df['MIDDLE WASTE FLOW'].max()
    Min_value = df['MIDDLE WASTE FLOW'].min()
    Order_max = int(math.log10(Max_value)) - 1
    Order_min = math.ceil(math.log10(Min_value))
    Delta = (Order_max - Order_min)/(nbins - 2)
    Bin_values = [Min_value - 10**(math.log10(Min_value) - 1)]
    Bin_values = Bin_values + [10**(Order_min + Delta*n) for n in range(nbins - 1)]
    Bin_values = Bin_values + [Max_value]
    Bin_values.sort()
    Bin_labels = [str(val) for val in range(1, len(Bin_values))]
    df['MIDDLE WASTE FLOW INTERVAL'] = pd.cut(df['MIDDLE WASTE FLOW'],
                                              bins = Bin_values)
    df['MIDDLE WASTE FLOW INTERVAL CODE'] = pd.cut(df['MIDDLE WASTE FLOW'],
                                                   bins = Bin_values,
                                                   labels = Bin_labels,
                                                   precision = 0)
    df_intervals = df[['MIDDLE WASTE FLOW INTERVAL', 'MIDDLE WASTE FLOW INTERVAL CODE']]
    df_intervals.drop_duplicates(keep = 'first', inplace = True)
    df_intervals['UNIT'] = 'kg/yr'
    df_intervals.sort_values(by = ['MIDDLE WASTE FLOW INTERVAL CODE'],
                            ascending =  True, inplace = True)
    df_intervals.to_csv(dir_path + '/Relationship_flow_interval_and_codes.csv',
                        sep = ',', index = False)
    del df_intervals
    df = df[df['CAS NUMBER'].isin(CAS_for_search)]
    df.drop(columns = ['MIDDLE WASTE FLOW INTERVAL', 'MIDDLE WASTE FLOW',
                       'REPORTING YEAR'],
            inplace = True)
    df.rename(columns = {'MIDDLE WASTE FLOW INTERVAL CODE': 'Waste flow'},
              inplace = True)
    return df


def Building_price_dataset(dir_path, Years, nbins, df_PCU):
    df =  pd.DataFrame()
    for Year in Years:
        df_year = pd.read_csv(dir_path + '/Chemical_price/Chemical_price_vs_PCU_{}.csv'.format(Year),
                     low_memory = False,
                     usecols = ['METHOD CODE - 2004 AND PRIOR', 'TRIFID',
                                'UNIT PRICE (USD/g)'])
        df_year = df_year.groupby(['METHOD CODE - 2004 AND PRIOR', 'TRIFID'],
                        as_index = False).mean()
        df = pd.concat([df, df_year], axis = 0,
                        sort = False,
                        ignore_index = True)
        del df_year
    df.drop_duplicates(keep = 'first', inplace = True)
    df.rename(columns = {'METHOD CODE - 2004 AND PRIOR': 'PCU'},
              inplace = True)
    df = pd.merge(df_PCU, df, how = 'left', on = ['PCU', 'TRIFID'])
    df['UNIT PRICE (USD/g)'] = df['UNIT PRICE (USD/g)'].fillna(df.groupby(['PCU'])\
                                                            ['UNIT PRICE (USD/g)'].transform('mean'))
    df = df.loc[pd.notnull(df).all(axis = 1)]
    Max_value = df['UNIT PRICE (USD/g)'].max()
    Min_value = df['UNIT PRICE (USD/g)'].min()
    Order_max = int(math.log10(Max_value)) - 1
    Order_min = math.ceil(math.log10(Min_value))
    Delta = (Order_max - Order_min)/(nbins - 2)
    Bin_values = [Min_value - 10**(math.log10(Min_value) - 1)]
    Bin_values = Bin_values + [10**(Order_min + Delta*n) for n in range(nbins - 1)]
    Bin_values = Bin_values + [Max_value]
    Bin_values.sort()
    Bin_labels = [str(val) for val in range(1, len(Bin_values))]
    df['PRICE INTERVAL'] = pd.cut(df['UNIT PRICE (USD/g)'],
                                  bins = Bin_values)
    df['PRICE INTERVAL CODE'] = pd.cut(df['UNIT PRICE (USD/g)'],
                                                   bins = Bin_values,
                                                   labels = Bin_labels,
                                                   precision = 0)
    df_intervals = df[['PRICE INTERVAL', 'PRICE INTERVAL CODE']]
    df_intervals.drop_duplicates(keep = 'first', inplace = True)
    df_intervals['UNIT'] = 'USD/g'
    df_intervals.sort_values(by = ['PRICE INTERVAL CODE'],
                            ascending =  True, inplace = True)
    df_intervals.to_csv(dir_path + '/Relationship_chemical_prices_and_codes.csv',
                        sep = ',', index = False)
    del df_intervals
    df.drop(columns = ['PRICE INTERVAL', 'UNIT PRICE (USD/g)', 'TRIFID'], inplace = True)
    df.rename(columns = {'PRICE INTERVAL CODE': 'Chemical price',
                        'PRIMARY NAICS CODE': 'NAICS code'},
              inplace = True)
    return df


def Building_PAOC_and_PACE_dataset(dir_path, nbins, ccccccc):
    df_PAOC = pd.read_csv(dir_path + '/PCU_expenditure_and_cost/PAOC.csv',
                          usecols = ['Activity', 'Media', 'Mean PAOC', 'NAICS code'],
                          dtype = {'NAICS code': 'object'})
    df_PAOC.rename(columns = {'Activity': 'Type of waste management',
                              'Media': 'Type of waste'},
                    inplace = True)
    df_PAOC = pd.merge(df_PAOC, df_PCU, how = 'right', on = ['Type of waste management',
                                                   'Type of waste',
                                                    'NAICS code'])
    df_PAOC['Mean PAOC'] = df_PAOC['Mean PAOC'].fillna(df_PAOC\
                                                    .groupby(['Type of waste', 'Type of waste management'])\
                                                    ['Mean PAOC']\
                                                    .transform('mean'))
    df_PAOC['Mean PAOC'] = df_PAOC['Mean PAOC'].fillna(df_PAOC\
                                                    .groupby(['Type of waste management'])\
                                                    ['Mean PAOC']\
                                                    .transform('mean'))
    df_PAOC = df_PAOC.loc[pd.notnull(df_PAOC).all(axis = 1)]
    Max_value = df_PAOC['Mean PAOC'].max()
    Min_value = df_PAOC['Mean PAOC'].min()
    Order_max = int(math.log10(Max_value)) - 1
    Order_min = math.ceil(math.log10(Min_value))
    Delta = (Order_max - Order_min)/(nbins - 2)
    Bin_values = [Min_value - 10**(math.log10(Min_value) - 1)]
    Bin_values = Bin_values + [10**(Order_min + Delta*n) for n in range(nbins - 1)]
    Bin_values = Bin_values + [Max_value]
    Bin_values.sort()
    Bin_labels = [str(val) for val in range(1, len(Bin_values))]
    df_PAOC['MEAN PAOC INTERVAL'] = pd.cut(df_PAOC['Mean PAOC'],
                                  bins = Bin_values)
    df_PAOC['MEAN PAOC INTERVAL CODE'] = pd.cut(df_PAOC['Mean PAOC'],
                                                   bins = Bin_values,
                                                   labels = Bin_labels,
                                                   precision = 0)
    df_intervals = df_PAOC[['MEAN PAOC INTERVAL', 'MEAN PAOC INTERVAL CODE']]
    df_intervals.drop_duplicates(keep = 'first', inplace = True)
    df_intervals['UNIT'] = 'USD/kg'
    df_intervals.sort_values(by = ['MEAN PAOC INTERVAL CODE'],
                            ascending =  True, inplace = True)
    df_intervals.to_csv(dir_path + '/Relationship_PAOC_and_codes.csv',
                        sep = ',', index = False)
    df_PAOC.rename(columns = {'MEAN PAOC INTERVAL CODE': 'PAOC'},
                    inplace = True)
    df_PAOC.drop(columns = ['Mean PAOC', 'MEAN PAOC INTERVAL'], inplace = True)
    df_PACE = pd.read_csv(dir_path + '/PCU_expenditure_and_cost/PACE.csv',
                          usecols = ['Activity', 'Media', 'Mean PACE', 'NAICS code'],
                          dtype = {'NAICS code': 'object'})
    df_PACE.rename(columns = {'Activity': 'Type of waste management',
                              'Media': 'Type of waste'},
                    inplace = True)
    df_PACE = pd.merge(df_PAOC, df_PACE, how = 'left', on = ['Type of waste management',
                                                   'Type of waste',
                                                    'NAICS code'])
    df_PACE['Mean PACE'] = df_PACE['Mean PACE'].fillna(df_PACE\
                                                    .groupby(['Type of waste', 'Type of waste management'])\
                                                    ['Mean PACE']\
                                                    .transform('mean'))
    df_PACE['Mean PACE'] = df_PACE['Mean PACE'].fillna(df_PACE\
                                                    .groupby(['Type of waste management'])\
                                                    ['Mean PACE']\
                                                    .transform('mean'))
    df_PACE = df_PACE.loc[pd.notnull(df_PACE).all(axis = 1)]
    Max_value = df_PACE['Mean PACE'].max()
    Min_value = df_PACE['Mean PACE'].min()
    Order_max = int(math.log10(Max_value)) - 1
    Order_min = math.ceil(math.log10(Min_value))
    Delta = (Order_max - Order_min)/(nbins - 2)
    Bin_values = [Min_value - 10**(math.log10(Min_value) - 1)]
    Bin_values = Bin_values + [10**(Order_min + Delta*n) for n in range(nbins - 1)]
    Bin_values = Bin_values + [Max_value]
    Bin_values.sort()
    Bin_labels = [str(val) for val in range(1, len(Bin_values))]
    df_PACE['MEAN PACE INTERVAL'] = pd.cut(df_PACE['Mean PACE'],
                                  bins = Bin_values)
    df_PACE['MEAN PACE INTERVAL CODE'] = pd.cut(df_PACE['Mean PACE'],
                                                   bins = Bin_values,
                                                   labels = Bin_labels,
                                                   precision = 0)
    df_intervals = df_PACE[['MEAN PACE INTERVAL', 'MEAN PACE INTERVAL CODE']]
    df_intervals.drop_duplicates(keep = 'first', inplace = True)
    df_intervals['UNIT'] = 'USD/kg'
    df_intervals.sort_values(by = ['MEAN PACE INTERVAL CODE'],
                            ascending =  True, inplace = True)
    df_intervals.to_csv(dir_path + '/Relationship_PACE_and_codes.csv',
                        sep = ',', index = False)
    df_PACE.rename(columns = {'MEAN PACE INTERVAL CODE': 'PACE'},
                    inplace = True)
    df_PACE.drop(columns = ['Mean PACE', 'MEAN PACE INTERVAL', 'NAICS code'], inplace = True)
    return df_PACE


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

    parser.add_argument('-N_Bins',
                        help = 'Number of bins to split the middle waste flow/prices/PAOC/PACE',
                        type = int,
                        default =  10,
                        required = False)

    args = parser.parse_args()

    CAS = args.CAS
    Years = args.Year
    N_Bins = args.N_Bins
    dir_path = os.path.dirname(os.path.realpath(__file__)) # Current directory

    df_chemicals = pd.read_csv(dir_path + '/Chemicals/Chemicals.csv',
                                usecols = ['CAS NUMBER'])
    CAS_for_search = list()
    for chem in CAS:
        Value = (df_chemicals['CAS NUMBER'] == chem).any()
        if Value:
            print('Chemical with CAS Number {} is in the TRI Program'.format(chem))
            CAS_for_search.append(chem)
        else:
            print('Chemical with CAS Number {} is not in the TRI Program'.format(chem))
    if CAS_for_search:
        Options = [['Yes', 'No'], ['1', '2', '3', '4','5'],
                ['W', 'S', 'L', 'A'],
                ['E1', 'E2','E3', 'E4', 'E5', 'E6']]
        Inputs_list = {'Manufacturing stage': Options[0], 'Produce the chemical': Options[0],
                    'Import the chemical' : Options[0], 'Sale or distribution of the chemical': Options[0],
                    'As a byproduct' : Options[0], 'As a manufactured impurity' : Options[0],
                    'Processing stage': Options[0], 'Used as an article component': Options[0],
                    'Added as a formulation component': Options[0],  'Used as a reactant': Options[0],
                    'As a process impurity': Options[0], 'Repackaging': Options[0],
                    'Use stage': Options[0], 'Used as a chemical processing aid': Options[0],
                    'Used as a manufacturing aid': Options[0], 'Ancillary or other use': Options[0],
                    'Concentration': Options[1], 'Type of waste': Options[2], 'Efficiency': Options[3],
                    'Type of waste management': ['Recycling', 'Energy recovery', 'Treatment']}
        edges = pd.read_csv(dir_path + '/Graph.csv')
        drawing_network(edges, dir_path)
        values = [val for val in set([col for row in  edges.values for col in row]) if 'stage' in val]
        values = {val: edges.loc[edges['Source'] == val, 'Target'].tolist()  for val in values}
        try:
            # PCU
            df_PCU = building_dataframe(dir_path, Years, values)
            # Waste flows
            df_PCU = Building_flows_dataset(dir_path, Years, N_Bins, df_PCU, CAS_for_search)
            # Chemical prices
            df_PCU = Building_price_dataset(dir_path, Years, N_Bins, df_PCU)
            # Pollution abatement operating cost (PAOC)
            # Pollution abatement capital expenditure(PACE)
            df_PCU = Building_PAOC_and_PACE_dataset(dir_path, N_Bins, df_PCU)
            df_PCU.to_csv(dir_path + '/DB_Bayesian_Network.csv',
                          sep = ',', index = False)
        except FileNotFoundError:
            df_PCU = pd.read_csv(dir_path + '/DB_Bayesian_Network.csv', low_memory = False)
        for chem in CAS_for_search:
            df_PCU_chem = df_PCU[df_PCU['CAS NUMBER'] == chem]
            df_PCU_chem.drop(columns = ['CAS NUMBER'], inplace = True)
            Levels_for_PCU = df_PCU_chem['PCU'].unique().tolist()
            df_1 = pd.read_csv(dir_path + '/Relationship_chemical_prices_and_codes.csv')
            Levels_for_chemical_prices = df_1['PRICE INTERVAL CODE'].tolist()
            Option_prices = [str(row['PRICE INTERVAL CODE']) + ': ' + str(row['PRICE INTERVAL']) for idx, row in df_1.iterrows()]
            del df_1
            df_2 = pd.read_csv(dir_path + '/Relationship_flow_interval_and_codes.csv')
            Levels_for_flow = df_2['MIDDLE WASTE FLOW INTERVAL CODE'].tolist()
            Option_flow = [str(row['MIDDLE WASTE FLOW INTERVAL CODE']) + ': ' + str(row['MIDDLE WASTE FLOW INTERVAL']) for idx, row in df_2.iterrows()]
            del df_2
            df_3 = pd.read_csv(dir_path + '/Relationship_PACE_and_codes.csv')
            Levels_for_PACE = df_3['MEAN PACE INTERVAL CODE'].tolist()
            Option_PACE = [str(row['MEAN PACE INTERVAL CODE']) + ': ' + str(row['MEAN PACE INTERVAL']) for idx, row in df_3.iterrows()]
            del df_3
            df_4 = pd.read_csv(dir_path + '/Relationship_PAOC_and_codes.csv')
            Levels_for_PAOC = df_4['MEAN PAOC INTERVAL CODE'].tolist()
            Option_PAOC = [str(row['MEAN PAOC INTERVAL CODE']) + ': ' + str(row['MEAN PAOC INTERVAL']) for idx, row in df_4.iterrows()]
            del df_4
            Inputs_list.update({'PCU': Levels_for_PCU, 'PAOC': Levels_for_PAOC,
                                'PACE': Levels_for_PACE, 'Waste flow': Levels_for_flow,
                                'Chemical price': Levels_for_chemical_prices})
            Probabilities = buidling_probabilities(df_PCU_chem, edges, Inputs_list)
            PCU_model = building_network(Probabilities, edges)
            Options = [['Yes', 'No'],
                      ['1: concentraton > 1%', '2: 0.01% < concentration <= 1%',
                      '3: 0.0001% < concentration <= 0.01%', '4: 0.0000001% < concentration <= 0.0001%',
                      '5: concentration <= 0.0000001%'],
                      ['W: wastewater', 'S: solid waste', 'L: liquid waste', 'A: gaseous waste'],
                      ['E1: efficiency > 99.9999%', 'E2: 99.99% < efficiency <= 99.9999%',
                      'E3: 99% < efficiency <= 99.99%', 'E4: 95% < efficiency <= 99%',
                      'E5: 50% < efficiency <= 95%', 'E6: 0% < efficiency <= 50%']]
            Inputs_list = {'Manufacturing stage': Options[0], 'Produce the chemical': Options[0],
                         'Import the chemical' : Options[0], 'Sale or distribution of the chemical': Options[0],
                         'As a byproduct' : Options[0], 'As a manufactured impurity' : Options[0],
                         'Processing stage': Options[0], 'Used as an article component': Options[0],
                         'Added as a formulation component': Options[0],  'Used as a reactant': Options[0],
                         'As a process impurity': Options[0], 'Repackaging': Options[0],
                         'Use stage': Options[0], 'Used as a chemical processing aid': Options[0],
                         'Used as a manufacturing aid': Options[0], 'Ancillary or other use': Options[0],
                         'Concentration': Options[1], 'Type of waste': Options[2], 'Efficiency': Options[3],
                         'Chemical price [USD/g]': Option_prices, 'Waste flow [kg/yr]': Option_flow,
                         'Pollution abatement capital expenditure [USD/kg]': Levels_for_PACE,
                         'Pollution Abatement operating cost [USD/kg]': Option_PAOC}
            # for chem in CAS_for_search:
            #     print('-'*120)
            #     print('\nFor the chemical with CAS Number {}, select an option for the following: '.format(chem))
            #     print()
            #     Inputs_dictionary = {Input: input('\n{}, options:\n\n{}\n\nInput: '.format(Input, '\n'.join(' - {}'.format(opt) for opt in option))) for Input, option in  Inputs_list.items()}
            #     print(Inputs_dictionary)
            #     Calculating_probabilities(Inputs_dictionary, chem, PCU_model, dir_path)
    else:
        print('It is not possible to build a Bayesian Network for the chemicals')
