# -*- coding: utf-8 -*-
# !/usr/bin/env python

# Importing libraries
from graphviz import Digraph
from chemical_flow_analysis.auxiliary import (estimating_val_with_log,
                                              emission_factor,
                                              method_of_moments,
                                              non_zero_output_streams,
                                              checking_outliers)
import pandas as pd
import numpy as np

def building_pau_black_box(chemical, Solubility, Tf, Tb, CoU, flow_inputs, waste_stream, PAU_code, efficiency, per_remaining, dir_path, Metal_indicator='NO'):
    Releases = pd.read_csv(dir_path + '/chemical_flow_analysis/tri_releases/TRI_releases.csv')
    Releases = Releases.loc[Releases['CAS NUMBER'] == chemical]
    for key, val in CoU.items():
        Releases = Releases.loc[Releases[key] == val]
    Releases.drop(columns=['CAS NUMBER', 'As a byproduct',
                           'As a manufactured impurity',
                           'As a process impurity'], inplace=True)
    Duplicated_columns = ['Reporting year', 'TRIFID',
                          'Maximum amount on-site',
                          'Total waste', 'Total release']
    Releases = Releases.loc[Releases.duplicated(subset=Duplicated_columns,
                            keep=False)]
    Releases.sort_values(by=Duplicated_columns, inplace=True)
    # Allocation
    Allocation = pd.read_csv(dir_path + '/Methods_TRI.csv')
    Allocation = Allocation.where(pd.notnull(Allocation), None)
    Carrier_for_PAU = Allocation.loc[Allocation['Code 2004 and prior'] == PAU_code, 'Output carrier'].iloc[0]
    Objective = Allocation.loc[Allocation['Code 2004 and prior'] == PAU_code, 'Objective'].iloc[0]
    Type_of_WM = Allocation.loc[Allocation['Code 2004 and prior'] == PAU_code, 'Type of waste management'].iloc[0]
    Type_of_treatment = Allocation.loc[Allocation['Code 2004 and prior'] == PAU_code, 'If it is treatment, what kind of?'].iloc[0]
    try:
        # Checking flows
        mu_flow_in, theta_2_flow_in = method_of_moments(flow_inputs)
        ## Emission factors
        Releases['Emission_factor'] = Releases.apply(lambda x: np.array([
                                                    emission_factor(
                                                    x['Flow to compartment'],
                                                    x['Maximum amount on-site'],
                                                    x['Total waste'],
                                                    x['Total release']) for i in range(100)]),
                                                    axis = 1)
        Fugitive_emissions_factors = Releases.loc[Releases['Compartment'] == 'Fugitive air emission', 'Emission_factor']
        n_rows = Fugitive_emissions_factors.shape[0]
        mu_Fugitive_emissions_factors, theta_2_Fugitive_emissions_factors = method_of_moments(pd.Series([v for val in Fugitive_emissions_factors for v in val]))
        Stack_emission_factors = Releases.loc[Releases['Compartment'] == 'Stack air emission', 'Emission_factor']
        mu_Stack_emission_factors, theta_2_Stack_emission_factors = method_of_moments(pd.Series([v for val in Stack_emission_factors for v in val]))
        Discharge_emission_factors = Releases.loc[Releases['Compartment'] == 'On-site surface water', 'Emission_factor']
        mu_Discharge_emission_factors, theta_2_Discharge_emission_factors = method_of_moments(pd.Series([v for val in Discharge_emission_factors for v in val]))
        Soil_emission_factors = Releases.loc[Releases['Compartment'] == 'On-site soil', 'Emission_factor']
        mu_Soil_emission_factors, theta_2_Soil_emission_factors = method_of_moments(pd.Series([v for val in Soil_emission_factors for v in val]))
        del Releases, Fugitive_emissions_factors, Stack_emission_factors, Discharge_emission_factors, Soil_emission_factors
        Flows_out = pd.DataFrame()
        for i in range(n_rows):
            flow_input = estimating_val_with_log(mu_flow_in, theta_2_flow_in)
            Fugitive_emission_factor = estimating_val_with_log(mu_Fugitive_emissions_factors, theta_2_Fugitive_emissions_factors)
            Stack_emission_factor = estimating_val_with_log(mu_Stack_emission_factors, theta_2_Stack_emission_factors)
            Discharge_emission_factor = estimating_val_with_log(mu_Discharge_emission_factors, theta_2_Discharge_emission_factors)
            Soil_emission_factor = estimating_val_with_log(mu_Soil_emission_factors, theta_2_Soil_emission_factors)
            ## Carrier
            if Objective == 'Removal':
                Carrier_flow = flow_input*efficiency/100
                Destroyed_converted_degraded_flow = 0.0
                # Carrier analysis
                if '/' in Carrier_for_PAU:
                    if Carrier_for_PAU == 'L/W':
                        if PAU_code in ['A03', 'A04', 'H103']:
                            if Solubility >= 1000: # Soluble in water if 1,000-10,000 and very solubla if >= 10,000mg/L
                                Carrier = 'W'
                            else:
                                Carrier = 'L'
                        elif PAU_code == 'P31':
                            Carrier = waste_stream
                        else:
                            if waste_stream == 'L':
                                Carrier = 'W'
                            else:
                                Carrier = 'L'
                    else:
                        if waste_stream in ['L', 'W']:
                            Carrier = 'S'
                        else:
                            Carrier = 'W'
                else:
                    Carrier = Carrier_for_PAU
            elif Objective == 'Reclamation':
                Carrier_flow = flow_input*efficiency/100
                Destroyed_converted_degraded_flow = 0.0
                if PAU_code in ['R11', 'R12', 'R13','R14', 'H20']: # Solvent recovery
                    Carrier = 'L'
                elif PAU_code in ['R21', 'R22', 'R23', 'R24', 'R26', 'R27', 'R28', 'R29', 'H10']: # Metal recovery
                    Carrier = 'S'
                else:
                    if Metal_indicator == 'YES':
                        Carrier = 'S'
                    else:
                        Carrier = 'W'
            elif (Objective == 'Safer disposal (RCRA Subtitle C)') | (Objective == 'Controlling flow rate and composition'):
                Carrier_flow = 0.0
                Carrier = None
                Destroyed_converted_degraded_flow = 0.0
            else:
                Carrier_flow = 0.0
                Carrier = None
                Destroyed_converted_degraded_flow = flow_input*efficiency/100
            if efficiency == 0.0:
                Carrier = None
            # Hidden stream
            if per_remaining:
                F_out_hidden = per_remaining/100*flow_input
            else:
                F_out_hidden = (1 - Fugitive_emission_factor)*(1 - efficiency/100)*flow_input
            ## By product and Remaining
            if (Type_of_WM == 'Energy recovery') | (Type_of_treatment == 'Incineration'):
                if PAU_code == 'A01':
                    if Tb <= 60: # deg C chemical is assessed as a gas or liquid (Knock-out Drum)
                        Remaining_flow = F_out_hidden
                        Remaining = 'A'
                        By_product_flow = 0.0
                        By_product = None
                    else:
                        Remaining_flow = Stack_emission_factor*F_out_hidden
                        Remaining = 'A'
                        By_product_flow = (1 - Stack_emission_factor)*F_out_hidden
                        if Tf <= 850: # Assesing normal operating temperature of flares
                            By_product = 'L'
                        else:
                            By_product = 'S'
                else:
                    # Assesing normal operating temperature of these devices
                    if (Tb <= 980) | (Tf <= 980): # deg C chemical is assessed as a gas or liquid
                        Remaining_flow = F_out_hidden
                        Remaining = 'A'
                        By_product_flow = 0.0
                        By_product = None
                    else: # deg C Chemical is assessed as a solid
                        Remaining_flow = Stack_emission_factor*F_out_hidden
                        Remaining = 'A'
                        By_product_flow = (1 - Stack_emission_factor)*F_out_hidden
                        By_product = 'S'
            else:
                if PAU_code in ['H20', 'R11', 'R12', 'R13', 'R14', 'R40', 'H39']:
                    if waste_stream == 'A':
                        Remaining_flow = Stack_emission_factor*F_out_hidden
                        Remaining = 'A'
                    else:
                        Remaining_flow = (1 - Discharge_emission_factor)*F_out_hidden
                        Remaining = waste_stream
                    By_product_flow = F_out_hidden - Remaining_flow
                    By_product = 'W'
                elif PAU_code in ['H10', 'R21', 'R22', 'R23', 'R24', 'R26']:
                    By_product_flow = Discharge_emission_factor*F_out_hidden
                    By_product = 'W'
                    Remaining_flow =  F_out_hidden - By_product_flow
                    Remaining = waste_stream
                elif PAU_code in ['R27', 'R28', 'R29']:
                    By_product_flow = Soil_emission_factor*F_out_hidden
                    By_product = 'S'
                    Remaining_flow =  F_out_hidden - By_product_flow
                    Remaining = waste_stream
                else:
                    Remaining_flow =  F_out_hidden
                    Remaining = waste_stream
                    By_product_flow = 0.0
                    By_product = None
            if per_remaining:
                Fugitive_flow = flow_input - Destroyed_converted_degraded_flow - Carrier_flow - Remaining_flow - By_product_flow
            else:
                Fugitive_flow = Fugitive_emission_factor*(1 - efficiency/100)*flow_input
            Flows_out_aux = df_result_aux = pd.DataFrame({'Destroyed/converted/degraded': [Destroyed_converted_degraded_flow],
                                                          'Carrier': [[Carrier_flow, Carrier]],
                                                          'Remaining': [[Remaining_flow, Remaining]],
                                                          'By-product': [[By_product_flow, By_product]],
                                                          'Fugitive emission': [Fugitive_flow]})
            Flows_out = pd.concat([Flows_out, Flows_out_aux], axis = 0,
                            sort = False,
                            ignore_index = True)
        Cols = ['PAU', 'Chemical', 'Type of stream', 'Mean quantity [kg/yr]', 'Std [kg/yr]', 'Number of data', 'Phase']
        df_result = pd.DataFrame(columns = Cols)
        for col in Flows_out.columns:
            Result = {'Type of stream': [col]}
            if col in ['Destroyed/converted/degraded', 'Fugitive emission']:
                Val =  Flows_out[col]
                if col == 'Destroyed/converted/degraded':
                    Phase = 'NA'
                else:
                    Phase = 'A'
            else:
                Val = pd.Series([val[0] for val in Flows_out[col]])
                Phase = Flows_out[col][0][1]
            Val = checking_outliers(Val)
            Mean = Val.mean()
            Number =  Val.shape[0]
            Std = Val.std()
            Result.update({'PAU': [PAU_code],
                          'Chemical': [chemical],
                          'Mean quantity [kg/yr]': [Mean],
                          'Std [kg/yr]': [Std],
                          'Number of data': [Number],
                          'Phase': [Phase]})
            df_result_aux = pd.DataFrame(Result)
            df_result = pd.concat([df_result, df_result_aux], axis = 0,
                            sort = False,
                            ignore_index = True)
        # Returning values
        return df_result, checking_outliers(pd.Series([val[0] for val in Flows_out['Remaining']]))
    except ZeroDivisionError:
        df_result = pd.DataFrame({'PAU': [PAU_code],
                                  'Chemical': [chemical],
                                  'Mean quantity [kg/yr]': ['NA'],
                                  'Std [kg/yr]': ['NA'],
                                  'Number of data': ['NA'],
                                  'Phase': ['NA']})
        return df_result,  pd.Series([0.0])



def picture(df_for_stream, stream, dir_path, drawing, paus_stream):
    # Organizing information
    n_paus = df_for_stream['PAU'].nunique()
    PAU_position = {row['Position']: row['PAU'] for idx, row in df_for_stream[['PAU', 'Position']].drop_duplicates(keep = 'first').iterrows()}
    Chemical_PAU = {value: df_for_stream.loc[df_for_stream['PAU'] == value, 'Chemical'].tolist() for value in PAU_position.values()}
    CoU_columns = ['As a manufactured impurity',
                   'As a process impurity',
                   'As a byproduct']
    Chemicals_in_stream = df_for_stream['Chemical'].unique().tolist()
    CoU_chemical = {chem['Chemical']: {col: chem[col] for col in CoU_columns} for idx, chem in df_for_stream[['Chemical'] + CoU_columns].drop_duplicates(keep = 'first').iterrows()}
    Solubility_chemical = {chem['Chemical']: chem['Solubility'] for idx, chem in df_for_stream[['Chemical', 'Solubility']].drop_duplicates(keep = 'first').iterrows()}
    Tf_chemical = {chem['Chemical']: chem['Tf'] for idx, chem in df_for_stream[['Chemical', 'Tf']].drop_duplicates(keep = 'first').iterrows()}
    Tb_chemical = {chem['Chemical']: chem['Tb'] for idx, chem in df_for_stream[['Chemical', 'Tb']].drop_duplicates(keep = 'first').iterrows()}
    Eff_chemical = {chem['Chemical']: chem['Efficiency [%]'] for idx, chem in df_for_stream[['Chemical', 'Efficiency [%]']].drop_duplicates(keep = 'first').iterrows()}
    Flow_in = {chem['Chemical']: chem['Chemical flow'] for idx, chem in df_for_stream[['Chemical', 'Chemical flow']].drop_duplicates(keep = 'first').iterrows()}
    Flow_in_aux = Flow_in.copy()
    waste_stream = df_for_stream['Type of waste'].iloc[0]
    # Creating graph
    PAU = Digraph()
    PAU.graph_attr['rankdir'] = 'LR'
    PAU.graph_attr['ranksep'] = '1'
    # Creating nodes
    ## Input
    PAU.node(str(0), 'Input stream', fontname = 'Times New Roman Bold', style = 'filled', color = 'lightsalmon', shape = 'ellipse')
    ## PAUs
    Flows = pd.DataFrame()
    for position in PAU_position.keys():
        PAU_code = PAU_position[position]
        PAU.node(str(position), PAU_code, fontname = 'Times New Roman Bold', shape = 'square')
        for chemical in Chemicals_in_stream:
            if chemical in Chemical_PAU[PAU_code]:
                efficiency = Eff_chemical[chemical]
                per_remaining = None
            else:
                pau_for_chem = [pos for paus_1, chems in Chemical_PAU.items() for pos, paus_2 in PAU_position.items() if (chemical in chems) and (paus_1 == paus_2)][0]
                efficiency = 0.0
                if pau_for_chem > position:
                    per_remaining = Eff_chemical[chemical]
                else:
                    per_remaining = None
            CoU = CoU_chemical[chemical]
            Solubility = Solubility_chemical[chemical]
            Tf = Tf_chemical[chemical]
            Tb = Tb_chemical[chemical]
            flow_input = Flow_in[chemical]
            Result, flow_effluent = building_pau_black_box(chemical, Solubility, Tf, Tb, CoU, flow_input, waste_stream, PAU_code, efficiency, per_remaining, dir_path)
            Result['Position'] = position
            Flows = pd.concat([Flows, Result], axis = 0,
                        sort = False,
                        ignore_index = True)
            Flow_in[chemical] = flow_effluent
    Flows['Stream'] = stream
    Flows.sort_values(by = ['Position', 'Chemical'],
                      ascending = [True, True],
                      inplace = True)
    Flows['Mean quantity [kg/yr]'] = Flows['Mean quantity [kg/yr]'].round(4)
    Flows['Std [kg/yr]'] = Flows['Std [kg/yr]'].round(4)
    ## Output
    PAU.node(str(n_paus + 1), 'Remaining stream', fontname = 'Times New Roman Bold', style = 'filled', color = 'lightblue', shape = 'ellipse')
    # Searching non zero streams for PAU
    Output_streams = dict()
    for Pos, PAU_code in PAU_position.items():
        PAU_flows = Flows.loc[Flows['PAU'] == PAU_code]
        Output_streams.update({Pos: non_zero_output_streams(PAU_flows)})
    # Creating edges
    n_streams = 0
    n_node = n_paus + 1
    for i in range(n_paus + 1):
        if i in Output_streams.keys():
            Flow_graph = Digraph()
            Flow_graph.graph_attr['rank'] = 'same'
            Flow_graph.graph_attr['ranksep'] = '2'
            for strm in Output_streams[i]:
                n_node = n_node + 1
                n_streams = n_streams + 1
                if strm in ['Fugitive emission', 'By-product']:
                    Flow_graph.node(str(n_node), strm, fontname = 'Times New Roman Bold',
                                    style = 'filled', color = 'lightgray', shape = 'ellipse')
                    Flow_graph.edge(str(n_node), str(i), label = str(n_streams),
                                    len = '1', arrowhead = 'normal', arrowsize = '1',
                                    dir = 'back', fontname = 'Comic Sans MS Bold')
                elif strm in ['Destroyed/converted/degraded', 'Carrier']:
                    Flow_graph.node(str(n_node), strm,
                                    fontname='Times New Roman Bold',
                                    style='filled', color='#b2df8a',
                                    shape='ellipse')
                    Flow_graph.edge(str(i), str(n_node), label=str(n_streams),
                                    len='1', arrowhead='normal', arrowsize='1',
                                    fontname='Comic Sans MS Bold')
            PAU.subgraph(Flow_graph)
        n_streams = n_streams + 1
        PAU.edge(str(i), str(i + 1), label=str(n_streams),
                 len='0.5', arrowhead='normal', arrowsize='1',
                 fontname='Comic Sans MS Bold')
    if drawing:
        PAU.view(filename=f'{dir_path}/chemical_flow_analysis/pau_draws/Pollution_abatement_for_{paus_stream}')
    return Flows
