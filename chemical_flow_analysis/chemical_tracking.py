# -*- coding: utf-8 -*-
# !/usr/bin/env python

# Importing libraries
from graphviz import Digraph
from chemical_flow_analysis.auxiliary import (estimating_mass,
                                              emission_factor,
                                              method_of_mements,
                                              non_zero_output_streams)
import pandas as pd


def building_pcu_black_box(chemical, Solubility, Tf, Tb, CoU, flow_influents, waste_stream, PCU_code, efficiency, dir_path, Metal_indicator='NO'):
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
    Carrier_for_PCU = Allocation.loc[Allocation['Code 2004 and prior'] == PCU_code, 'Output carrier'].iloc[0]
    Objective = Allocation.loc[Allocation['Code 2004 and prior'] == PCU_code, 'Objective'].iloc[0]
    Type_of_WM = Allocation.loc[Allocation['Code 2004 and prior'] == PCU_code, 'Type of waste management'].iloc[0]
    Type_of_treatment = Allocation.loc[Allocation['Code 2004 and prior'] == PCU_code, 'If it is treatment, what kind of?'].iloc[0]
    try:
        # Checking flows
        mu_flow_in, theta_2_flow_in = method_of_mements(flow_influents)
        ## Emission factors
        Releases['Emission_factor'] = Releases.apply(lambda x:
                                                    emission_factor(
                                                    x['Flow to compartment'],
                                                    x['Maximum amount on-site'],
                                                    x['Total waste'],
                                                    x['Total release']),
                                                    axis = 1)

        Fugitive_emissions_factors = pd.DataFrame(Releases.loc[Releases['Compartment'] == 'Fugitive air emission', 'Emission_factor'])
        Fugitive_emissions_factors.reset_index(inplace = True)
        n_rows = Fugitive_emissions_factors.shape[0]
        Stack_emission_factors = pd.DataFrame(Releases.loc[Releases['Compartment'] == 'Stack air emission', 'Emission_factor'])
        Stack_emission_factors.reset_index(inplace = True)
        Discharge_emission_factors = pd.DataFrame(Releases.loc[Releases['Compartment'] == 'On-site surface water', 'Emission_factor'])
        Discharge_emission_factors.reset_index(inplace = True)
        Soil_emission_factors = pd.DataFrame(Releases.loc[Releases['Compartment'] == 'On-site soil', 'Emission_factor'])
        Soil_emission_factors.reset_index(inplace = True)
        del Releases
        Flows_out = pd.DataFrame()
        for i in range(n_rows):
            flow_influent = estimating_mass(mu_flow_in, theta_2_flow_in)
            Fugitive_emission_factor = Fugitive_emissions_factors.loc[i, 'Emission_factor']
            Stack_emission_factor = Stack_emission_factors.loc[i, 'Emission_factor']
            Discharge_emission_factor = Discharge_emission_factors.loc[i, 'Emission_factor']
            Soil_emission_factor = Soil_emission_factors.loc[i, 'Emission_factor']
            ## Carrier
            if Objective == 'Removal':
                Carrier_flow = (1 - Fugitive_emission_factor)*flow_influent*efficiency/100
                Destroyed_converted_degradated_flow = 0.0
                # Carrier analysis
                if '/' in Carrier_for_PCU:
                    if Carrier_for_PCU == 'L/W':
                        if PCU_code in ['A03', 'A04', 'H103']:
                            if Solubility >= 1000: # Soluble in water if 1,000-10,000 and very solubla if >= 10,000mg/L
                                Carrier = 'W'
                            else:
                                Carrier = 'L'
                        elif PCU_code == 'P31':
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
                    Carrier = Carrier_for_PCU
            elif Objective == 'Reclamation':
                Carrier_flow = (1 - Fugitive_emission_factor)*flow_influent*efficiency/100
                Destroyed_converted_degradated_flow = 0.0
                if PCU_code in ['R11', 'R12', 'R13','R14', 'H20']: # Solvent recovery
                    Carrier = 'L'
                elif PCU_code in ['R21', 'R22', 'R23', 'R24', 'R26', 'R27', 'R28', 'R29', 'H10']: # Metal recovery
                    Carrier = 'S'
                else:
                    if Metal_indicator == 'YES':
                        Carrier = 'S'
                    else:
                        Carrier = 'W'
            elif (Objective == 'Safer disposal (RCRA Subtitle C)') | (Objective == 'Controlling flow rate and composition'):
                Carrier_flow = 0.0
                Carrier = None
                Destroyed_converted_degradated_flow = 0.0
            else:
                Carrier_flow = 0.0
                Carrier = None
                Destroyed_converted_degradated_flow = (1 - Fugitive_emission_factor)*flow_influent*efficiency/100
            if efficiency == 0.0:
                Carrier = None
            F_out_hidden = (1 - Fugitive_emission_factor)*(1 - efficiency/100)*flow_influent
            ## By product and Remanent
            if (Type_of_WM == 'Energy recovery') | (Type_of_treatment == 'Incineration'):
                if PCU_code == 'A01':
                    if Tb <= 60: # deg C chemical is assessed as a gas or liquid (Knock-out Drum)
                        Remanent_flow = F_out_hidden
                        Remanent = 'A'
                        By_product_flow = 0.0
                        By_product = None
                    else:
                        Remanent_flow = Stack_emission_factor*F_out_hidden
                        Remanent = 'A'
                        By_product_flow = (1 - Stack_emission_factor)*F_out_hidden
                        if Tf <= 850: # Assesing normal operating temperature of flares
                            By_product = 'L'
                        else:
                            By_product = 'S'
                else:
                    # Assesing normal operating temperature of these devices
                    if (Tb <= 980) | (Tf <= 980): # deg C chemical is assessed as a gas or liquid
                        Remanent_flow = F_out_hidden
                        Remanent = 'A'
                        By_product_flow = 0.0
                        By_product = None
                    else: # deg C Chemical is assessed as a solid
                        Remanent_flow = Stack_emission_factor*F_out_hidden
                        Remanent = 'A'
                        By_product_flow = (1 - Stack_emission_factor)*F_out_hidden
                        By_product = 'S'
            else:
                if PCU_code in ['H20', 'R11', 'R12', 'R13', 'R14', 'R40', 'H39']:
                    if waste_stream == 'A':
                        Remanent_flow = Stack_emission_factor*F_out_hidden
                        Remanent = 'A'
                    else:
                        Remanent_flow = (1 - Discharge_emission_factor)*F_out_hidden
                        Remanent = waste_stream
                    By_product_flow = F_out_hidden - Remanent_flow
                    By_product = 'W'
                elif PCU_code in ['H10', 'R21', 'R22', 'R23', 'R24', 'R26']:
                    By_product_flow = Discharge_emission_factor*F_out_hidden
                    By_product = 'W'
                    Remanent_flow =  F_out_hidden - By_product_flow
                    Remanent = waste_stream
                elif PCU_code in ['R27', 'R28', 'R29']:
                    By_product_flow = Soil_emission_factor*F_out_hidden
                    By_product = 'S'
                    Remanent_flow =  F_out_hidden - By_product_flow
                    Remanent = waste_stream
                else:
                    Remanent_flow =  F_out_hidden
                    Remanent = waste_stream
                    By_product_flow = 0.0
                    By_product = None
            Flows_out_aux = df_result_aux = pd.DataFrame({'Destroyed/converted/degradated': [Destroyed_converted_degradated_flow],
                                                          'Carrier': [[Carrier_flow, Carrier]],
                                                          'Remanent': [[Remanent_flow, Remanent]],
                                                          'By-product': [[By_product_flow, By_product]],
                                                          'Fugitive emission': [Fugitive_emission_factor*flow_influent]})
            Flows_out = pd.concat([Flows_out, Flows_out_aux], axis = 0,
                            sort = False,
                            ignore_index = True)
        Cols = ['PCU', 'Chemical', 'Type of stream', 'Mean quantity [kg/yr]', 'Std [kg/yr]', 'Number of data', 'Phase']
        df_result = pd.DataFrame(columns = Cols)
        for col in Flows_out.columns:
            Result = {'Type of stream': [col]}
            if col in ['Destroyed/converted/degradated', 'Fugitive emission']:
                Val =  Flows_out[col]
                if col == 'Destroyed/converted/degradated':
                    Phase = 'NA'
                else:
                    Phase = 'A'
            else:
                Val = pd.Series([val[0] for val in Flows_out[col]])
                Phase = Flows_out[col][0][1]
            Mean = Val.mean()
            Number =  Val.shape[0]
            Std =  100*Val.std()
            Result.update({'PCU': [PCU_code],
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
        return df_result, pd.Series([val[0] for val in Flows_out['Remanent']])
    except ZeroDivisionError:
        df_result = pd.DataFrame({'PCU': [PCU_code],
                                  'Chemical': [chemical],
                                  'Mean quantity [kg/yr]': ['NA'],
                                  'Std [kg/yr]': ['NA'],
                                  'Number of data': ['NA'],
                                  'Phase': ['NA']})
        return df_result,  pd.Series([0.0])



def picture(df_for_stream, stream, dir_path):
    # Organizing information
    n_pcus = df_for_stream['PCU'].nunique()
    PCU_position = {row['Position']: row['PCU'] for idx, row in df_for_stream[['PCU', 'Position']].drop_duplicates(keep = 'first').iterrows()}
    Chemical_PCU = {value: df_for_stream.loc[df_for_stream['PCU'] == value, 'Chemical'].tolist() for value in PCU_position.values()}
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
    PCU = Digraph()
    PCU.graph_attr['rankdir'] = 'LR'
    PCU.graph_attr['ranksep'] = '1'
    # Creating nodes
    ## Input
    PCU.node(str(0), 'Influent stream', fontname = 'Times New Roman Bold', style = 'filled', color = 'lightsalmon', shape = 'ellipse')
    ## PCUs
    Flows = pd.DataFrame()
    for position in PCU_position.keys():
        PCU_code = PCU_position[position]
        PCU.node(str(position), PCU_code, fontname = 'Times New Roman Bold', shape = 'square')
        for chemical in Chemicals_in_stream:
            if chemical in Chemical_PCU[PCU_code]:
                efficiency = Eff_chemical[chemical]
            else:
                efficiency = 0.0
            CoU = CoU_chemical[chemical]
            Solubility = Solubility_chemical[chemical]
            Tf = Tf_chemical[chemical]
            Tb = Tb_chemical[chemical]
            flow_influent = Flow_in[chemical]
            Result, flow_effluent = building_pcu_black_box(chemical, Solubility, Tf, Tb, CoU, flow_influent, waste_stream, PCU_code, efficiency, dir_path)
            Result['Position'] = position
            Flows = pd.concat([Flows, Result], axis = 0,
                        sort = False,
                        ignore_index = True)
            Flow_in[chemical] = flow_effluent
    Flows['Stream'] = stream
    ## Adding uncertainty in Influent stream
    Max_position = Flows['Position'].max()
    Flows.loc[(Flows['Type of stream'] == 'Remanent') &\
              (Flows['Position'] == Max_position),\
              'Type of stream'] = 'Effluent'
    Chemicals = Flows['Chemical'].unique().tolist()
    for chem in Chemicals:
        df_chem = Flows.loc[Flows['Chemical'] == chem]
        df_chem = df_chem.loc[df_chem['Type of stream'] != 'Remanent']
        df_chem['PCU'] = None
        df_chem['Type of stream'] = 'Influent'
        df_chem['Position'] = 0
        df_chem['Mean quantity [kg/yr]'] = df_chem['Mean quantity [kg/yr]'].sum()
        df_chem['Std [kg/yr]'] = df_chem['Std [kg/yr]'].sum()
        df_chem['Phase'] = waste_stream
        df_chem.drop_duplicates(keep = 'first', inplace = True)
        df_chem['The influent is in the interval?'] = df_chem.apply(lambda x: 'Yes' if (Flow_in_aux[chem] <= x['Mean quantity [kg/yr]'] + x['Std [kg/yr]'] and\
                                                                                        Flow_in_aux[chem] >= x['Mean quantity [kg/yr]'] - x['Std [kg/yr]'])
                                                                                    else 'No',
                                                                    axis = 1)
        Flows = pd.concat([df_chem, Flows], axis = 0,
                    sort = False,
                    ignore_index = True)
    Flows.sort_values(by = ['Position', 'Chemical'],
                      ascending = [True, True],
                      inplace = True)
    Flows['Mean quantity [kg/yr]'] = Flows['Mean quantity [kg/yr]'].round(4)
    Flows['Std [kg/yr]'] = Flows['Std [kg/yr]'].round(4)
    ## Output
    PCU.node(str(n_pcus + 1), 'Effluent stream', fontname = 'Times New Roman Bold', style = 'filled', color = 'lightblue', shape = 'ellipse')
    # Searching non zero streams for PCU
    Output_streams = dict()
    for Pos, PCU_code in PCU_position.items():
        PCU_flows = Flows.loc[Flows['PCU'] == PCU_code]
        Output_streams.update({Pos: non_zero_output_streams(PCU_flows)})
    # Creating edges
    n_streams = 0
    n_node = n_pcus + 1
    for i in range(n_pcus + 1):
        if i in Output_streams.keys():
            Flow_graph = Digraph()
            Flow_graph.graph_attr['rank'] = 'same'
            Flow_graph.graph_attr['ranksep'] = '2'
            for strm in Output_streams[i]:
                n_node = n_node + 1
                n_streams = n_streams + 1
                if strm in ['Fugitive emission', 'By-product']:
                    Flow_graph.node(str(n_node), strm, fontname = 'Times New Roman Bold', style = 'filled', color = 'lightgray', shape = 'ellipse')
                    Flow_graph.edge(str(n_node), str(i), label = str(n_streams),
                                    len = '1', arrowhead = 'normal', arrowsize = '1',
                                    dir = 'back', fontname = 'Comic Sans MS Bold')
                elif strm in ['Destroyed/converted/degradated', 'Carrier']:
                    Flow_graph.node(str(n_node), strm,
                                    fontname='Times New Roman Bold',
                                    style='filled', color='#b2df8a',
                                    shape='ellipse')
                    Flow_graph.edge(str(i), str(n_node), label=str(n_streams),
                                    len='1', arrowhead='normal', arrowsize='1',
                                    fontname='Comic Sans MS Bold')
            PCU.subgraph(Flow_graph)
        n_streams = n_streams + 1
        PCU.edge(str(i), str(i + 1), label=str(n_streams),
                 len='0.5', arrowhead='normal', arrowsize='1',
                 fontname='Comic Sans MS Bold')
    PCU.view(filename=f'{dir_path}/chemical_flow_analysis/Pollution_abatement_for_stream_{stream}')
    return Flows
