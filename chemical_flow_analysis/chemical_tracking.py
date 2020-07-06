# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Importing libraries
from graphviz import Digraph
from chemical_flow_analysis.auxiliary import *
import pandas as pd

def building_pcu_black_box(chemical, Solubility, Tf, Tb, CoU, flow_influents, waste_stream, PCU_code, efficiency, dir_path):
    Releases = pd.read_csv(dir_path + '/chemical_flow_analysis/tri_releases/TRI_releases.csv')
    Releases = Releases.loc[Releases['CAS NUMBER'] == chemical]
    for key, val in CoU.items():
        Releases = Releases.loc[Releases[key] == val]
    Releases.drop(columns = ['CAS NUMBER', 'As a byproduct',
    	                      'As a manufactured impurity',
                              'As a process impurity'], inplace = True)
    Duplicated_columns = ['Maximum amount on-site',
    	                  'Total waste', 'Total release']
    Releases = Releases.loc[Releases.duplicated(subset = Duplicated_columns, keep = False)]
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
        Fugitive_emissions_factors = Releases.loc[Releases['Compartment'] == 'Indoor air', 'Emission_factor']
        n_rows = Fugitive_emissions_factors.shape[0]
        Stack_emission_factors = Releases.loc[Releases['Compartment'] == 'Outdoor air', 'Emission_factor']
        del Releases
        Flows_out = pd.DataFrame()
        for i in range(n_rows):
            flow_influent = estimating_mass(mu_flow_in, theta_2_flow_in)
            Fugitive_emission_factor = Fugitive_emissions_factors.iloc[i]
            Stack_emission_factor = Stack_emission_factors.iloc[i]
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
                if PCU_code in ['R11', 'R12', 'R13','R14', 'R19', 'H20']: # Solvent recovery
                    Carrier = 'L'
                elif PCU_code in ['R21', 'R22', 'R23', 'R24', 'R26', 'R27', 'R28', 'R29', 'R30', 'H10']: # Metal recovery
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
                if PCU_code in ['H20', 'R11', 'R12', 'R13', 'R14', 'R19']:
                    if waste_stream == 'A':
                        Remanent_flow = Stack_emission_factor*F_out_hidden
                        Remanent = 'A'
                    else:
                        Remanent_flow = 0.0
                        Remanent = None
                    By_product_flow = F_out_hidden - Remanent_flow
                    By_product = 'S'
                else:
                    Remanent_flow =  F_out_hidden
                    Remanent = waste_stream
                    By_product_flow = 0.0
                    By_product = None
            Flows_out_aux = df_result_aux = pd.DataFrame({'Destroyed_converted_degradated': [Destroyed_converted_degradated_flow],
                                                          'Carrier': [[Carrier_flow, Carrier]],
                                                          'Remanent': [[Remanent_flow, Remanent]],
                                                          'By-product': [[By_product_flow, By_product]],
                                                          'Fugitive emission': [Fugitive_emission_factor*flow_influent]})
            Flows_out = pd.concat([Flows_out, Flows_out_aux], axis = 0,
                            sort = False,
                            ignore_index = True)
        Cols = ['PCU', 'Chemical', 'Type of stream', 'Mean quantity [kg/yr]', 'CV [%]', 'Number of data', 'Phase']
        df_result = pd.DataFrame(columns = Cols)
        for col in Flows_out.columns:
            Result = {'Type of stream': [col]}
            if col in ['Destroyed_converted_degradated', 'Fugitive emission']:
                Val =  Flows_out[col]
                if col == 'Destroyed_converted_degradated':
                    Phase = 'NA'
                else:
                    Phase = 'A'
            else:
                Val = pd.Series([val[0] for val in Flows_out[col]])
                Phase = Flows_out[col][0][1]
            Mean = Val.mean()
            Number =  Val.shape[0]
            if Mean != 0.0:
                CV = 100*Val.std()/Mean
            else:
                CV = 0.0
            Result.update({'PCU': [PCU_code],
                          'Chemical': [chemical],
                          'Mean quantity [kg/yr]': [Mean],
                          'CV [%]': [CV],
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
                                  'CV [%]': ['NA'],
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
    waste_stream = df_for_stream['Type of waste'].iloc[0]
    # Creating graph
    PCU = Digraph()
    PCU.graph_attr['rankdir'] = 'LR'
    PCU.graph_attr['ranksep'] = '2'
    # Creating nodes
    ## Input
    PCU.node(str(0), 'Influent stream', style = 'filled', color = 'lightsalmon', shape = 'ellipse')
    ## PCUs
    Flows = pd.DataFrame()
    for position in PCU_position.keys():
        PCU_code = PCU_position[position]
        PCU.node(str(position), PCU_code, shape = 'square')
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
    ## Output
    PCU.node(str(n_pcus + 1), 'Effluent stream', style = 'filled', color = 'lightblue', shape = 'ellipse')
    # Creating edges
    for i in range(n_pcus + 1):
        PCU.edge(str(i), str(i + 1),
                len = '1', arrowhead = 'normal', arrowsize = '1')
    PCU.view(filename = dir_path + '/chemical_flow_analysis/Pollution_abatement_for_stream_{}'.format(stream))
    return Flows
# def picture(F_id, Name_chem, waste_stream, flow, PCU_name, Objective, Efficiency, Fugitive_emission, Carrier_flow, Carrier, Destroyed_converted_degradated_flow, Remanent_flow, Remanent, By_product_flow, By_product):
#     # Creating graph
#     PCU = Digraph(filename = 'PCU_' + PCU_name + '_' + F_id +'_' + Name_chem)
#     PCU.graph_attr['rankdir'] = 'LR'
#     PCU.graph_attr['ranksep'] = '2'
#
#
#     # Creating nodes
#     PCU.node('B', PCU_name + '\n(Objective: ' + Objective + ')\n(' + str(Efficiency) + '%)', shape = 'square')
#     PCU.attr('node', shape = 'ellipse', style = 'filled', color = 'lightblue')
#     PCU.node('A', 'Influent', pos = '30,30!')
#     PCU.attr('node', style = 'filled', color = 'lightgray', shape = 'ellipse')
#     PCU.node('C', 'Fugitive emission')
#     PCU.node('E', 'Remanent')
#     PCU.node('F', 'By-product')
#     PCU.attr('node', style = 'filled', color = 'lightsalmon', shape = 'ellipse')
#     PCU.node('D', 'Carrier', pos = '30,30!')
#     PCU.node('G', 'Destroyed/Converted/Degradated')
#     PCU.attr(label = r'\n\nFlow allocation of ' + PCU_name.lower() + ':\n' + Name_chem)
#     PCU.attr(fontsize = '20')
#
#     PCU.attr(rank = 'same'; 'A', 'C')
#
#     ## Creating edges
#     PCU.edge('A', 'B', label = waste_stream + ': ' + f'{flow:n}' + ' kg/yr',
#              len = '2.00', arrowhead = 'vee', arrowsize = '2')
#     PCU.edge('B', 'E', label = Remanent + ': '+ f'{Remanent_flow:n}' + ' kg/yr' ,
#              len = '2.00', arrowhead = 'vee', arrowsize = '2')
#     PCU.edge('B', 'F', label = str(By_product) + ': '+ f'{By_product_flow:n}' + ' kg/yr',
#              len = '2.00', arrowhead = 'vee', arrowsize = '2')
#     PCU.edge('B', 'C', label = 'A: '+ f'{Fugitive_emission:n}' + ' kg/yr',
#              len = '2.00', style = 'dashed', tailport = 'n', headport = 'w',
#              arrowhead = 'vee', arrowsize = '2')
#     PCU.edge('B', 'D', label = str(Carrier) + ': '+ f'{Carrier_flow:n}' + ' kg/yr',
#              len = '2.00', arrowhead = 'vee', arrowsize = '2')
#     PCU.edge('B', 'G', label = f'{Destroyed_converted_degradated_flow:n}' + ' kg/yr',
#              len = '2.00', style = 'dashed', color = 'red', tailport = 's', headport = 'w',
#              arrowhead = 'vee', arrowsize = '2')
#     PCU.view()
