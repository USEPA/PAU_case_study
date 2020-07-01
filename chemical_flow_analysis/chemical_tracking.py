# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Importing libraries
from graphviz import Digraph
from chemical_flow_analysis.auxiliary import *

def building_PCU_black_box(chemical, flow, waste_stream, PCU, dir_path):
    Releases = pd.read_csv(dir_path + '/chemical_flow_analysis/tri_releases/TRI_releases.csv')
    Releases = Releases.loc[Releases['CAS NUMBER'] == chemical]
    # Allocation
    Allocation = pd.read_csv(dir_path + '/Methods_TRI.csv')
    Allocation = Allocation.where(pd.notnull(Allocation), None)
    Carrier = Allocation.loc[Allocation['Code 2004 and prior'] == PCU, 'Output carrier'].iloc[0]
    # Checking flows
    ## Emission factors
    ## Carrier
    if Objective == 'Removal':
        Carrier_flow = int((flow - Fugitive_emission)*Efficiency/100)
        Destroyed_converted_degradated_flow = 0
        # Carrier analysis
        if '/' in Carrier:
            if Carrier == 'L/W':
                if PCU in ['A03', 'A04', 'H103']:
                    if Solubility >= 1000: # Soluble in water if 1,000-10,000 and very solubla if >= 10,000mg/L
                        Carrier = 'W'
                    else:
                        Carrier = 'L'
                elif PCU == 'P31':
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
    elif Objective == 'Reclamation':
        Destroyed_converted_degradated_flow = 0
        Carrier_flow = int((flow - Fugitive_emission)*Efficiency/100)
        #if PCU in ['R11', 'R12', 'R13','R14', 'R19', 'H20']: # Solvent recovery
        #    Carrier = 'L'
        #elif PCU in ['R21', 'R22', 'R23', 'R24', 'R26', 'R27', 'R28', 'R29', 'R30', 'H10']: # Metal recovery
        #    Carrier = 'S'
        #else:
        #    if Metal_indicator == 'YES':
        #        Carrier = 'S'
        #    else:
        #        Carrier = 'W'
    elif (Objective == 'Safer disposal (RCRA Subtitle C)') | (Objective == 'Controlling flow rate and composition'):
        Carrier_flow = 0
        Destroyed_converted_degradated_flow = 0
    else:
        Carrier_flow = 0
        Destroyed_converted_degradated_flow = int((flow - Fugitive_emission)*Efficiency/100)
    if Efficiency == 0.0:
        Carrier = None
    ## By product and Remanent
    if (Type_of_WM == 'Energy recovery') | (Type_of_treatment == 'Incineration'):
        if PCU == 'A01':
            if Tb <= 60: # deg C chemical is assessed as a gas or liquid (Knock-out Drum)
                Remanent_flow = flow - (Destroyed_converted_degradated_flow + Fugitive_emission)
                Remanent = 'A'
                By_product_flow = 0
                By_product = None
            else:
                Remanent_flow = int((flow - (Destroyed_converted_degradated_flow + Fugitive_emission))*Total_stack_emission/Production_volume)
                Remanent = 'A'
                By_product_flow = int(flow - (Remanent_flow + Carrier_flow + Fugitive_emission))
                if Tf <= 850: # Assesing normal operating temperature of flares
                    By_product = 'L'
                else:
                    By_product = 'S'
        else:
            # Assesing normal operating temperature of these devices
            if (Tb <= 980) | (Tf <= 980): # deg C chemical is assessed as a gas or liquid
                Remanent_flow = int(flow - (Destroyed_converted_degradated_flow + Fugitive_emission))
                Remanent = 'A'
                By_product_flow = 0
                By_product = None
            else: # deg C Chemical is assessed as a solid
                Remanent_flow = int((flow - (Destroyed_converted_degradated_flow + Fugitive_emission))*Total_stack_emission/Production_volume)
                Remanent = 'A'
                By_product_flow = int(flow - (Remanent_flow + Carrier_flow + Fugitive_emission))
                By_product = 'S'
    else:
        if PCU in ['H20', 'R11', 'R12', 'R13', 'R14', 'R19']:
            if waste_stream == 'A':
                Remanent_flow = int((flow - (Carrier_flow + Fugitive_emission))*Total_stack_emission/Production_volume)
                Remanent = 'A'
            else:
                Remanent_flow = 0
                Remanent = None
            By_product_flow = int(flow - (Remanent_flow + Carrier_flow + Fugitive_emission))
            By_product = 'S'
        else:
            if (Objective == 'Safer disposal (RCRA Subtitle C)') | (Objective == 'Controlling flow rate and composition'):
                Remanent_flow = int(flow - Fugitive_emission)
            else:
                Remanent_flow = int(flow - (Carrier_flow + Fugitive_emission))
            Remanent = waste_stream
            By_product_flow = 0
            By_product = None
    # Returning values
    return Fugitive_emission, Carrier_flow, Carrier, Destroyed_converted_degradated_flow, Remanent_flow, Remanent, By_product_flow, By_product


def picture(df_for_stream, stream, dir_path):
    # Creating graph
    PCU = Digraph()
    PCU.graph_attr['rankdir'] = 'LR'
    PCU.graph_attr['ranksep'] = '2'
    # Creating nodes
    n_pcus = df_for_stream['PCU'].nunque()
    PCU_position = {row['Position']: row['PCU'] for idx, row in df_for_stream[['PCU', 'Position']].drop_duplicates(keep = 'first').iterrows()}
    Chemical_PCU = {value: df_for_stream.loc[df_for_stream['PCU'] == value, 'Chemical'].tolist() for value in PCU_position.values()}
    Chemicals_in_stream = df_for_stream['Chemical'].unique().tolist()
    Flows = {}
    waste_stream = df_for_stream['Type of waste'].iloc[0]
    ## Input
    PCU.node(str(0), 'Influent stream', style = 'filled', color = 'lightsalmon', shape = 'ellipse')
    ## PCUs
    for position in PCU_position.keys():
        PCU_code = PCU_position[position]
        PCU.node(str(position), PCU_code, shape = 'square')
        for chemical in Chemicals_in_stream:
            Result = building_PCU_black_box(chemical, flow, waste_stream, PCU_code, dir_path)
    ## Output
    PCU.node(str(n_pcus + 1), 'Effluent stream', style = 'filled', color = 'lightblue', shape = 'ellipse')
    # Creating edges
    for i in range(n_pcus + 1):
        PCU.edge(str(i), str(i + 1),
                len = '1', arrowhead = 'normal', arrowsize = '1')
    PCU.view(filename = dir_path + '/chemical_flow_analysis/Pollution_abatement_for_stream_{}'.format(stream))
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
