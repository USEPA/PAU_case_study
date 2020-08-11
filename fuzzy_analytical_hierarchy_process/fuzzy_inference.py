# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Importing libraries
import numpy as np
import pandas as pd
import os
import networkx as nx

def comparison_matrix(df, cols):
    Some_criteri = ['WHM importance', 'T Flammability',
                    'T Instability', 'T Corrosiveness']
    # Based on WMH
    m_criteria = 0
    n = df.shape[0]
    for col in cols:
        m_criteria = m_criteria + 1
        val = df[col].tolist()
        N_aux = np.empty((n, n))
        N_aux[:] = np.nan
        for i in range(n):
            for j in range(i, n):
                diff = val[i] - val[j]
                if col == 'Position database':
                    diff = round(4*((diff - (n - 1))/(n - 1) + 1))
                elif not col in Some_criteri:
                    diff = round(4*(diff - 1) + 4)
                N_aux[i][j] = diff
        if m_criteria == 1:
            N = N_aux
        else:
            N =  np.concatenate((N, N_aux), axis = 0)
    return N


def analysing_intersection(df, n_concerning_constituents):
    PCUs = df['PCU'].unique().tolist()
    df['P-based_on_intersection'] = 1/n_concerning_constituents
    for PCU in PCUs:
        n = len(df.loc[df['PCU'] == PCU, 'Chemical'].tolist())
        if n != 1:
            df.loc[df['PCU'] == PCU, 'P-based_on_intersection'] = n/n_concerning_constituents
        else:
            chem = df.loc[df['PCU'] == PCU, 'Chemical'].values[0]
            n_pcu = len(df.loc[df['Chemical'] == chem, 'PCU'].tolist())
            n_pcu_false = len(df.loc[(df['Chemical'] == chem) & (~ df['Intersection']), 'PCU'].tolist())
            df.loc[(df['PCU'] == PCU) & (df['Chemical'] == chem), 'P-based_on_intersection'] =  n_pcu_false/n_pcu
    return df

def analysing_position_based_on_PCU_database(df):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    df_order = pd.read_csv(dir_path + '/PCU_positions/Positions.csv')
    PCUs = df['PCU'].unique().tolist()
    df_order = df_order.loc[(df_order['First'].isin(PCUs)) & (df_order['Second'].isin(PCUs))]
    if not df_order.empty:
        G = nx.from_pandas_edgelist(df_order, 'First', 'Second',
                                     create_using = nx.DiGraph())
        if nx.is_directed_acyclic_graph(G):
            Order = list(reversed(list(nx.topological_sort(G))))
            Position = {val: idx for idx, val in enumerate(Order)}
            df['Position database'] = df['PCU'].apply(lambda x: Position[x])
        else:
            df['Position database'] = 1
    else:
        df['Position database'] = 1
    return df


def pairwise_comparison(df, objective):
    # In the current work for on-site pollution abatement, only treatment, energy recovery,
    # and recycling are included. However, we added the others to clarify
    WHM_importance = {'Source reduction': 5,
                      'Recycling': 4,
                      'Energy recovery': 3,
                      'Treatment': 2,
                      'Disposal': 1}
    Based_on_objective = {'Controlling flow rate and composition': 1,
                          'Removal': 2,
                          'Reclamation': 3,
                          'Conversion': 4,
                          'Degradation': 5,
                          'Destruction': 6,
                          'Energy': 6} # The order is inverse
    # objective:
    ## pcu: for selecting PCU
    ## seq: for selecting sequence
    if objective == 'pcu':
        df['WHM importance'] = df['Type_of_waste_management'].apply(lambda x: WHM_importance[x])
        n_concerning_constituents = df['Chemical'].nunique()
        if n_concerning_constituents == 1:
            n_probable_pcu = len(df['PCU'].tolist())
            if n_probable_pcu == 1:
                # Only a concerning chemical in the stream
                df['Selected'] = 'Yes'
                df.drop(columns = ['WHM importance'], inplace = True)
            else:
                # Call fuzzy function to make decision based on WMH, PCU-probability, and WMH probability
                criteria_on_columns = ['WHM importance', 'PCU-probability', 'Type_of_waste_management-probability']
                df = fahp(n_probable_pcu, criteria_on_columns, df)
                df['Selected'] = 'No'
                df.loc[df['Weight'].idxmax(), 'Selected'] = 'Yes'
                df.drop(columns = ['Weight', 'WHM importance'], inplace = True)
        else:
            df['Intersection'] = df.duplicated(subset = ['PCU'], keep = False)
            # Analyzing importance based on intersections
            df = analysing_intersection(df, n_concerning_constituents)
            chemicals = df['Chemical'].unique().tolist()
            criteria_on_columns = ['WHM importance', 'PCU-probability', 'Type_of_waste_management-probability', 'P-based_on_intersection']
            df_aux = df.copy()
            df = pd.DataFrame()
            for chemical in chemicals:
                df_chem = df_aux.loc[df_aux['Chemical'] == chemical]
                if df_chem.shape[0] == 1:
                    df_chem = pairwise_comparison(df_chem, objective = 'pcu')
                    df_chem.drop(columns = ['P-based_on_intersection', 'Intersection'], inplace = True)
                else:
                    n_probable_pcu =  len(df_chem['PCU'].tolist())
                    df_chem = fahp(n_probable_pcu, criteria_on_columns, df_chem)
                    df_chem['Selected'] = 'No'
                    df_chem.loc[df_chem['Weight'].idxmax(), 'Selected'] = 'Yes'
                    df_chem.drop(columns = ['Weight', 'WHM importance', 'P-based_on_intersection', 'Intersection'], inplace = True)
                df = pd.concat([df, df_chem], axis = 0, sort = False, ignore_index = True)
            del df_chem
    else:
        n_pcus = df['PCU'].nunique()
        if n_pcus == 1:
            df['Position'] = 1
        else:
            # Based on the objective
            df['Sorted objective'] = df['Objective'].apply(lambda x: Based_on_objective[x])
            df.sort_values(by = ['Sorted objective', 'PCU'], inplace = True, ascending = [True, False])
            Position_dict = {val: idx + 1 for idx, val in enumerate(df['PCU'].unique())}
            df['Position'] = df['PCU'].apply(lambda x: Position_dict[x])
            # Checking removal treatments
            df_removal = df.loc[df['Objective'] == 'Removal']
            df = df.loc[df['Objective'] != 'Removal']
            if not df_removal.empty:
                if df_removal.shape[0] != 1:
                    min = df_removal['Position'].min()
                    # Summing removal flows and searching highest flammability, instability, and corrosiveness
                    Aux_dict = {'Chemical flow': ['T Chemical flow', 'sum'],
                                'Flammability': ['T Flammability', 'max'],
                                'Instability': ['T Instability', 'max'],
                                'Corrosiveness': ['T Corrosiveness', 'max']}
                    for key, val in Aux_dict.items():
                        df_removal[val[0]] = df_removal.groupby('PCU',
                                                                as_index = False)\
                                                                [key].transform(val[1])
                        if key == 'Chemical flow':
                            df_removal[val[0]] = df_removal[val[0]]/df_removal[val[0]].max()
                    cols_for_fuzziness = ['PCU', 'T Flammability', 'T Instability',
                                          'T Corrosiveness', 'T Chemical flow']
                    df_removal_reduce = df_removal[cols_for_fuzziness]
                    df_removal.drop(columns = ['T Flammability', 'T Instability',
                                               'T Corrosiveness', 'T Chemical flow',
                                               'Position'],
                                    inplace = True)
                    df_removal_reduce.drop_duplicates(keep = 'first', inplace = True)
                    df_removal_reduce = analysing_position_based_on_PCU_database(df_removal_reduce)
                    criteria_on_columns = ['T Flammability', 'T Instability',
                                           'T Corrosiveness', 'Position database',
                                           'T Chemical flow']
                    n_probable_pcu = df_removal_reduce['PCU'].nunique()
                    df_removal_reduce = fahp(n_probable_pcu, criteria_on_columns, df_removal_reduce)
                    df_removal_reduce.sort_values(by = ['Weight'], inplace = True, ascending = False)
                    df_removal_reduce['Position'] = pd.Series(np.arange(min, df_removal_reduce.shape[0] + min))
                    df_removal_reduce.drop(columns = criteria_on_columns + ['Weight'],
                                          inplace = True)
                    df_removal = pd.merge(df_removal, df_removal_reduce, on = 'PCU', how = 'inner')
                    del df_removal_reduce
                    df = pd.concat([df, df_removal], axis = 0, sort = False, ignore_index = True)
                else:
                    df = pd.concat([df, df_removal], axis = 0, sort = False, ignore_index = True)
            df.drop(columns = ['Sorted objective'],
                    inplace = True)
    return df


def fahp(n, cols_criteri, df):
    m_criteria = len(cols_criteri)
    N = comparison_matrix(df, cols_criteri)
    # Definition of variables
    W = np.zeros((1, n)) # initial value of weights vector (desired output)
    w = np.zeros((m_criteria, n)) # initial value of weithts matrix
    phi = np.zeros((m_criteria, 1)) # diversification degree
    eps = np.zeros((m_criteria, 1)) # entropy
    theta = np.zeros((1, m_criteria)) # criteria' uncertainty degrees
    sumphi = 0 # Initial value of diversification degree
    # Triangular fuzzy numbers (TFN)
    # n is number of PCUs
    # TFN is a vector with n segments and each one has 3 different numbers
    # The segments are the linguistic scales
    # The 3 differente numbers in the segments are a triangular fuzzy number (l,m,u)
    # the first segment: equal importance; the second one: moderate importance of one over another;
    # the third one: strong importance of one over another; the fourth one: very strong importance of one over another
    # the fifth one: Absolute importance of one over another
    TFN = np.array([1, 1, 1 ,2/3, 1, 3/2, 3/2, 2, 5/2, 5/2, 3, 7/2, 7/2, 4, 9/2])
    for k in range(1, m_criteria + 1):
        a = np.zeros((1, n*n*3))  # Comparison matrix (In this case is a vector because of computational memory
        for i in range(k*n-(n-1), k*n + 1):
            for j in range(i-n*(k-1), n + 1):
                # This is the position of the third element of the segment for
	            # a*(i,j) (upper triangular part)
                jj = 3*(n*((i-n*(k-1))-1)+j)
                # This is the position of the thrid element of the segment for
	            # a*(j,i) (lower triangular part)
                jjj = 3*(n*(j-1) + i-n*(k-1))
                if N[i - 1][j - 1] == -4:
                    a[0][jjj-3:jjj] =  TFN[12:15]
                    a[0][jj-3:jj] = np.array([TFN[14]**-1, TFN[13]**-1, TFN[12]**-1])
                elif N[i - 1][j - 1] == -3:
                    a[0][jjj-3:jjj] =  TFN[9:12]
                    a[0][jj-3:jj] = np.array([TFN[11]**-1, TFN[10]**-1, TFN[9]**-1])
                elif N[i - 1][j - 1] == -2:
                    a[0][jjj-3:jjj] = TFN[6:9]
                    a[0][jj-3:jj] = np.array([TFN[8]**-1, TFN[7]**-1, TFN[6]**-1])
                elif N[i - 1][j - 1] == -1:
                    a[0][jjj-3:jjj] =  TFN[3:6]
                    a[0][jj-3:jj] = np.array([TFN[5]**-1, TFN[4]**-1, TFN[3]**-1])
                elif N[i - 1][j - 1] == 0:
                    a[0][jj-3:jj] =  TFN[0:3]
                    a[0][jjj-3:jjj] = TFN[0:3]
                elif N[i - 1][j - 1] == 1:
                    a[0][jj-3:jj] =  TFN[3:6]
                    a[0][jjj-3:jjj] = np.array([TFN[5]**-1, TFN[4]**-1, TFN[3]**-1])
                elif N[i - 1][j - 1] == 2:
                    a[0][jj-3:jj] = TFN[6:9]
                    a[0][jjj-3:jjj] = np.array([TFN[8]**-1, TFN[7]**-1, TFN[6]**-1])
                elif N[i - 1][j - 1] == 3:
                    a[0][jj-3:jj] =  TFN[9:12]
                    a[0][jjj-3:jjj] = np.array([TFN[11]**-1, TFN[10]**-1, TFN[9]**-1])
                elif N[i - 1][j - 1] == 4:
                    a[0][jj-3:jj] =  TFN[12:15]
                    a[0][jjj-3:jjj] = np.array([TFN[14]**-1, TFN[13]**-1, TFN[12]**-1])
        # (2) fuzzy synthetic extension
        A = np.zeros((n,3))
        B = np.zeros((1,3))
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                jj = 3*(n*(i-1)+j)
                A[i - 1][:] = A[i - 1][:] + a[0][jj-3:jj]
            B = B + A[i - 1][:]
        BB = np.array([B[0][2]**-1, B[0][1]**-1, B[0][0]**-1])
        S = A*BB
        # (3) Degree of possibility
        for i in range(n):
            V = np.zeros(n)
            for j in range(n):
                if S[i][1] >= S[j][1]:
                    V[j] = 1
                elif S[j][0] >= S[i][2]:
                    V[j] = 0
                else:
                    V[j] = (S[j][0] - S[i][2])/((S[i][1] - S[i][2]) - (S[j][1] - S[j][0]))
            w[k - 1][i] = np.min(V)
        # (4) Weight of each flow for a criterium
        w[k - 1][:] = (np.sum(w[k - 1][:])**-1)*w[k - 1][:]
        # (5) Criteria' uncertainty degrees
        for i in range(n):
            if w[k - 1][i] != 0:
                eps[k - 1][0] = eps[k - 1][0] - w[k - 1][i]*np.log(w[k - 1][i])
        eps[k - 1][0] = eps[k - 1][0]/np.log(n)
        phi[k - 1][0] = 1 + eps[k - 1]
        sumphi = sumphi + phi[k - 1][0]
    # (6) Final weight of all flows
    for i in range(n):
        for k in range(m_criteria):
            theta[0][k] = phi[k][0]/sumphi
            W[0][i] = W[0][i] + w[k][i]*theta[0][k]
    W = (np.sum(W)**-1)*W # Final weights
    W = W.T
    df = df.assign(Weight = W)
    return df
