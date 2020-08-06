# -*- coding: utf-8 -*-
# !/usr/bin/env python

# Importing libraries
import numpy as np
import pandas as pd
from scipy.stats import lognorm


def maximum_on_site(x):
    if x == 1:
        return 0.453592*99
    elif x == 2:
        return 0.453592*999
    elif x == 3:
        return 0.453592*9999
    elif x == 4:
        return 0.453592*99999
    elif x == 5:
        return 0.453592*999999
    elif x == 6:
        return 0.453592*9999999
    elif x == 7:
        return 0.453592*49999999
    elif x == 8:
        return 0.453592*99999999
    elif x == 9:
        return 0.453592*499999999
    elif x == 10:
        return 0.453592*999999999
    elif x == 11:
        return 0.453592*10000000000
    elif x == 12:
        return 0.001*0.099
    elif x == 13:
        return 0.001*0.99
    elif x == 14:
        return 0.001*9.99
    elif x == 15:
        return 0.001*99
    elif x == 16:
        return 0.001*999
    elif x == 17:
        return 0.001*9999
    elif x == 18:
        return 0.001*99999
    elif x == 19:
        return 0.001*999999
    elif x == 20:
        return 0.001*100000000


def annual_change(Max_onsite, Total_releases_from_facility,
                  Total_waste_at_facility):
    min_value = max([-Max_onsite, Total_releases_from_facility
                    - Total_waste_at_facility])
    max_value = Max_onsite
    return np.random.uniform(min_value, max_value)


def emission_factor(Release_to_compartment, Max_onsite_code, Tota_waste,
                    Total_release):
    Max_onsite = maximum_on_site(Max_onsite_code)
    Emission_factor = Release_to_compartment/(annual_change(Max_onsite,
                                              Total_release, Tota_waste)
                                              + Tota_waste)
    return Emission_factor


def method_of_mements(Flow_vals):
    if (isinstance(Flow_vals, pd.Series)):
        Mean = Flow_vals.mean()
        StD = Flow_vals.std()
    else:
        Mean = Flow_vals
        # CV less than 1 (low variance)
        StD = Mean*0.01
    mu = np.log(Mean**2/(StD**2 + Mean**2)**0.5)
    theta_2 = np.log(StD**2/Mean**2 + 1)
    return mu, theta_2


def estimating_mass(mu, theta_2):
    try:
        Flow = lognorm.rvs(s=theta_2**0.5,
                           scale=np.exp(mu))
    except ValueError:
        Flow = lognorm.rvs(s=10**-9,
                           scale=np.exp(mu))
    return Flow


def non_zero_output_streams(PCU_flows):
    PCU_flows = PCU_flows.loc[~PCU_flows['Type of stream'].isin(['Effluent',
                              'Remanent'])]
    PCU_flows = PCU_flows[['Type of stream', 'Mean quantity [kg/yr]']]\
        .groupby('Type of stream', as_index=False).sum()
    List_of_non_zero = list()
    for idx, row in PCU_flows.iterrows():
        # Flow = row['Mean quantity [kg/yr]']
        Type = row['Type of stream']
        if row['Mean quantity [kg/yr]'] != 0.0:
            List_of_non_zero.append(Type)
    return List_of_non_zero
