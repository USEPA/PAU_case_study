# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Importing libraries
import numpy as np

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


def annual_change(max_annual_change, Total_releases_from_facility, Total_waste_at_facility):
    min_value = max([-max_annual_change, Total_releases_from_facility - Total_waste_at_facility])
    max_value = max_annual_change
    return np.random.uniform(min_value, max_value)


def cal(keys, m, w, r):
    Max_onsite = maximum_on_site(m)
    results = {key: 1/(annual_change(Max_onsite, r, w) + w) for key in keys}
    return results
