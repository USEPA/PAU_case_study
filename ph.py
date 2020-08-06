# -*- coding: utf-8 -*-
# !/usr/bin/env python

# Importing libraries
import math


def corrosiveness_score(pKa, pKb, w_w, pho, MW):
    pH = ph_calculator(pKa, pKb, w_w, pho, MW)
    if (2 <= pH) and (pH >= 12.5):
        return 4
    elif ((2 < pH) and (pH <= 4)) and ((9 <= pH) and (pH < 12.5)):
        return 3
    elif ((4 < pH) and (pH <= 5.5)) and ((8 <= pH) and (pH < 9)):
        return 2
    elif ((5.5 < pH) and (pH < 7)) and ((7 < pH) and (pH < 8)):
        return 1
    elif pH == 7:
        return 0


def ph_calculator(pKa, pKb, w_w, pho, MW):
    # density of solution
    pho_sln = f_pho_sln(w_w, pho)
    molarity = f_molarity(w_w, pho_sln, MW)
    if pKa:
        # HA <-> A- + H+
        Ka = 10**-pKa
        a = 1
        b = Ka
        c = -Ka*molarity
        pH = - math.log10(quadratic_formula(a, b, c, molarity))
        if not pH:
            pH = 7
    elif pKb:
        # B + H2O ⇌ HB+ + OH−
        Kb = 10**-pKb
        a = 1
        b = Kb
        c = -Kb*molarity
        pOH = - math.log10(quadratic_formula(a, b, c, molarity))
        if not pOH:
            pOH = 7
        pH = 14 - pOH
    else:
        return 7
    return pH


def f_pho_sln(w_w, pho):
    # density in g/cm3
    # % w/w
    # pho_sln in g/L
    pho_water = 0.997
    pho_sln = 10*((1 - w_w)*pho_water + w_w*pho)
    return pho_sln


def f_molarity(w_w, pho_sln, MW):
    # pho_sln in g/L
    molarity = (w_w/100)*(1/MW)*pho_sln
    return molarity


def quadratic_formula(a, b, c, molarity):
    solution_1 = (-b - (b**2 - 4*a*c)**0.5)/(2*a)
    solution_2 = (-b + (b**2 - 4*a*c)**0.5)/(2*a)
    if solution_2 == solution_1:
        return solution_1
    elif (solution_2 < 0) and (solution_1 < 0):
        return None
    elif solution_2 < 0:
        return solution_1
    elif solution_1 < 0:
        return solution_2
    elif isinstance(solution_1, complex):
        return None
    else:
        return max([s for s in [solution_1, solution_2] if s <= molarity])
