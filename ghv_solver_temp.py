import pandas as pd
import numpy as np
from thermo import ChemicalConstantsPackage, PRMIX, CEOSLiquid, CEOSGas, FlashVLN
import thermo
import chemicals
from thermo.interaction_parameters import IPDB
import pint
import copy
import fluids
import config
import timeit
from scipy.optimize import newton
import correlations

start_time1 = timeit.default_timer()

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

ureg = pint.UnitRegistry()

settings = {
    "T_STANDARD": 288.70555,  # Temperature in Kelvin, 60F
    "P_STANDARD": 101325.0,  # Pressure in Pascal, 1 atm
    "R": 8.31446261815324,
}

config.update_config(settings)


def normalize_composition(comp_dict):
    """
    :param comp_dict: un-normalized dictionary of composition. {"CH4": 3, "C2H6", 6}
    :return: normalized dictionary of composition. {"CH4": 0.3333, "C2H6", 0.66666666}
    """
    total_comp = sum(comp_dict.values())
    if total_comp > 1:
        comp_dict = {k: v / total_comp for k, v in comp_dict.items()}
    return comp_dict


def check_properties_exists(constants):
    """
    :param constants: constants object of the thermo library
    Molecular weight and normal boiling points are minimum information needed to characterize a fluid, incase of any
    missing properties. For example, docosane is missing Heat of combustion (J/mol), but it can be correlated.
    """
    rhol_60Fs_mass = constants.rhol_60Fs_mass
    mws = constants.MWs
    Tbs = constants.Tbs
    names = constants.names
    Hcs = constants.Hcs

    # Todo: implement a method to add missing data.
    for rhol_60F_mass, mw, Tb, Hc, name in zip(rhol_60Fs_mass, mws, Tbs, Hcs, names):
        if rhol_60F_mass is None:
            raise ValueError(
                "Chemical name '%s' is recognized but missing a required data (rhol_60F_mass, liquid mass density at 60F)." % name)
        if mw is None:
            raise ValueError(
                "Chemical name '%s' is recognized but missing a required data (mw, molecular weight)." % name)
        if Tb is None:
            raise ValueError(
                "Chemical name '%s' is recognized but missing a required data (Tb, normal boiling temperature)." % name)
        if Hc is None:
            raise ValueError(
                "Chemical name '%s' is recognized but missing a required data (Hc, heat of combustion [J/mol])." % name)


def ideal_gas_molar_volume():
    """
    PV=nRT, where number of moles n=1. Rearranging -> V=RT/P
    R = 8.31446261815324 ((m^3-Pa)/(mol-K))
    T = 288.7056 K, 60F, standard temperature
    P = 101325 Pa, 1 atm, standard pressure
    :return: ideal gas molar volume in a standard condition (m^3/mol). 0.0236 m^3/mol at standard conditions for all compounds
    """
    return config.constants['R'] * config.constants['T_STANDARD'] / config.constants['P_STANDARD']


def is_fraction(s):
    return s.lower() in ['fraction', 'fractions']


df = pd.read_pickle("GPA 2145-16 Compound Properties Table - English.pkl")

# df = df[df['Compound'] != 'nitrogen']

comp_dict = dict([
    ('methane', 3),
    ('ethane', 3),
    ('propane', 3),
    ('i-butane', 3),
    ('n-pentane', 3),
    ('hexane', 3),
    ('heptane', 3),
    ('octane', 3),
    ('decane', 3),
    ('hydrogen sulfide', 3),
    ('cyclohexane', 3),
    ('nitrogen', 3),
    ('heptadecane', 3),
    # ('fractions', 0.7),
])

comp_dict = normalize_composition(comp_dict)
comps = list(comp_dict.keys())
zs = list(comp_dict.values())

constants = ChemicalConstantsPackage.constants_from_IDs(comps)

# check if the compounds have molecular weight and normal boiling T data
check_properties_exists(constants)

# fixed 0.0236 m^3/mol at standard conditions for all compounds
V_molar = ideal_gas_molar_volume()

ghvs = []
nhvs = []

for cas, name, Hc, mw, rhol_60F_mass in zip(constants.CASs, constants.names, constants.Hcs, constants.MWs,
                                            constants.rhol_60Fs_mass):

    matching_row = df[df['CAS'] == cas]
    print(name, '--------------------------------------')

    # chemical is found in the GPA data table
    if not matching_row.empty:
        ghv_ideal = matching_row['Gross Heating Value Ideal Gas [Btu/ft^3]'].iloc[0]

        if pd.isna(ghv_ideal):
            if Hc != 0:  # chemically inert or contains to combustible energy. Skip them. Ex: nitrogen, argon, helium
                ghv_ideal = Hc / V_molar
                ghv_ideal = ureg('%.15f joule/m^3' % ghv_ideal).to('Btu/ft^3')._magnitude * -1
            else:
                ghv_ideal = 0

    # chemical is NOT identified in the GPA datatable
    else:
        print('No match found in GPA table')

        if Hc == 0:  # chemically inert or contains to combustible energy.
            ghv_ideal = 0
        elif Hc is None:
            """
            sg_liq = correlations.sg_liq(rhol_60F_mass)
            API = newton(lambda API: correlations.API_sg_liq(API, sg_liq), x0=45, maxiter=500)

            # ghv_liquid correlation usig API.

            """
            pass  # Todo: implement a handler that can correlate API to ghv for a working range of the model. Show a warning sign if outside range. Prompt the user to activate ghv_correlate=True
        else:
            ghv_ideal = Hc / V_molar  # Hc = Heat of combustio, J/mol. V_molar = molar volume, m^3/mol
            ghv_ideal = ureg('%.15f joule/m^3' % ghv_ideal).to('Btu/ft^3')._magnitude * -1

    print('ghv_gas (BTU/scf) (Hc/V_molar):    %.1f' % ghv_ideal)

# sg_gas = newton(lambda sg: correlations.mw_sg_gas(mw, sg), x0=0.65, maxiter=50)
# sg_liq = correlations.sg_liq(rhol_60F_mass)
# API = newton(lambda API: correlations.API_sg_liq(API, sg_liq), x0=45, maxiter=500)


comp = dict([
    # ('methane', 1),
    # ('methylamine', 1),
    # ('ethane', 1),
    # ('methylcyclohexane', 1),
    ('docosane', 1),
])

"""
Thermo has way to calculate specific gravity of liquid and gas
"""

