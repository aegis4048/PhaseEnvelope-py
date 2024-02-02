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


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

ureg = pint.UnitRegistry()

settings = {
    "T_STANDARD": 288.70555,  # Temperature in Kelvin, 60F
    "P_STANDARD": 101325.0,    # Pressure in Pascal, 1 atm
    "R": 8.31446261815324,
}
config.update_config(settings)


def normalize_composition(comp_dict):
    """
    :param comp_dict: un-normalized dictionary of composition. {"CH4": 3, "C2H6", 6}
    :return: normalized dictionary of composition. {"CH4": 0.3333, "C2H6", 0.66666666}
    """
    total_comp = sum(comp_dict.values())
    if total_comp > 0:
        keys = list(comp_dict.keys())
        last_key = keys[-1]
        normalized_values = [v / total_comp for v in comp_dict.values()]

        # Normalize all but the last element
        comp_dict = {keys[i]: normalized_values[i] for i in range(len(keys) - 1)}

        # Adjust the last element so that the sum is exactly 1
        comp_dict[last_key] = 1 - sum(comp_dict.values())

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
            raise ValueError("Chemical name '%s' is recognized but missing a required data (rhol_60F_mass, liquid mass density at 60F)." % name)
        if mw is None:
            raise ValueError("Chemical name '%s' is recognized but missing a required data (mw, molecular weight)." % name)
        if Tb is None:
            raise ValueError("Chemical name '%s' is recognized but missing a required data (Tb, normal boiling temperature)." % name)
        if Hc is None:
            raise ValueError("Chemical name '%s' is recognized but missing a required data (Hc, heat of combustion [J/mol])." % name)


def ideal_gas_molar_volume():
    """
    PV=nRT, where number of moles n=1. Rearranging -> V=RT/P
    R = 8.31446261815324 ((m^3-Pa)/(mol-K))
    T = 288.7056 K, 60F, standard temperature
    P = 101325 Pa, 1 atm, standard pressure
    :return: ideal gas molar volume in a standard condition (m^3/mol)
    """
    return config.constants['R'] * config.constants['T_STANDARD'] / config.constants['P_STANDARD']


def is_fraction(s):
    """
    string detector for petroleum fractions. Specialized codes for fractions are triggered if detected.
    """
    return s.lower() in ['fraction', 'fractions']

def get_ghvs_pure_compounds(constants, df_GPA):
    """
    :param constants: thermo's constants object
    :param GPA_data: pandas dataframe of the GPA 2145-16 Table
    :return:
    """
    ghvs_ideal_gas = []
    df = df_GPA
    V_molar = ideal_gas_molar_volume()  # fixed 0.0236 m^3/mol at standard conditions for all compounds
    for cas, name, Hc, rhol_60F_mass in zip(constants.CASs, constants.names, constants.Hcs, constants.rhol_60Fs_mass):

        matching_row = df[df['CAS'] == cas]

        # chemical is found in the GPA data table
        if not matching_row.empty:
            ghv_ideal_gas = matching_row['Gross Heating Value Ideal Gas [Btu/ft^3]'].iloc[0]

            if pd.isna(ghv_ideal_gas):

                # chemically inert or contains to combustible energy. Skip them. Ex: nitrogen, argon, helium
                if Hc != 0:
                    ghv_ideal_gas = Hc / V_molar
                    ghv_ideal_gas = ureg('%.15f joule/m^3' % ghv_ideal_gas).to('Btu/ft^3')._magnitude * -1
                else:
                    ghv_ideal_gas = 0

        # chemical is NOT identified in the GPA datatable
        else:

            if Hc == 0:  # chemically inert or contains to combustible energy.
                ghv_ideal_gas = 0
            elif Hc is None:
                """
                sg_liq = correlations.sg_liq(rhol_60F_mass)
                API = newton(lambda API: correlations.API_sg_liq(API, sg_liq), x0=45, maxiter=500)

                # ghv_liquid correlation usig API.

                """
                pass  # Todo: implement a handler that can correlate API to ghv for a working range of the model. Show a warning sign if outside range. Prompt the user to activate ghv_correlate=True
            else:
                ghv_ideal_gas = Hc / V_molar
                ghv_ideal_gas = ureg('%.15f joule/m^3' % ghv_ideal_gas).to('Btu/ft^3')._magnitude * -1

        ghvs_ideal_gas.append(ghv_ideal_gas)
    return np.array(ghvs_ideal_gas)


df = pd.read_pickle("GPA 2145-16 Compound Properties Table - English.pkl")

#df = df[df['Compound'] != 'nitrogen']

# Brazos Gas
comp_dict = dict([
    ('hydrogen sulfide', 0.001),
    ('nitrogen', 2.304),
    ('carbon dioxide', 1.505),
    ('methane', 71.432),
    ('ethane', 11.732),
    ('propane', 7.595),
    ('isobutane', 0.827),
    ('n-butane', 2.540),
    ('i-pentane', 0.578),
    ('n-pentane', 0.597),
    ('fractions', 0.889),
])

# StateCordell VRU
comp_dict = dict([
    ('hydrogen sulfide', 0.001),
    ('nitrogen', 0.1590),
    ('carbon dioxide', 0.6910),
    ('methane', 51.0870),
    ('ethane', 19.9110),
    ('propane', 14.8830),
    ('isobutane', 2.6940),
    ('n-butane', 6.0390),
    ('isopentane', 1.4840),
    ('n-pentane', 1.5790),
    ('fractions', 1.4720),
])

# Thurmond
comp_dict = dict([
    ('carbon dioxide', 0.261),
    ('nitrogen', 1.295),
    ('methane', 86.878),
    ('ethane', 5.659),
    ('propane', 3.237),
    ('isobutane', 0.405),
    ('n-butane', 1.051),
    ('neopentane', 0),
    ('isopentane', 0.306),
    ('n-pentane', 0.381),
    ('fractions', 0.527),
])

# Fesco - Combs, HC Liq - HT
comp_dict = dict([
    ('hydrogen sulfide', 0.001),
    ('nitrogen', 0.103),
    ('carbon dioxide', 1.485),
    ('methane', 34.352),
    ('ethane', 24.949),
    ('propane', 24.893),
    ('isobutane', 3.633),
    ('n-butane', 7.661),
    ('neopentane', 0.017),
    ('isopentane', 1.240),
    ('n-pentane', 1.108),
    ('hexane', 0.409),
    ('fractions', 0.15),
])

# Kendall - Mirror VRU Discharge
comp_dict = dict([
    ('hydrogen sulfide', 0.01),
    ('nitrogen', 0.0810),
    ('carbon dioxide', 0.3810),
    ('methane', 18.5830),
    ('ethane', 22.6350),
    ('propane', 30.2270),
    ('isobutane', 5.5680),
    ('n-butane', 14.5620),
    ('isopentane', 3.2510),
    ('n-pentane', 2.7770),
    ('fractions', 1.9250),
])

# Fesco - Combs, Sep Gas
comp_dict = dict([
    ('hydrogen sulfide', 0),
    ('nitrogen', 0.135),
    ('carbon dioxide', 1.707),
    ('methane', 47.167),
    ('ethane', 20.530),
    ('propane', 18.116),
    ('isobutane', 2.928),
    ('n-butane', 6.460),
    ('neopentane', 0.020),
    ('isopentane', 1.088),
    ('n-pentane', 0.965),
    ('hexane', 0.540),
    ('fractions', 0.344),
])

# Storey Ranch
comp_dict = dict([
    ('hydrogen sulfide', 0.001),
    ('nitrogen', 0),
    ('carbon dioxide', 0.428),
    ('methane', 1.078),
    ('ethane', 8.520),
    ('propane', 25.730),
    ('isobutane', 6.757),
    ('n-butane', 23.3999),
    ('neopentane', 0.056),
    ('i-pentane', 8.021),
    ('n-pentane', 9.569),
    ('hexane', 9.467),
    ('fractions', 6.975),
])

comp_dict = normalize_composition(comp_dict)

comp_dict_pure = {}
comp_dict_fraction = {}
for key, value in comp_dict.items():
    if is_fraction(key):
        comp_dict_fraction[key] = value
    else:
        comp_dict_pure[key] = value

comps_pure = list(comp_dict_pure.keys())
zs_pure = np.array(list(comp_dict_pure.values()))
zs_fraction = list(comp_dict_fraction.values())[0]

constants_pure = ChemicalConstantsPackage.constants_from_IDs(comps_pure)

# check if the compounds have molecular weight and normal boiling T data
check_properties_exists(constants_pure)

ghvs_pure = get_ghvs_pure_compounds(constants_pure, df_GPA=df)
wghtd_ghvs_pure = ghvs_pure * zs_pure

ghv = 3417

ghv_fraction = round((ghv - sum(wghtd_ghvs_pure)) / zs_fraction, 1)
whgtd_ghv_fraction = ghv_fraction * zs_fraction

###########

fraction_row = ['fraction', None, ghv_fraction, zs_fraction, whgtd_ghv_fraction]
sum_row = ['Total', None, None, (sum(zs_pure) + zs_fraction) * 100, sum(wghtd_ghvs_pure) + whgtd_ghv_fraction]

df_ghvs_table = pd.DataFrame(data=[comps_pure, constants_pure.CASs, ghvs_pure, zs_pure * 100, wghtd_ghvs_pure]).T
df_ghvs_table.columns = ['Compound Name', 'CAS', 'Ideal Gas GHV [Btu/scf]', 'Mole Frac. [%]', 'Wghtd. Ideal Gas GVH [Btu/scf]']
df_ghvs_table.loc[len(df_ghvs_table)] = fraction_row
df_ghvs_table.loc[len(df_ghvs_table)] = sum_row
print(df_ghvs_table.to_string())
