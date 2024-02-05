import numpy as np
import config


def Tb_mw(Tb, mw):
    """
    source: [2] (eq 2), or [3] (eq 2.42 + constants from Table 4.5)
    """
    return -mw + ((6.97996 - np.log(1080 - Tb))/0.01964) ** (3/2)

def liq_ghv_sg(ghv, API):
    """
    notes: gross heating value (ghv, also known has high heating value) vs. API gravity for fuel liquids
    source: [4]
    units: ghv (Btu/lb)
    working range: 0 < API < 60
    https://www.cheresources.com/invision/blog/4/entry-297-heats-of-combustion-correlations/#:~:text=HHV%20%2F%20LHV%20%3D%20Higher%20%2F%20Lower%20Heating,molecular%20weight%20of%20the%20fuel%20gas%20%28%3D%20SG%2A28.96%29

    """
    return -ghv + 17721 + 89.08 * API - 0.348 * API**2 + 0.009518 * API**3

def gas_ghv_sg(ghv, sg):
    """
    notes: gross heating value (ghv, also known has high heating value) vs. gas specific gravity for fuel gases
    source: [3]
    units: ghv (Btu/scf)
    working range: < 2.0 sg
    https://www.cheresources.com/invision/blog/4/entry-297-heats-of-combustion-correlations/#:~:text=HHV%20%2F%20LHV%20%3D%20Higher%20%2F%20Lower%20Heating,molecular%20weight%20of%20the%20fuel%20gas%20%28%3D%20SG%2A28.96%29

    """
    return -ghv + 229.60 + 1321 * sg + 207.97 * sg**2 - 57.084 * sg**3


def gas_nhv_sg(nhv, sg):
    """
    notes: net heating value (nhv, also known has low heating value) vs. gas specific gravity for fuel gases
    source: [3]
    units: nhv (Btu/scf)
    """
    return -nhv + 186.37 + 1219.3 * sg + 206.93 * sg**2 - 56.936 * sg**3

def Tb_mw_sg(Tb, mw, sg_liq):
    """
    source: [1] (eq 2.51)
    notes:
    working range: mw 70~700, Tb 300~850K (90-1050F), API 14.4~93.
    errors: 3.5% for mw < 300, 4.7% for mw > 300.
    """
    return -mw + 42.965 * (np.exp(2.097e-4 * Tb - 7.78712 * sg_liq + 2.08476e-3 * Tb * sg_liq)) * Tb**1.26007 * sg_liq**4.983098


def mw_sg_liq(mw, sg_liq):
    """
    source: [2] (eq 3), or [3] (eq 2.42 + constants from Table 4.5). Valid upto C7 ~ C100. Off by 11% for C6.
    """
    return -sg_liq + 1.07 - np.exp(3.56073 - 2.93886 * mw ** 0.1)


def API_sg_liq(API, sg_liq):
    """
    source: [1] (eq 2.4)
    notes: valid for all liquid
    """
    return -API + 141.5 / sg_liq - 131.5


def mw_sg_gas(mw, sg_gas):
    """
    source: [1] (eq 2.6)
    notes: valid for all gas
    """
    sg_air = config.constants['MW_AIR']

    return mw / sg_air - sg_gas


def sg_liq(rhol_60F_mass):
    """
    rhol_60F_mass = Liquid mass densities at 60 °F, [kg/m^3].
    Water liquid density at 60F is assumed to be 999.0170125317171 kg/m^3 by default
    """
    return rhol_60F_mass / config.constants['RHO_WATER']

"""
.. [1] Riazi, M. R.: "Characterization and Properties of Petroleum Fractions," first edition (1985), West Conshohocken, Pennsylvania: ASTM International`
.. [2] Nourozieh, H., Kariznovi,  M., and Abedi, J.: "Measurement and Modeling of Solubility and Saturated - Liquid Density and Viscosity for Methane / Athabasca - Bitumen Mixtures," paper SPE-174558-PA (2016). `(link) <https://onepetro.org/SJ/article/21/01/180/205922/Measurement-and-Modeling-of-Solubility-and>`__
.. [3] API Technical Databook GPA Publication 2145-82
.. [4] Maxwell's Databook on Hydrocarbons
"""