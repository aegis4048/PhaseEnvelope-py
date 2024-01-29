import numpy as np


def Tb_mw(Tb, mw):
    """
    source: [2] (eq 2), or [3] (eq 2.42 + constants from Table 4.5)
    """
    return -mw + ((6.97996 - np.log(1080 - Tb))/0.01964) ** (3/2)

def gas_ghv_sg(ghv, sg):
    """
    notes: gross heating value (ghv, also known has high heating value) vs. gas specific gravity for fuel gases
    source: [3]
    """
    return -ghv + 229.60 + 1321 * sg + 207.97 * sg**2 - 57.084 * sg**3

def gas_nhv_sg(nhv, sg):
    """
    notes: net heating value (nhv, also known has low heating value) vs. gas specific gravity for fuel gases
    source: [3]
    """
    return -nhv + 186.37 + 1219.3 * sg + 206.93 * sg**2 - 56.936 * sg**3

def Tb_mw_sg(Tb, mw, sg_liq):
    return Tb - (mw + 0.5 * sg_liq)


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
    return mw / 28.97 - sg_gas


"""
.. [1] Riazi, M. R.: "Characterization and Properties of Petroleum Fractions," first edition (1985), West Conshohocken, Pennsylvania: ASTM International`
.. [2] Nourozieh, H., Kariznovi,  M., and Abedi, J.: "Measurement and Modeling of Solubility and Saturated - Liquid Density and Viscosity for Methane / Athabasca - Bitumen Mixtures," paper SPE-174558-PA (2016). `(link) <https://onepetro.org/SJ/article/21/01/180/205922/Measurement-and-Modeling-of-Solubility-and>`__
.. [3] API Technical Databook GPA Publication 2145-82
"""