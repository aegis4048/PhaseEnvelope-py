from fluids.constants import R

# Configuration dictionary with default values
constants = {
    "T_STANDARD": 288.7056,  # Temperature in Kelvin
    "P_STANDARD": 101325.0,    # Pressure in Pascal
    "R": 8.31446261815324,
    "MW_AIR": 28.97,  # molecular weight of air at standard conditions, g/mol
    "RHO_WATER": 999.0170125317171,  # density of water @60F, 1atm (kg/m^3) according to IAPWS-95 standard. Calculate rho at different conditions by:  chemicals.iapws95_rho(288.706, 101325) (K, pascal)
}

def update_config(user_config):
    """
    Update configuration values using a user-provided dictionary.
    :param user_config: A dictionary containing configuration keys and their new values
    """
    constants.update(user_config)


