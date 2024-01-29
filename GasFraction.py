import numpy as np
from scipy.optimize import newton
import correlations

class GasFraction(object):

    def __init__(self, mw=None, sg=None, VABP=None, ghv=None, nhv=None, Pc=None, Tc=None, omega=None, Tb=None):
        # Note that 'sg' is assumed to be 'sg_gas' and there's no 'api' attribute
        self.attributes = {
            'mw': mw,
            'sg_gas': sg,  # Renamed for clarity within the class that this is sg for gas
            'VABP': VABP,
            'ghv': ghv,
            'nhv': nhv,
            'Pc': Pc,
            'Tc': Tc,
            'omega': omega,
            'Tb': Tb,
            '_sg_liq': None  # Internal attribute for liquid specific gravity, calculated later
        }

        self.correlations = {
            correlations.mw_sg_gas: ['mw', 'sg_gas'],
            correlations.mw_sg_liq: ['mw', '_sg_liq'],
            correlations.Tb_mw_sg: ['Tb', 'mw', '_sg_liq'],
            correlations.gas_ghv_sg: ['ghv', 'sg_gas'],
            correlations.gas_nhv_sg: ['nhv', 'sg_gas'],
        }
        self.resolve_dependencies()

        # Round numerical attributes to 5 decimal places
        n = 5
        self.attributes = {key: round(value, n) if isinstance(value, float) else value for key, value in self.attributes.items()}

    def resolve_dependencies(self):
        resolved = set([attr for attr, value in self.attributes.items() if value is not None])

        # Calculate _sg_liq from mw if mw is provided but _sg_liq is not
        if 'mw' in self.attributes and self.attributes['mw'] is not None and '_sg_liq' not in resolved:
            try:
                self.attributes['_sg_liq'] = newton(lambda sg_liq: correlations.mw_sg_liq(self.attributes['mw'], sg_liq), x0=self.get_initial_guess('_sg_liq'))
                resolved.add('_sg_liq')
            except RuntimeError as e:
                print("Error in calculating _sg_liq from mw: {}".format(e))
                return

        # Resolve other dependencies
        while len(resolved) < len(self.attributes):
            resolved_this_iteration = False
            for correlation_func, variables in self.correlations.items():
                unresolved_vars = [var for var in variables if var not in resolved]
                if len(unresolved_vars) == 1:
                    unresolved_var = unresolved_vars[0]
                    resolved_vars = [self.attributes[var] for var in variables if var in resolved]
                    try:
                        self.attributes[unresolved_var] = newton(lambda x: correlation_func(*self.prepare_args(correlation_func, x, resolved_vars)), x0=self.get_initial_guess(unresolved_var))
                        resolved.add(unresolved_var)
                        resolved_this_iteration = True
                    except RuntimeError as e:
                        print("Error in calculating {}: {}".format(unresolved_var, e))
                        return

            if not resolved_this_iteration:
                break

    def prepare_args(self, correlation_func, x, resolved_vars):
        arg_order = self.correlations[correlation_func]
        args = []
        for arg in arg_order:
            if arg in self.attributes and self.attributes[arg] is not None:
                args.append(self.attributes[arg])
            else:
                args.append(x)
        return args

    def get_initial_guess(self, variable):
        initial_guesses = {'mw': 100, 'api': 30, 'sg_liq': 0.8, 'sg_gas': 0.6, 'Tb': 300, 'ghv': 3000, 'nhv': 3000}
        return initial_guesses.get(variable, 1.0)

    


# Example usage

# Example usage
print('---------------------- gas ------------------------')
phase = 'gas'
a = GasFraction(mw=175, sg=None)
print(a.attributes)

a = GasFraction(mw=None, sg=6.04073)
print(a.attributes)

a = GasFraction(mw=175.1, sg=6.04)
print(a.attributes)

print('------------ Upton Axis C7+ ------------')
a = GasFraction(mw=None, sg=3.464)
print(a.attributes)

a = GasFraction(mw=96.82, sg=None)
print(a.attributes)

print('------------ Upton Axis Whole ------------')
a = GasFraction(mw=59.44, sg=2.126)
print(a.attributes)

print('------------ Brazos Gas ------------')
a = GasFraction(mw=90.161)
print(a.attributes)

a = GasFraction(mw=None, sg=3.1228)
print(a.attributes)

a = GasFraction(mw=90.161, sg=3.1228)
print(a.attributes)

a = GasFraction(sg=1.035)
print(a.attributes)

print('------------ Colorado Facility ------------')
a = GasFraction(sg=1.1)
print(a.attributes)

print('------------ Test ------------')
a = GasFraction(sg=0.5537)
print(a.attributes)



"""
file:///C:/Users/EricKim/Downloads/Upton_-_Axis_-_Storey_Ranch_Breathing_Vapors_Final%20(4).pdf

Compared the analysis file for C6+ fractions for mw vs. sg correlation. Verify whether I gotta use
(eq 4) of this paper: file:///C:/Users/EricKim/Documents/AP42/AP42-Emissions-py/spe-174558-pa.pdf

Also compare with Brazos. There seems to be a big understanding I need to make for converting sg_gas to sg_liq

or 

the above implemented correlation from the Riazi book
"""
