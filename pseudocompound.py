import numpy as np
from scipy.optimize import newton

class PseudoComponent(object):

    def __init__(self, mw=None, sg_gas=None, sg_liq=None, VABP=None, api=None, ghv=None, lhv=None, Pc=None, Tc=None, omega=None, Tb=None, phase='liquid'):
        if phase not in ['liquid', 'gas']:
            raise TypeError("Unsupported phase type '{}'. Available phase types are ['liquid', 'gas']".format(phase))
        if phase == 'gas' and api is not None:
            raise ValueError("api value is not applicable for the gas phase. Do not input api, or set api=None.")

        self.attributes = {
            'mw': mw,
            'sg_gas': sg_gas,
            'sg_liq': sg_liq,
            'VABP': VABP,
            'api': api,
            'ghv': ghv,
            'lhv': lhv,
            'Pc': Pc,
            'Tc': Tc,
            'omega': omega,
            'Tb': Tb,
            'phase': phase
        }

        self.correlations = {
            self.obj_func_correlation_mw_sg_liq: ['mw', 'sg_liq'],
            self.obj_func_correlation_api_sg_liq: ['api', 'sg_liq'],
            self.obj_func_correlation_mw_sg_gas: ['mw', 'sg_gas'],
            self.obj_func_correlation_Tb_mw_sg: ['Tb', 'mw', 'sg_liq']
        }
        self.resolve_dependencies()

        n = 5
        self.attributes = {key: round(value, n) if isinstance(value, float) else value for key, value in self.attributes.items()}

    def resolve_dependencies(self):
        resolved = set([attr for attr, value in self.attributes.items() if value is not None])

        # Liquid phase: Calculate sg_liq from api if api is provided and sg_liq is not
        if self.attributes['phase'] == 'liquid':
            if 'api' in self.attributes and self.attributes['api'] is not None and 'sg_liq' not in resolved:
                try:
                    self.attributes['sg_liq'] = newton(lambda sg_liq: self.obj_func_correlation_api_sg_liq(self.attributes['api'], sg_liq), x0=self.get_initial_guess('sg_liq'))
                    resolved.add('sg_liq')
                except RuntimeError as e:
                    print("Error in calculating sg_liq from api: {}".format(e))
                    return

        # Resolve other dependencies
        while len(resolved) < len(self.attributes):
            resolved_this_iteration = False
            for correlation_func, variables in self.correlations.items():
                unresolved_vars = [var for var in variables if var not in resolved]
                if len(unresolved_vars) == 1:
                    unresolved_var = unresolved_vars[0]

                    # Skip api calculation for gas phase
                    if self.attributes['phase'] == 'gas' and unresolved_var == 'api':
                        continue

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
        initial_guesses = {'mw': 100, 'api': 30, 'sg_liq': 0.8, 'sg_gas': 0.6, 'Tb': 300}
        return initial_guesses.get(variable, 1.0)

    def obj_func_correlation_Tb_mw_sg(self, Tb, mw, sg_liq):
        return Tb - (mw + 0.5 * sg_liq)

    def obj_func_correlation_mw_sg_liq(self, mw, sg_liq):
        """
        source: [2] (eq 3), or [3] (eq 2.42 + constants from Table 4.5).
        """
        return -sg_liq + 1.07 - np.exp(3.56073 - 2.93886 * mw ** 0.1)

    def obj_func_correlation_api_sg_liq(self, api, sg_liq):
        """
        source: [1] (eq 2.4)
        notes: valid for all liquid
        """
        return -api + 141.5 / sg_liq - 131.5

    def obj_func_correlation_mw_sg_gas(self, mw, sg_gas):
        """
        source: [1] (eq 2.6)
        notes: valid for all gas
        """
        return mw / 28.97 - sg_gas

# Example usage
a = PseudoComponent(mw=175, sg_gas=None, phase='liquid')
print(a.attributes)

a = PseudoComponent(mw=None, sg_gas=None, sg_liq=0.8146543700286001, phase='liquid')
print(a.attributes)

a = PseudoComponent(mw=None, sg_gas=None, sg_liq=None, api=42.193292770322145, phase='liquid')
print(a.attributes)

a = PseudoComponent(mw=175.1, sg_gas=None, sg_liq=0.815, phase='liquid')
print(a.attributes)


# Example usage
print('---------------------- gas ------------------------')
phase = 'gas'
a = PseudoComponent(mw=175, sg_liq=None, sg_gas=None, phase=phase)
print(a.attributes)

a = PseudoComponent(mw=None, sg_liq=None, sg_gas=6.04073, phase=phase)
print(a.attributes)

a = PseudoComponent(mw=175.1, sg_liq=None, sg_gas=6.04, phase=phase)
print(a.attributes)

print('------------ Upton Axis ------------')
a = PseudoComponent(mw=None, sg_liq=None, sg_gas=3.464, phase=phase)
print(a.attributes)

a = PseudoComponent(mw=96.82, sg_liq=None, sg_gas=None, phase=phase)
print(a.attributes)

print('------------ Brazos Gas ------------')
a = PseudoComponent(mw=None, sg_liq=0.7142, sg_gas=None, phase=phase)
print(a.attributes)

a = PseudoComponent(mw=90.161, sg_liq=None, sg_gas=None, phase=phase)
print(a.attributes)

a = PseudoComponent(mw=None, sg_liq=None, api=66.6, phase='liquid')
print(a.attributes)

a = PseudoComponent(mw=None, sg_gas=3.1228, api=None, phase='liquid')
print(a.attributes)



"""
file:///C:/Users/EricKim/Downloads/Upton_-_Axis_-_Storey_Ranch_Breathing_Vapors_Final%20(4).pdf

Compared the analysis file for C6+ fractions for mw vs. sg correlation. Verify whether I gotta use
(eq 4) of this paper: file:///C:/Users/EricKim/Documents/AP42/AP42-Emissions-py/spe-174558-pa.pdf

Also compare with Brazos. There seems to be a big understanding I need to make for converting sg_gas to sg_liq

or 

the above implemented correlation from the Riazi book
"""
