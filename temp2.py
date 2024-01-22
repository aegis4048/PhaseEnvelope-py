import numpy as np
from scipy.optimize import newton


class PseudoComponent(object):

    def __init__(self, mw=None, sg=None, VABP=None, API=None, Pc=None, Tc=None, omega=None, Tb=None, phase='liquid'):
        if phase not in ['liquid', 'gas']:
            raise TypeError("Unsupported phase type '{}'. Available phase types are ['liquid', 'gas']".format(phase))
        if phase == 'gas' and API is not None:
            raise ValueError("API value is not applicable for the gas phase.")

        self.attributes = {
            'mw': mw,
            'sg': sg,
            'VABP': VABP,
            'API': API,
            'Pc': Pc,
            'Tc': Tc,
            'omega': omega,
            'Tb': Tb,
            'phase': phase
        }

        self.correlations = {
            self.obj_func_correlation_mw_sg_liq: ['mw', 'sg'],
            self.obj_func_correlation_API_sg_liq: ['API', 'sg'],
            self.obj_func_correlation_mw_sg_gas: ['mw', 'sg'],
            self.obj_func_correlation_Tb_mw_sg: ['Tb', 'mw', 'sg']
        }
        self.resolve_dependencies()

        n = 5
        self.attributes = {key: round(value, n) if isinstance(value, float) else value for key, value in self.attributes.items()}

    def resolve_dependencies(self):
        resolved = set([attr for attr, value in self.attributes.items() if value is not None])

        if self.attributes['phase'] == 'liquid':
            if 'API' in self.attributes and self.attributes['API'] is not None and 'sg' not in resolved:
                try:
                    self.attributes['sg'] = newton(lambda sg: self.obj_func_correlation_API_sg_liq(self.attributes['API'], sg), x0=self.get_initial_guess('sg'))
                    resolved.add('sg')
                except RuntimeError as e:
                    print("Error in calculating sg from API: {}".format(e))
                    return

        while len(resolved) < len(self.attributes):
            resolved_this_iteration = False
            for correlation_func, variables in self.correlations.items():
                unresolved_vars = [var for var in variables if var not in resolved]
                if len(unresolved_vars) == 1:
                    unresolved_var = unresolved_vars[0]
                    if self.attributes['phase'] == 'gas' and unresolved_var == 'API':
                        continue  # Skip API calculation for gas phase
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
        initial_guesses = {'mw': 100, 'API': 30, 'sg': 0.8, 'Tb': 300}
        return initial_guesses.get(variable, 1.0)

    def obj_func_correlation_Tb_mw_sg(self, Tb, mw, sg):
        return Tb - (mw + 0.5 * sg)

    def obj_func_correlation_mw_sg_liq(self, mw, sg):
        return -sg + 1.07 - np.exp(3.56073 - 2.93886 * mw ** 0.1)

    def obj_func_correlation_API_sg_liq(self, API, sg):
        """
        source: [1] (eq 2.4)
        """
        return -API + 141.5 / sg - 131.5

    def obj_func_correlation_mw_sg_gas(self, mw, sg):
        """
        source: [1] (eq 2.6)
        """
        return mw / 28.97 - sg


# Example usage
a = PseudoComponent(mw=175, sg=None, phase='liquid')
print(a.attributes)

a = PseudoComponent(mw=None, sg=0.8146543700286001, phase='liquid')
print(a.attributes)

a = PseudoComponent(mw=None, sg=None, API=42.193292770322145, phase='liquid')
print(a.attributes)

a = PseudoComponent(mw=175.1, sg=0.815, phase='liquid')
print(a.attributes)


# Example usage
print('---------------------- gas ------------------------')
phase = 'gas'
a = PseudoComponent(mw=175, sg=None, phase=phase)
print(a.attributes)

a = PseudoComponent(mw=None, sg=0.8146543700286001, phase=phase)
print(a.attributes)

#a = PseudoComponent(mw=None, sg=None, API=42.193292770322145, phase=phase)
#print(a.attributes)

a = PseudoComponent(mw=175.1, sg=0.815, phase=phase)
print(a.attributes)

