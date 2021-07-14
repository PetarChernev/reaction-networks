import re
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)

from reaction import Reaction


class ReactionNetwork:
    # abstract
    INTERNAL_REACTANTS_NUM = None

    def __init__(self, formula, rates, sep=';\n', external_reactants=None):
        self.external_reactants = set(external_reactants) if external_reactants is not None else set()

        lines = [line.strip() for line in re.split(f"[{sep}]", formula) if line.strip()]
        assert len(lines) == len(rates)
        self.reactions = [Reaction(lines[i], rates[i], self.external_reactants)
                          for i in range(len(lines))]

        self.reactant_set = set().union(*[r.reactant_set for r in self.reactions])
        self.internal_reactants = list(sorted(self.reactant_set - self.external_reactants))
        if self.INTERNAL_REACTANTS_NUM is not None:
            assert len(self.internal_reactants) == self.INTERNAL_REACTANTS_NUM

        self.left_stoichiometry = [r.get_left_stoichiometry_matrix_row(self.internal_reactants) for r in self.reactions]
        self.left_stoichiometry = np.vstack(self.left_stoichiometry).T
        
        self.right_stoichiometry = [r.get_right_stoichiometry_matrix_row(self.internal_reactants) for r in self.reactions]
        self.right_stoichiometry = np.vstack(self.right_stoichiometry).T
        
        self.stoichiometry = [r.get_stoichiometry_matrix_row(self.internal_reactants) for r in self.reactions]
        self.stoichiometry = np.vstack(self.stoichiometry).T
        
        self.rates = rates
        self.rate_laws = [r.build_rate_law(self.internal_reactants) for r in self.reactions]

        self.ode_rhs = self.build_ode_rhs()

    def build_ode_rhs(self):
        def rhs(t, masses):
            flow_rates = [rate_law(masses) for rate_law in self.rate_laws]
            return self.stoichiometry.dot(flow_rates)
        return rhs

    def vectorized_ode_rhs(self, t_arr, mass_arr):
        return np.array([self.ode_rhs(t, masses) for t, masses in zip(t_arr, mass_arr)])

    def solve_ivp(self, initial_masses, max_t=10, max_t_step=0.05):
        return solve_ivp(self.ode_rhs, [0, max_t], initial_masses, max_step=max_t_step)

    def plot(self, solution, show=True, **kwargs):
        for i in range(len(self.internal_reactants)):
            plt.plot(solution.t, solution.y[i, :], label=self.internal_reactants[i], **kwargs)
        plt.legend()
        if show:
            plt.show()

    def initial_cond_plots(self, a, b, bins=100, max_t=10, max_t_step=0.05):
        step = (b - a) / bins
        results = []
        solution = None
        for i in range(1, bins):
            x_0 = a + i * step
            initial_masses = [x_0, 1 - x_0]
            solution = self.solve_ivp(initial_masses, max_t, max_t_step)
            x_t = solution.y[1, :]
            results.append(x_t)
        if solution is None:
            raise ValueError('Empty initial masses loop.')
        results = np.array(results)
        initial = results[:, 1]
        X, Y = np.meshgrid(initial, solution.t)
        ax = plt.gca(projection='3d')
        ax.plot_wireframe(X, Y, results.T, color='black')
        ax.set_xlabel(r'$x_0$')
        ax.set_ylabel(r'$t$')
        ax.set_zlabel(r'$x$')
        plt.tight_layout()
        plt.show()

    def constant(self, masses):
        raise NotImplementedError()

    @staticmethod
    def string_from_stoichiometry(left_stoichiometry, right_stoichiometry, variables):
        reactions = []
        for row in range(len(left_stoichiometry)):
            left_side = []
            for i in range(len(variables)):
                term_str = ''
                if left_stoichiometry[row][i]:
                    if left_stoichiometry[row][i] != 1:
                        term_str += str(left_stoichiometry[row][i])
                    term_str += variables[i].upper()
                if term_str:
                    left_side.append(term_str)
            left_side = ' + '.join(left_side)

            right_side = []
            for i in range(len(variables)):
                term_str = ''
                if right_stoichiometry[row][i]:
                    if right_stoichiometry[row][i] != 1:
                        term_str += str(right_stoichiometry[row][i])
                    term_str += variables[i].upper()
                if term_str:
                    right_side.append(term_str)
            if right_side:
                right_side = ' + '.join(right_side)
            else:
                right_side = 'Q'
            reactions.append(left_side + '->' + right_side)
        reactions = '\n'.join(reactions)
        return reactions

    @staticmethod
    def from_stoichiometry(left_stoichiometry, right_stoichiometry, variables, external_reactants, rates):
        string = ReactionNetwork.string_from_stoichiometry(left_stoichiometry, right_stoichiometry, variables)
        return ReactionNetwork(string, rates, external_reactants=external_reactants)

    def print_ode_system(self, rate_names=None):
        if rate_names is None:
            rate_names = [f"k_{i}" for i in range(len(self.rates))]
        equations = [r.lower() + "'=" for r in self.internal_reactants]
        for i, reaction in enumerate(self.reactions):
            term = rate_names[i]
            for reactant in self.internal_reactants:
                order = reaction.left.stoichiometry[reactant]
                if order == 0:
                    continue
                term += reactant.lower()
                if order > 1:
                    term += '^' + str(order)
            if term == rate_names[i]:
                raise ValueError(f"No reactants in term for reaction {reaction.__repr__()}")
            for j in range(len(equations)):
                multiplicity = reaction.stoichiometry[self.internal_reactants[j]]
                if multiplicity == 0:
                    continue
                equations[j] += '+' if multiplicity > 0 else '-'
                if abs(multiplicity) != 1:
                    equations[j] += str(abs(multiplicity))
                equations[j] += term
        return '\n'.join(equations)

    def __eq__(self, other):
        for r in self.reactions:
            matching_reactions = [other_r for other_r in other.reactions if other_r == r]
            if len(matching_reactions) != 1:
                return False
        return True


class XSReactionNetwork(ReactionNetwork):
    """
    Implements a reaction network with 2 internal reactants - X and S. The X reactant is the one we are most interested
    in, while the S one is some sort of a control over the X one (food source, catalyst, etc.)
    """
    RATE_NAMES = None

    def __init__(self, formula, rates, sep=';\n', external_reactants=None):
        if self.RATE_NAMES is not None:
            assert len(rates) == len(self.RATE_NAMES)
        super().__init__(formula, rates, sep, external_reactants)

    def get_s_mass(self, x_mass, c=1):
        raise NotImplementedError()

    def solve_ivp(self, initial_x_mass, c=1, max_t=10, max_t_step=0.05):
        s_mass = self.get_s_mass(initial_x_mass, c)
        return ReactionNetwork.solve_ivp(self, [s_mass, initial_x_mass], max_t, max_t_step)

    def plot_x(self, solution, show=True, **kwargs):
        plt.plot(solution.t, solution.y[1, :], label=self.label, **kwargs)
        if show:
            plt.show()

    @property
    def label(self):
        name = type(self).__name__.replace('ReactionNetwork', '')
        if self.RATE_NAMES is not None:
            rates = {self.RATE_NAMES[i]: self.rates[i] for i in range(len(self.rates))}
            rates_str = f""" (${', '.join(rate_name + '=' + str(rate_value)
                                          for rate_name, rate_value in rates.items())}$)"""
            rates_str = re.sub(r'0+[1-9]', '', rates_str)
        else:
            rates_str = ''
        return name + rates_str
