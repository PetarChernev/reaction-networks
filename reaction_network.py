import re
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from reaction import Reaction


class ReactionNetwork:
    def __init__(self, formula, rates, sep=';\n', external_reactants=None):
        self.external_reactants = set(external_reactants) if external_reactants is not None else set()

        lines = re.split(f"[{sep}]", formula)
        assert len(lines) == len(rates)
        self.reactions = [Reaction(lines[i].strip(), rates[i], self.external_reactants)
                          for i in range(len(lines))]

        self.reactant_set = set().union(*[r.reactant_set for r in self.reactions])
        self.internal_reactants = list(sorted(self.reactant_set - self.external_reactants))

        self.stoichiometry = [r.get_stoichiometry_matrix_row(self.internal_reactants) for r in self.reactions]
        self.stoichiometry = np.vstack(self.stoichiometry).T
        self.rate_laws = [r.build_rate_law(self.internal_reactants) for r in self.reactions]

    def build_ode_rhs(self):
        def rhs(t, masses):
            flow_rates = [rate_law(masses) for rate_law in self.rate_laws]
            return self.stoichiometry.dot(flow_rates)
        return rhs

    def solve_ivp(self, initial_masses, max_t=10, max_t_step=0.02):
        return solve_ivp(self.build_ode_rhs(), [0, max_t], initial_masses, max_step=max_t_step)

    def plot(self, solution):
        for i in range(len(self.internal_reactants)):
            plt.plot(solution.t, solution.y[i, :], label=self.internal_reactants[i])
        plt.legend()
        plt.show()