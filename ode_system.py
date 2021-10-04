"""
Implements functionality to build a reaction network from a system of ODEs.
"""
from typing import List
import re

import matplotlib
from multiset import FrozenMultiset

from plot import plot_xs_state_space_trajectories
from reaction_network.abstract import ReactionNetwork


class ODESystem:
    def __init__(self, equations: str, variables: List[str]):
        self.equations_str = equations.split('\n')
        self.equations = [VariableDynamics(*line.split('='), variables) for line in self.equations_str]

    def to_reaction_network(self, reaction_rates, reaction_rate_names=None):
        unique_terms = list(set([t for eq in self.equations for t in eq.terms]))
        left_stoichiometry = []
        right_stoichiometry = []
        total_stoichiometry = []
        for term in unique_terms:
            left_row = []
            right_row = []
            total_row = []
            for eq in self.equations:
                variable_order = term.variables[eq.variable]
                left_row.append(variable_order)
                if term not in eq.terms:
                    total_row.append(0)
                    right_row.append(variable_order)
                else:
                    total_row.append(eq.term_multiplicity[term])
                    right_row.append(variable_order + eq.term_multiplicity[term])
            left_stoichiometry.append(left_row)
            right_stoichiometry.append(right_row)
            total_stoichiometry.append(total_row)
        variables = [eq.variable for eq in self.equations]

        reaction_rates = [reaction_rates[term.coef] for term in unique_terms]
        if reaction_rate_names is not None:
            reaction_rate_names = [reaction_rate_names[term.coef] for term in unique_terms]
        else:
            reaction_rate_names = [f'k_{i}' for i in range(1, len(unique_terms) + 1)]
        return ReactionNetwork.from_stoichiometry(
            left_stoichiometry,
            right_stoichiometry,
            variables,
            external_reactants={'Q'},
            rates=reaction_rates,
            rate_names=reaction_rate_names
        )


class VariableDynamics:
    r"""
    Represents an ODE of the form \dot a_i(t) = \sum_{j=1}^p k_j \prod_{m=1}^n a_m(t).
    For this implementation, the variables a_i and coefficients k_i are single lowercase characters ('a', 'b', etc.)
    Example:
        \dot a(t) = kba - nc

    TODO: Currently does not support higher orders powers of the variables (eg. no kb^2a).
    """
    def __init__(self, variable, right_hand_side, variables_names):
        self.variable = variable.strip("'")
        self.right_hand_side = right_hand_side
        terms = list(re.findall(r"([-+]?\d*)?(\w+)", right_hand_side))
        self.terms = []
        self.term_multiplicity = {}
        for t in terms:
            multiplicity, letters = t
            if multiplicity == '-':
                multiplicity = -1
            elif multiplicity == '+' or not multiplicity:
                multiplicity = 1
            multiplicity = int(multiplicity)

            variables = []
            coef = None
            for letter in letters:
                if letter not in variables_names:
                    if coef is not None:
                        raise ValueError(f"More than one coefficient found in term {''.join(t)}.")
                    coef = letter
                else:
                    variables.append(letter)
            term = Term(coef, variables)
            self.terms.append(term)
            self.term_multiplicity[term] = multiplicity

    def __repr__(self):
        return f"{self.variable}'={self.right_hand_side}"


class Term:
    def __init__(self, coef, variables):
        self.coef = coef
        self.variables = FrozenMultiset(variables)

    def __repr__(self):
        return f"{self.coef}{''.join(self.variables)}"

    def __hash__(self):
        return hash((self.coef, self.variables))

    def __eq__(self, other):
        return hash(self) == hash(other)


if __name__ == '__main__':
    from printer import latex_reaction_network

    matplotlib.rcParams.update({'font.size': 22})
    self = ODESystem("s'=-ksx+ns\nx'=ksx-mx", ['s', 'x'])
    network = self.to_reaction_network({'k': 1, 'n': 0.3, 'm': 0.7},
                                       {'k': 'k', 'n': r'\nu', 'm': r'\mu'})
    #plot_xs_state_space_trajectories(network, [(2, 0.1 * i) for i in range(1, 21)])
    solution = network.solve_ivp([0.7, 0.3])
    network.plot(solution, show=True)
