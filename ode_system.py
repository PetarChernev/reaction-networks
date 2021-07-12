from typing import List
import re

from multiset import FrozenMultiset
import numpy as np

from reaction_network.abstract import ReactionNetwork


class ODESystem:
    def __init__(self, equations: str, variables: List[str]):
        self.equations_str = equations.split('\n')
        self.equations = [VariableDynamics(*line.split('='), variables) for line in self.equations_str]

    def to_reaction_network(self):
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

        return ReactionNetwork.from_stoichiometry(
            left_stoichiometry,
            right_stoichiometry,
            variables,
            external_reactants={'Q'},
            rates=['_' for _ in range(len(unique_terms))]
        )



class VariableDynamics:
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
    self = ODESystem("s'=-ksx+ns\nx'=ksx-mx", ['s', 'x'])
    print(self.to_reaction_network())