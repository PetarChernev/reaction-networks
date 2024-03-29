import re
from copy import copy

from hand_sides import LeftHandSide, RightHandSide


class Reaction:
    def __init__(self, formula: str, rate: float, external_reactants=None):
        external_reactants = external_reactants or set()
        self.formula = formula
        self.rate = rate

        left, right = formula.split('->')
        self.left = LeftHandSide(left)
        if external_reactants.intersection(self.left.reactant_set):
            raise ValueError("External reactant in left hand side.")
        self.right = RightHandSide(right)

        self.reactant_set = set(self.left.reactants).union(set(self.right.reactants))

        self.stoichiometry = self.right.stoichiometry - self.left.stoichiometry

    def __repr__(self):
        return self.formula

    def __str__(self):
        return self.formula

    def build_rate_law(self, all_reactants):
        def rate_law(masses):
            result = self.rate
            for r in all_reactants:
                if r in self.left.reactant_set:
                    result *= masses[all_reactants.index(r)] ** float(self.left.stoichiometry[r])
            return result
        return rate_law

    def get_stoichiometry_matrix_row(self, all_reactants):
        return [self.stoichiometry[r] for r in all_reactants]

    def get_left_stoichiometry_matrix_row(self, all_reactants):
        return [self.left.stoichiometry[r] for r in all_reactants]

    def get_right_stoichiometry_matrix_row(self, all_reactants):
        return [self.right.stoichiometry[r] for r in all_reactants]

    def __eq__(self, other):
        return self.left == other.left and self.right == other.right


if __name__ == "__main__":
    reaction = Reaction("A + B -> C + B", 1)
