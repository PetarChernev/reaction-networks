from collections import Counter


class HandSide:
    def __init__(self, formula: str):
        self.reactants = [r.strip() for r in formula.strip().split('+')]
        self.reactant_set = set(self.reactants)
        self.stoichiometry = Counter(self.reactants)
        self.formula = formula

    def __repr__(self):
        return self.formula

    def __str__(self):
        return self.formula

    def __eq__(self, other):
        return self.stoichiometry == other.stoichiometry


class LeftHandSide(HandSide):
    pass


class RightHandSide(HandSide):
    pass