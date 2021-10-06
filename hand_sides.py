import re

from stoichiometry import Stoichiometry


class HandSide:
    def __init__(self, formula: str):
        self.stoichiometry = Stoichiometry()
        for term in formula.strip().split('+'):
            term = term.strip()
            match = re.match(r'(\d*\.?\d*)([A-Z])', term)
            coef = match.group(1)
            if not coef:
                coef = 1
            self.stoichiometry[match.group(2)] = float(coef)
        self.reactants = list(self.stoichiometry.keys())
        self.reactant_set = set(self.reactants)
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