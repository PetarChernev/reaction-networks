from reaction_network.abstract import XSReactionNetwork
import numpy as np


class GompertzReactionNetwork(XSReactionNetwork):
    REACTIONS = """
        S + X -> 2X + S
        S -> Q
    """

    RATE_NAMES = ['k', r'\nu']

    def __init__(self, rates):
        self.gamma = rates[0] / rates[1]
        super().__init__(self.REACTIONS, rates, external_reactants={'Q'})

    def constant(self, masses):
        return self.gamma * masses[0] + np.log(masses[1])

    def get_s_mass(self, x_mass, c=1):
        return (c - np.log(x_mass)) * self.rates[1] / self.rates[0]
