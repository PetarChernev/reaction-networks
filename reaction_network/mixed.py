from reaction_network.abstract import XSReactionNetwork
import numpy as np


class MixedReactionNetwork(XSReactionNetwork):
    REACTIONS = """
        S + X -> 2X + S
        S -> X
    """
    RATE_NAMES = ['k', r'\nu']

    def __init__(self, rates):
        super().__init__(self.REACTIONS, rates)
        self.gamma = self.rates[0] / self.rates[1]

    def constant(self, masses):
        return self.gamma * masses[0] + np.log(masses[1] + 1 / self.gamma)

    def get_s_mass(self, x_mass, c=1):
        return (c - np.log(x_mass + 1 / self.gamma)) / self.gamma