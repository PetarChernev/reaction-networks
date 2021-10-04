import numpy as np

from reaction_network.abstract import ReactionNetwork


class GompertzBatemanReactionNetwork(ReactionNetwork):
    REACTIONS = """
           P + X -> 2X + P
           S -> P
           P -> Q
       """

    RATE_NAMES = ['k', r'k_1', r'k_2']

    def __init__(self, rates):
        self.gamma = self.rates[0]/self.rates[1]
        super().__init__(self.REACTIONS, rates, external_reactants={'Q'})

    def constant(self, masses):
        return self.gamma * (masses[0] + masses[1]) + np.log(masses[2])

    def get_x_0_for_asymptote_1(self, s_0, p_0):
        return np.exp(-(s_0 + p_0) * self.gamma)


class ModifiedGompertzBatemanReactionNetwork(ReactionNetwork):
    REACTIONS = """
           P + X -> 2X + P
           S -> P
           P -> X
       """

    RATE_NAMES = ['k', r'k_1', r'k_2']

    def __init__(self, rates):
        self.gamma = rates[0] / rates[1]
        super().__init__(self.REACTIONS, rates)

    def constant(self, masses):
        return masses[0] + masses[1] + np.log(self.gamma * masses[2] + 1) / self.gamma

    def get_x_0_for_asymptote_1(self, p_0, s_0):
        return ((self.gamma + 1) / np.exp(self.gamma * (p_0 + s_0)) - 1) / self.gamma
