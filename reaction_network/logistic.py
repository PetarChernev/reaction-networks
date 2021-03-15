from reaction_network.abstract import XSReactionNetwork


class LogisticReactionNetwork(XSReactionNetwork):
    REACTIONS = """
        S + X -> 2X
    """
    RATE_NAMES = ['k']

    def __init__(self, rates):
        super(LogisticReactionNetwork, self).__init__(self.REACTIONS, rates)

    def constant(self, masses):
        return sum(masses)

    def get_s_mass(self, x_mass, c=1):
        return c - x_mass