from reaction_network.abstract import ReactionNetwork


class SSystemReactionNetwork(ReactionNetwork):
    REACTIONS = "2.34S + X -> 2X"
    RATE_NAMES = ['k']

    def __init__(self, rate):
        ReactionNetwork.__init__(self, self.REACTIONS, [rate])


if __name__ == '__main__':
    from printer import latex_reaction_network_to_ode_system
    network = SSystemReactionNetwork(1)
    print(latex_reaction_network_to_ode_system(network))
    solution = network.solve_ivp([2, 1])
    network.plot(solution)