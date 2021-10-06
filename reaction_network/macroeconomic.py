from reaction_network.abstract import ReactionNetwork


class MacroeconomicReactionNetwork(ReactionNetwork):
    REACTIONS = """
        {alpha}K + {beta}C + P -> Y + H + P
        {gamma}Y + {delta}H + P -> C + P
        {sigma}C -> K
    """

    def __init__(self, alpha, beta, gamma, delta, rates):
        assert 0 <= alpha <= 1
        reactions = self.REACTIONS.format(alpha=alpha,
                                          beta=beta,
                                          gamma=gamma,
                                          delta=delta,
                                          sigma=1-beta)
        ReactionNetwork.__init__(self, reactions, rates, external_reactants={'Q'})


if __name__ == '__main__':
    network = MacroeconomicReactionNetwork(1, 0.5, 1, 1, [1, 1, 1])
    from printer import latex_reaction_network_to_ode_system
    print(latex_reaction_network_to_ode_system(network))
    solution = network.solve_ivp([1, 1, 1, 1, 1], max_t=20)
    network.plot(solution)
