from reaction_network import ReactionNetwork

network = ReactionNetwork("S + X -> Y + S; Y -> S", rates=[2, 1])
solution = network.solve_ivp(initial_masses=[0, 10, 0.01])
network.plot(solution)