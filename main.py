from reaction_network import ReactionNetwork

network = ReactionNetwork("S + X -> 2X", [2], external_reactants=['Q'])
solution = network.solve_ivp([1, 0.01])
network.plot(solution)