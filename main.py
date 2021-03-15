from reaction_network.logistic import LogisticReactionNetwork
from reaction_network.gompertz import GompertzReactionNetwork
from reaction_network.mixed import MixedReactionNetwork

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)

max_t = 10
logistic = LogisticReactionNetwork(rates=[1])
solution = logistic.solve_ivp(0.9, max_t)
plt.plot(solution.t, solution.y[1, :], label='Logistic')

mixed = MixedReactionNetwork(rates=[1, 1])
solution = mixed.solve_ivp(0.9, max_t)
plt.plot(solution.t, solution.y[1, :], label='Mixed')

gompertz = GompertzReactionNetwork(rates=[1, 1])
solution = gompertz.solve_ivp(0.9, max_t)
plt.plot(solution.t, solution.y[1, :], label='Gompertz')

plt.legend()
plt.show()