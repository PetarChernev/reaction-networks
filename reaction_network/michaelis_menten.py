from reaction_network.abstract import ReactionNetwork


class MichaelisMensenReactionNetwork(ReactionNetwork):
    REACTIONS = """
        E + S -> C
        C -> E + S
        C -> E + P
    """

    RATE_NAMES = ('k_f', 'k_r', 'k_c')

    def __init__(self, rates):
        ReactionNetwork.__init__(self, self.REACTIONS, rates)


if __name__ == '__main__':
    from printer import latex_reaction_network_to_ode_system
    import matplotlib.pyplot as plt
    from scipy.optimize import least_squares
    import numpy as np

    params = [2.38420773, 3.83905267, 2.95569299]
    network = MichaelisMensenReactionNetwork(params)
    solution = network.solve_ivp([0, 1, 0, 1], max_t=10)
    error_sigma = 0.02
    data = np.abs(solution.y + np.random.normal(0, error_sigma, solution.y.shape))
    for i, row in enumerate(data):
        if i == 0:
            continue
        plt.scatter(solution.t, row, label=network.internal_reactants[i])
    plt.legend()
    plt.show()
    data = np.delete(data, 0, axis=0)


    def f(k):
        network = MichaelisMensenReactionNetwork(k)
        result = network.solve_ivp([0, 1, 0, 1], max_t=10)
        result_data = np.delete(result.y, 0, axis=0)
        min_size = min(data.shape[1], result_data.shape[1])
        return (data[:, :min_size] - result_data[:, :min_size]).flatten()


    optimized_population = []
    max_runs = 500
    for i in range(max_runs):
        print(f'Run {i+1}/{max_runs}')
        optimized = least_squares(f, x0=np.random.rand(3))
        optimized_population.append(optimized.x)

    results = np.array(optimized_population)
    fig, axs = plt.subplots(1, 3)
    for i in range(3):
        axs[i].hist(results[:, i], label=f'${network.rate_names[i]}$', bins=30)
    plt.show()
