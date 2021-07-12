from time import time
from datetime import timedelta
import pickle

from matplotlib.lines import Line2D

from reaction_network.abstract import XSReactionNetwork
from reaction_network.gompertz import GompertzReactionNetwork
from reaction_network.mixed import MixedReactionNetwork
from reaction_network.logistic import LogisticReactionNetwork

import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)


TIME_FORMAT = '%H:%M:%S'
COMPARISON_DELTA = 0.0005


def plot_logistic_gompertz_mixed_asymptote_1(k_logistic,
                                             k_gompertz,
                                             nu_gompertz,
                                             k_mixed,
                                             nu_mixed,
                                             initial_x,
                                             max_t,
                                             x_asymptote):
    logistic = LogisticReactionNetwork(rates=[k_logistic])
    solution = logistic.solve_ivp(initial_x, c=logistic.constant([0, x_asymptote]), max_t=max_t)
    logistic.plot_x(solution, show=False)
    
    gompertz = GompertzReactionNetwork(rates=[k_gompertz, nu_gompertz])
    solution = gompertz.solve_ivp(initial_x, c=gompertz.constant([0, x_asymptote]), max_t=max_t)
    gompertz.plot_x(solution, show=False)
    
    mixed = MixedReactionNetwork(rates=[k_mixed, nu_mixed])
    solution = mixed.solve_ivp(initial_x, c=mixed.constant([0, x_asymptote]), max_t=max_t)
    mixed.plot_x(solution, show=False)
    
    plt.legend()
    plt.xlabel('$t$')
    plt.ylabel('$x(t)$')
    plt.show()


def plot_networks(*networks: XSReactionNetwork, initial_x: float = 0.1, max_t: float = 10):
    for network in networks:
        solution = network.solve_ivp(initial_x, c=network.constant([0, 1]), max_t=max_t)
        network.plot_x(solution, show=False)
    plt.ylabel('$x(t)$')
    plt.legend()
    plt.show()


def plot_state_space(network, x_min, x_max, s_min, s_max, resolution):
    x = np.linspace(x_min, x_max, resolution)
    s = np.linspace(s_min, s_max, resolution)

    states = np.array(np.meshgrid(s, x)).reshape(2, -1)
    velocities = network.vectorized_ode_rhs(np.zeros(resolution**2), states.T).T
    vel_magnitude = np.hypot(*velocities)
    normalized_velocities = velocities / vel_magnitude
    s_states, x_states = states
    s_vel, x_vel = normalized_velocities

    fig, ax = plt.subplots()
    ax.set_xlabel('x')
    ax.set_ylabel('s')

    q = ax.quiver(x_states, s_states, x_vel, s_vel, vel_magnitude, angles='xy')
    cb = plt.colorbar(q)
    cb.set_label('Velocity magnitude')
    plt.title(r'Vector field for the Gompertz model ($k=1, \nu=1$)')
    plt.show()


def plot_xsnetwork_crosses_by_last_rate(xsnetwork_cls_1,
                                        xsnetwork_cls_2,
                                        rate_1_min,
                                        rate_1_max,
                                        rate_2_min,
                                        rate_2_max,
                                        resolution_1,
                                        resolution_2,
                                        x_0=0.1,
                                        x_inf=1,
                                        max_t=10,
                                        t_step=0.05):
    rate_1_array = np.linspace(rate_1_min, rate_1_max, resolution_1)
    rate_2_array = np.linspace(rate_2_min, rate_2_max, resolution_2)
    rate_combinations = np.array(np.meshgrid(rate_1_array, rate_2_array)).reshape(2, -1)

    network_1_other_rates = [1] * (len(xsnetwork_cls_1.RATE_NAMES) - 1)
    network_2_other_rates = [1] * (len(xsnetwork_cls_2.RATE_NAMES) - 1)

    start_time = time()
    counter = 0
    print(f"Integrating {rate_combinations.shape[1]} pairs of reaction networks...")

    x_1s = []
    x_2s = []

    for rate_1, rate_2 in rate_combinations.T:
        network_1 = xsnetwork_cls_1(network_1_other_rates + [rate_1])
        network_1_c = network_1.constant([0, x_inf])
        solution_1 = network_1.solve_ivp(x_0, c=network_1_c, max_t=max_t, max_t_step=t_step)
        x_1s.append(solution_1.y[1, :])
        
        network_2 = xsnetwork_cls_2(network_2_other_rates + [rate_2])
        network_2_c = network_2.constant([0, x_inf])
        solution_2 = network_2.solve_ivp(x_0, c=network_2_c, max_t=max_t, max_t_step=t_step)
        x_2s.append(solution_2.y[1, :])

        counter += 1
        if not counter % 100:
            elapsed_time_seconds = time() - start_time
            elapsed_time = timedelta(seconds=elapsed_time_seconds)

            expected_time_seconds = elapsed_time_seconds * (rate_combinations.shape[1] / counter - 1)
            expected_time = timedelta(seconds=expected_time_seconds)
            print(f"({counter}/{rate_combinations.shape[1]}) |"
                  f" Elapsed time: {elapsed_time} |"
                  f" Expected time to finish: {expected_time}")

    with open('colors_бк', 'wb') as file:
        pickle.dump({'x1': x_1s, 'x2': x_2s}, file)

    colors = []
    for i in range(len(x_1s)):
        x_1 = x_1s[i]
        x_2 = x_2s[i]
        min_length = min([len(x_1), len(x_2)])
        x_1 = x_1[:min_length]
        x_2 = x_2[:min_length]
        difference = x_1 - x_2
        #difference = difference[abs(difference) > COMPARISON_DELTA]
        if np.all(difference >= 0):
            colors.append([0, 150, 0])
        elif np.all(difference <= 0):
            colors.append([0, 0, 150])
        else:
            colors.append([150, 0, 0])

    colors = np.array(colors).reshape(resolution_2, resolution_1, 3)
    plt.imshow(colors, origin='lower', extent=[rate_1_min, rate_1_max, rate_2_min, rate_2_max])
    plt.xlabel('$k_L$')
    plt.ylabel(r'$\nu$')
    legend_elements = [Line2D([0], [0], lw=0, marker='s', color=[0, 150 / 255, 0], markerfacecolor=[0, 150 / 255, 0],
                              label=r'$x_L(t) > x_G(t), \forall t > 0$', markersize=3),
                       Line2D([0], [0], lw=0, marker='s', color=[0, 0, 150 / 255], markerfacecolor=[0, 0, 150 / 255],
                              label=r'$x_L(t) < x_G(t), \forall t > 0$', markersize=3),
                       Line2D([0], [0], lw=0, marker='s', color=[150 / 255, 0, 0], markerfacecolor=[150 / 255, 0, 0],
                              label=r'$\exists t: x_L(t) = x_G(t)$', markersize=3),
                       ]
    plt.legend(handles=legend_elements, bbox_to_anchor=(0.8, -0.3))
    plt.title(r'Logistic ($k_L$) vs. Gompertz ($k_G=1, \nu$) for different $k_L$ and $\nu$ ($x(0)=0.1, x(\infty)=1$)')
    plt.show()


if __name__ == '__main__':
    # network = LogisticReactionNetwork([1])
    # plot_state_space(network, 1, 10, 1, 10, 10)

    xsnetwork_cls_1 = LogisticReactionNetwork
    xsnetwork_cls_2 = GompertzReactionNetwork
    rate_1_min = 0.1
    rate_1_max = 20
    rate_2_min = 0.1
    rate_2_max = 20
    resolution_1 = 256
    resolution_2 = 256
    x_0 = 0.1
    x_inf = 1
    max_t = 10
    t_step = 0.05
    # plot_xsnetwork_crosses_by_last_rate(xsnetwork_cls_1,
    #                                     xsnetwork_cls_2,
    #                                     rate_1_min,
    #                                     rate_1_max,
    #                                     rate_2_min,
    #                                     rate_2_max,
    #                                     resolution_1,
    #                                     resolution_2,
    #                                     x_0,
    #                                     x_inf,
    #                                     max_t,
    #                                     t_step)

    with open('colors_бк', 'rb') as file:
        data = pickle.load(file)

    x_1s = data['x1']
    x_2s = data['x2']
    print('laoded data')
    colors = []
    for i in range(len(x_1s)):
        x_1 = x_1s[i]
        x_2 = x_2s[i]
        min_length = min([len(x_1), len(x_2)])
        x_1 = x_1[:min_length]
        x_2 = x_2[:min_length]
        difference = x_1 - x_2
        difference = difference[abs(difference) > COMPARISON_DELTA]
        if np.all(difference >= 0):
            colors.append([0, 150, 0])
        elif np.all(difference <= 0):
            colors.append([0, 0, 150])
        else:
            colors.append([150, 0, 0])
    print('determined colors')
    colors = np.array(colors).reshape(resolution_2, resolution_1, 3)
    plt.imshow(np.transpose(colors, [1,0,2]), origin='lower', extent=[rate_1_min, rate_1_max, rate_2_min, rate_2_max])
    plt.xlabel('$k_L$')
    plt.ylabel(r'$\nu$')
    legend_elements = [Line2D([0], [0], lw=0, marker='s', color=[0, 150 / 255, 0], markerfacecolor=[0, 150 / 255, 0],
                              label=r'$x_L(t) > x_G(t), \forall t > 0$', markersize=3),
                       Line2D([0], [0], lw=0, marker='s', color=[0, 0, 150 / 255], markerfacecolor=[0, 0, 150 / 255],
                              label=r'$x_L(t) < x_G(t), \forall t > 0$', markersize=3),
                       Line2D([0], [0], lw=0, marker='s', color=[150 / 255, 0, 0], markerfacecolor=[150 / 255, 0, 0],
                              label=r'$\exists t: x_L(t) = x_G(t)$', markersize=3),
                       ]
    plt.legend(handles=legend_elements, loc='lower right')
    plt.title(r'Logistic ($k_L$) vs. Gompertz ($k_G=1, \nu$) for different $k_L$ and $\nu$ ($x(0)=0.1, x(\infty)=1$)')
    plt.show()

