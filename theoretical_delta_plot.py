import numpy as np

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rc('text', usetex=True)

def f(t, k, nu):
    e = np.exp(-k*t)
    x_pow = np.exp(-nu * t)
    return 0.1 - e * 0.1 ** x_pow - (1 - e) * 0.1 ** (1 + x_pow)

x = np.meshgrid(np.linspace(0, 10, 1000), np.linspace(0.1, 20, 256), np.linspace(0.1, 20, 256), indexing='ij')

result = f(*x)

def color(x):
    x = x[abs(x) > 0.00005]
    if all(x > 0):
        return 1
    if all(x < 0):
        return -1
    return 0


colors = np.apply_along_axis(color, 0, result)
colors[colors == -1] = 2

c = np.zeros([256,256,3])
for i in range(256):
    for j in range(256):
        c[i, j, colors[i, j]] = 150

plt.imshow(c / 255, origin='lower', extent=[0.1, 20, 0.1, 20])
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