#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
from numpy import *
from math import gamma
import matplotlib.pyplot as plt
from math import pi
import matplotlib.gridspec as gridspec

def cdf(x0, k):
    xb = 10
    nn = 200
    dx = (xb - x0) / nn
    s = 0 
    for i in range(nn) :
        x = x0 + dx * (i+0.5)
        s = s + pow(x, k-1) * exp(-x*x*0.5) * dx

    return s / (pow(2, k * 0.5 - 1) * gamma(k * 0.5))

plt.figure(figsize=(8, 4))
axs = plt.gca()

N = int (raw_input('N='))
n_dof = N * (N-1) / 2

max_size_s = 1.0
n_size_s = 1000
bin_width = max_size_s / n_size_s 

x_vec = np.linspace( bin_width * 0.5, max_size_s - bin_width*0.5, n_size_s, endpoint=True )

cdf_vec = [ 1.0 - cdf(sqrt(N / x * 0.5), n_dof) for x in x_vec]

axs.plot(x_vec, cdf_vec, color='k', linewidth=1.5, linestyle ='-', label=r'estimate') 

axs.set_xlim(0, max_size_s)
axs.set_ylim(0.01, 1.0)
#axs.set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
#axs.set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
#axs.set_yticks([0.1, 0.2, 0.3])
axs.tick_params(axis='h', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
#axs.set_title(r'$\Theta$', fontsize=22)
axs.legend(frameon=False, fontsize=18, bbox_to_anchor=(1.06, 1.00))

plt.show()

fig_file_name = '../fig/ex4_rate_estimate_%d.eps' % N

savefig(fig_file_name, bbox_inches='tight')


