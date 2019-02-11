#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
from numpy import *
import matplotlib.pyplot as plt
from math import pi
import matplotlib.gridspec as gridspec


plt.figure(figsize=(8, 4))
axs = plt.gca()

scheme_id = 0

data_file = open('../data/counter_%d.txt' % scheme_id, 'r')

n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = 2 * pi / n_bin 
  
counter_data = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]

angle_vec = np.linspace( bin_width * 0.5, 2*pi-bin_width*0.5, n_bin, endpoint=True )

axs.plot(angle_vec, counter_data, color='k',linewidth=1.5, label=r'empirical') 
axs.axhline(1.0/(2*pi), 0, 2*pi, color='k', linestyle=':', label=r'density of $\mu_1$')
axs.set_xlim(0, 2 * pi)
#axs.set_ylim(0.05, 0.35)
axs.set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
axs.set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
axs.set_yticks([0.1, 0.2, 0.3])
axs.tick_params(axis='x', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
axs.set_title(r'$\Theta$', fontsize=22)
axs.legend(frameon=False, fontsize=18, bbox_to_anchor=(1.06, 1.00))

plt.show()

fig_file_name = '../fig/dist_%d.eps' % scheme_id
savefig(fig_file_name, bbox_inches='tight')
