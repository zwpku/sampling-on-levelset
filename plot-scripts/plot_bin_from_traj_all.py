#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
from numpy import *
import matplotlib.pyplot as plt
from math import pi
import matplotlib.gridspec as gridspec

c = 3
fig, axs = plt.subplots(1,3, figsize=(12,3), sharey=True)

scheme_id = 0 
data_file = open('../data/counter_%d.txt' % (scheme_id) , 'r')

n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = 2 * pi / n_bin 
  
counter_data = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]

angle_vec = np.linspace( bin_width * 0.5, 2*pi-bin_width*0.5, n_bin, endpoint=True )

axs[0].plot(angle_vec, counter_data, color='k',linewidth=1.5, label=r'empirical') 
axs[0].axhline(1.0/(2*pi), 0, 2*pi, color='k', linestyle=':', label=r'density of $\mu_1$')
axs[0].set_xlim(0, 2 * pi)
axs[0].set_ylim(0.05, 0.35)
axs[0].set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
axs[0].set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
axs[0].set_yticks([0.1, 0.2, 0.3])
axs[0].tick_params(axis='x', labelsize=20)
axs[0].tick_params(axis='y', labelsize=20)
axs[0].set_title(r'$\Theta$', fontsize=22)
axs[0].legend(frameon=False, fontsize=18, bbox_to_anchor=(1.06, 1.00))

scheme_id = 1 
data_file = open('../data/counter_%d.txt' % (scheme_id) , 'r')

n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = 2 * pi / n_bin 
  
counter_data = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]

angle_vec = np.linspace( bin_width * 0.5, 2*pi-bin_width*0.5, n_bin, endpoint=True )
axs[1].plot(angle_vec, counter_data, color='k',linewidth=1.5, label=r'empirical') 

axs[1].axhline(1.0/(2*pi), 0, 2*pi, color='k', linestyle=':', label=r'density of $\mu_1$')

axs[1].set_xlim(0, 2 * pi)
axs[1].set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
axs[1].set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
axs[1].tick_params(axis='x', labelsize=20)
axs[1].set_title(r'$\Theta^A$', fontsize=22)
axs[1].legend(frameon=False, fontsize=18, bbox_to_anchor=(1.06, 1.00))

scheme_id = 2 
data_file = open('../data/counter_%d.txt' % (scheme_id) , 'r')

n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = 2 * pi / n_bin 
  
counter_data = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]

angle_vec = np.linspace( bin_width * 0.5, 2*pi-bin_width*0.5, n_bin, endpoint=True )
axs[2].plot(angle_vec, counter_data, color='k',linewidth=1.5, label='empirical') 

rho = [np.sqrt(c*c* np.sin(x)*np.sin(x) + np.cos(x)*np.cos(x)) for x in angle_vec]
norm = sum(rho) * bin_width 
axs[2].plot(angle_vec, rho / norm, color='k', linestyle=':',linewidth=1.5, label=r'density of $\mu_2$') 

axs[2].set_xlim(0, 2 * pi)
axs[2].set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
axs[2].set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
axs[2].tick_params(axis='x', labelsize=20)
axs[2].set_title(r'$\Pi$', fontsize=22)
axs[2].legend(frameon=False, fontsize=18, bbox_to_anchor=(1.05, 1.07))

plt.show()

fig_file_name = '../fig/dist.eps' 
savefig(fig_file_name, bbox_inches='tight')
