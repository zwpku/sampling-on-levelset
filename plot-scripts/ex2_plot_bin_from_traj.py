#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
from numpy import *
import matplotlib.pyplot as plt
from math import pi
from scipy.special import erf
import matplotlib.gridspec as gridspec

plt.figure(figsize=(8, 4))
axs = plt.gca()

id_mat_a_flag = 1
stiff_eps = 0.09

data_file = open('../data/phi_counter_%d.txt' % id_mat_a_flag , 'r')

n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = 2 * pi / n_bin 
counter_data = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]
angle_vec = np.linspace( bin_width * 0.5, 2*pi-bin_width*0.5, n_bin, endpoint=True )
axs.plot(angle_vec, counter_data, color='k',linewidth=1.5, label=r'empirical') 
axs.axhline(1.0/(2*pi), 0, 2*pi, color='k', linestyle=':', label=r'density of $\varphi$')
axs.set_xlim(0, 2 * pi)
#axs.set_ylim(0.05, 0.35)
axs.set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
axs.set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
#axs.set_yticks([0.1, 0.2, 0.3])
axs.tick_params(axis='x', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
axs.set_title(r'$\varphi$', fontsize=22)
axs.legend(frameon=False, fontsize=18, bbox_to_anchor=(1.06, 1.00))

fig_file_name = '../fig/ex2_phi_dist_%d.eps' % id_mat_a_flag
savefig(fig_file_name, bbox_inches='tight')

plt.clf()
axs = plt.gca()

data_file = open('../data/theta_counter_%d.txt' % id_mat_a_flag , 'r')

n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = pi / n_bin 
counter_data = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]
angle_vec = np.linspace( -pi * 0.5 + bin_width * 0.5, pi * 0.5 - bin_width*0.5, n_bin, endpoint=True )
axs.plot(angle_vec, counter_data, color='k',linewidth=1.5, label=r'empirical') 

normal_z = sqrt(2 * pi * stiff_eps) * erf( pi/(2 * sqrt(2*stiff_eps)) )
print 'normalization constant Z=%.4e\n' % normal_z
theta_density = [ exp(-x*x*0.5/stiff_eps) / normal_z for x in angle_vec] 
axs.plot(angle_vec, theta_density, color='k',linestyle=':', label=r'density of $\theta$') 
axs.set_xlim(-pi * 0.5, pi * 0.5)
#axs.set_ylim(0.05, 0.35)
axs.set_xticks([-pi/2, -pi/4, 0, pi/4, pi/2]) 
axs.set_xticklabels([r'$-\pi/2$', r'$-\pi/4$', r'$0$', r'$\pi/4$', r'$\pi/2$'])
#axs.set_yticks([0.1, 0.2, 0.3])
axs.tick_params(axis='x', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
axs.set_title(r'$\theta$', fontsize=22)
axs.legend(frameon=False, fontsize=18, bbox_to_anchor=(1.06, 1.00))

fig_file_name = '../fig/ex2_theta_dist_%d.eps' % id_mat_a_flag
savefig(fig_file_name, bbox_inches='tight')

