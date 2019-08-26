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

N = int (raw_input('N='))
mcmc_flag = int (raw_input('mcmc? (0/1)'))

if mcmc_flag == 1:
    data_file = open('../working_dir/data/ex2_mcmc_traj_%d.txt' % N, 'r')
else :
    data_file = open('../working_dir/data/ex2_no_mcmc_traj_%d.txt' % N, 'r')

n_traj, = [ int(x) for x in data_file.readline().split() ]

trace_data = [ float(x) for x in data_file.readline().split() ]

x_vec = np.linspace( 1, n_traj, n_traj, endpoint=True )

axs.plot(x_vec, trace_data, color='k', linewidth=1.5, linestyle ='-') 

axs.set_xlim(0, n_traj)
#axs.set_ylim(0.05, 0.35)
#axs.set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
#axs.set_xticklabels(['$0$', r'$\pi/2$', r'$\pi$', r'$3\pi/2$', r'$2\pi$'])
#axs.set_yticks([0.1, 0.2, 0.3])
axs.tick_params(axis='x', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
#axs.set_title(r'$\Theta$', fontsize=22)
#axs.legend(frameon=False, fontsize=18, bbox_to_anchor=(1.06, 1.00))

plt.show()

if mcmc_flag == 1:
    fig_file_name = '../fig/ex2_mcmc_trace_traj_%d.eps' % N
else :
    fig_file_name = '../fig/ex2_no_mcmc_trace_traj_%d.eps' % N

savefig(fig_file_name, bbox_inches='tight')

