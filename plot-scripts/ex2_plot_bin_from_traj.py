#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
from numpy import *
import matplotlib.pyplot as plt
from math import pi
from scipy.special import erf
import matplotlib.gridspec as gridspec

#compute the normalization constant
def int_norm():
  nn = 5000
  dx = pi / nn 
  s = 0 
  for i in range(nn):
      x = -pi * 0.5 + dx * (i + 0.5)
      s += exp(-x * x * 0.5 / stiff_eps) * cos(x)
  return s * dx

fig, [ax1,ax2] = plt.subplots(1,2, figsize=(9,2.5) )

stiff_eps = 0.005

# plot the pdf of angles for a=id
id_mat_a_flag = 1 

# read data corresponding to small step-size
data_file = open('../working_dir/data/ex2_theta_counter_%d_0.txt' % id_mat_a_flag , 'r')
n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = pi / n_bin 
counter_data_0 = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]

# read data corresponding to large step-size
data_file = open('../working_dir/data/ex2_theta_counter_%d_1.txt' % id_mat_a_flag , 'r')
n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = pi / n_bin 
counter_data_1 = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]

angle_vec = np.linspace( -pi * 0.5 + bin_width * 0.5, pi * 0.5 - bin_width*0.5, n_bin, endpoint=True )
ax1.plot(angle_vec, counter_data_0, color='r',linewidth=1.5, linestyle='-.', label=r'$h=0.0002$') 
ax1.plot(angle_vec, counter_data_1, color='k', linestyle='-', linewidth=1.5, label=r'$h=0.005$') 

# plot the exact density
normal_z = sqrt(2 * pi * stiff_eps) * erf( pi/(2 * sqrt(2*stiff_eps)) )
normal_z = int_norm()
print 'normalization constant Z=%.4e\n' % normal_z
theta_density = [ exp(-x*x*0.5/stiff_eps) / normal_z * cos(x) for x in angle_vec] 
ax1.plot(angle_vec, theta_density, color='k',linewidth=1.0, linestyle=':', label=r'exact') 
ax1.set_xlim(-pi * 0.5, pi * 0.5)
ax1.set_ylim(-0.5, 6.0)
ax1.set_xticks([-pi/2, -pi/4, 0, pi/4, pi/2]) 
ax1.set_xticklabels([r'$-\frac{\pi}{2}$', r'$-\frac{\pi}{4}$', r'$0$',
    r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$'])
ax1.set_yticks([0, 2, 4])
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=18)
ax1.set_title(r'$\theta$', fontsize=22)
ax1.legend(frameon=False, fontsize=12, bbox_to_anchor=(0.51, 0.40))

# data for angle varphi, small step-size
data_file = open('../working_dir/data/ex2_phi_counter_%d_0.txt' % id_mat_a_flag , 'r')
n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = 2 * pi / n_bin 
counter_data_0 = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]

# data for angle varphi, large step-size
data_file = open('../working_dir/data/ex2_phi_counter_%d_1.txt' % id_mat_a_flag , 'r')
n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = 2 * pi / n_bin 
counter_data_1 = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]

# exact solution
angle_vec = np.linspace( bin_width * 0.5, 2*pi-bin_width*0.5, n_bin, endpoint=True )

ax2.plot(angle_vec, counter_data_0, color='r',linewidth=1.5, linestyle='-.', label=r'$h=0.0002$') 
ax2.plot(angle_vec, counter_data_1, color='k',linewidth=1.0, linestyle='-', label=r'$h=0.005$') 
ax2.axhline(1.0/(2*pi), 0, 2*pi, color='k', linewidth=1.0, linestyle=':', label=r'exact')
ax2.set_xlim(0, 2 * pi)
ax2.set_ylim(-0.05, 0.22)
ax2.set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
ax2.set_xticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
ax2.set_yticks([0, 0.1, 0.2])
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=18)
ax2.set_title(r'$\varphi$', fontsize=22)
ax2.legend(frameon=False, fontsize=12, bbox_to_anchor=(0.90, 0.55))

fig_file_name = '../fig/ex2_theta_phi_dist_%d.eps' % id_mat_a_flag
savefig(fig_file_name, bbox_inches='tight')

plt.clf()
fig, [ax1,ax2] = plt.subplots(1,2, figsize=(9,2.5) )

# plot when matrix a is non-identity
id_mat_a_flag = 0 

# data for angle theta
data_file = open('../working_dir/data/ex2_theta_counter_%d_1.txt' % id_mat_a_flag , 'r')
n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = pi / n_bin 
counter_data_0 = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]

angle_vec = np.linspace( -pi * 0.5 + bin_width * 0.5, pi * 0.5 - bin_width*0.5, n_bin, endpoint=True )
ax1.plot(angle_vec, counter_data_0, color='r',linewidth=1.5, linestyle='-.', label=r'$h=0.01$') 

normal_z = sqrt(2 * pi * stiff_eps) * erf( pi/(2 * sqrt(2*stiff_eps)) )
normal_z = int_norm()
print 'normalization constant Z=%.4e\n' % normal_z
theta_density = [ exp(-x*x*0.5/stiff_eps) / normal_z * cos(x) for x in angle_vec] 
ax1.plot(angle_vec, theta_density, color='k',linewidth=1.0, linestyle=':', label=r'exact') 
ax1.set_xlim(-pi * 0.5, pi * 0.5)
ax1.set_ylim(-0.5, 6.0)
ax1.set_xticks([-pi/2, -pi/4, 0, pi/4, pi/2]) 
ax1.set_xticklabels([r'$-\frac{\pi}{2}$', r'$-\frac{\pi}{4}$', r'$0$',
    r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$'])
ax1.set_yticks([0, 2, 4])
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(axis='y', labelsize=18)
ax1.set_title(r'$\theta$', fontsize=22)
ax1.legend(frameon=False, fontsize=12, bbox_to_anchor=(0.53, 0.75))

# data for angle varphi
data_file = open('../working_dir/data/ex2_phi_counter_%d_1.txt' % id_mat_a_flag , 'r')
n, n_bin = [ int (x) for x in data_file.readline().split() ]
bin_width = 2 * pi / n_bin 
counter_data_0 = [ (float(x) * 1.0 / (n * bin_width)) for x in data_file.readline().split() ]

angle_vec = np.linspace( bin_width * 0.5, 2*pi-bin_width*0.5, n_bin, endpoint=True )

ax2.plot(angle_vec, counter_data_0, color='r',linewidth=1.5, linestyle='-.', label=r'$h=0.01$') 
ax2.axhline(1.0/(2*pi), 0, 2*pi, color='k', linewidth=1.0, linestyle=':', label=r'exact')
ax2.set_xlim(0, 2 * pi)
ax2.set_ylim(0.05, 0.22)
ax2.set_xticks([0, pi/2, pi, pi*1.5, 2*pi]) 
ax2.set_xticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
ax2.set_yticks([0, 0.1, 0.2])
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=18)
ax2.set_title(r'$\varphi$', fontsize=22)
ax2.legend(frameon=False, fontsize=12, bbox_to_anchor=(0.90, 0.55))

fig_file_name = '../fig/ex2_theta_phi_dist_%d.eps' % id_mat_a_flag
savefig(fig_file_name, bbox_inches='tight')

