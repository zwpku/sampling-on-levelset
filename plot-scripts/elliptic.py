#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')

from pylab import *
from numpy import *
from math import pi
import numpy as np
from numpy import linalg as LA
from scipy.integrate import odeint
import matplotlib.patches as mpatches

def xi(x):
    return 0.5 * ( x[0]*x[0]/(c*c) + x[1]*x[1] - 1 )

def f(x,t, AA):
    return  [-xi(x) * (x[0] / (c*c) - AA * x[1]), -xi(x)*(AA * x[0] / (c*c) + x[1]) ] 

fig = plt.figure()

lc = ['b', 'r', 'k', 'c', 'm', 'y']

# range of x, [-xb, xb] 
xb = 5
yb = 2
# number of discrete interval
n = 200
c = 3

plt.figure(figsize=(8, 4))
ax = plt.gca()

theta = np.linspace( 0, 2*pi, n, endpoint=True )
plt.plot( c*np.cos(theta), np.sin(theta), color='k' )
#plt.grid( color='lightgray', linestyle='--' )
plt.show() 

ax.spines['top'].set_color('none')
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')

angle = pi / 5 
ll = 1.2
pt = [ c*np.cos(angle), np.sin(angle) ]
normal_dir = [ pt[0] / (c * c), pt[1] ] 
norm = LA.norm(normal_dir) 
normal_dir = normal_dir / norm 

# pt1 is the point x to project
pt1 = [ pt[0] + ll * normal_dir[0],  pt[1] + ll * normal_dir[1] ]

#orthogonal map projection \Pi along a straight line
line=Line2D( [pt[0], pt1[0]], [pt[1], pt1[1]], color='k', linewidth=1.5, linestyle=':')
ax.add_line(line)

scatter(pt1[0], pt1[1], s=13, color='k')
ax.text( pt1[0], pt1[1]+0.08, r'$x$', fontsize=18)

scatter(pt[0], pt[1], s=13, color='k')
ax.annotate( r'$\Pi(x)$', xy=(pt[0]-0.07, pt[1]+0.06), xytext=(pt[0]-0.3,
    pt[1]+1.0), textcoords='data', arrowprops=dict(facecolor='black',
        shrink=0.05, width=0.1, headwidth=5),fontsize=18,
        horizontalalignment='right',verticalalignment='center' )

#compute the flow map \varphi, by solving the ode
nt = 101
t_list = np.linspace(0, 30, nt)
sol = odeint(f, pt1, t_list, args=(0.0,))
plot(sol[:,0], sol[:,1], linewidth=2,color='k')

idx = 3 

plt.arrow(sol[idx,0], sol[idx,1], sol[idx+1,0]-sol[idx,0],
        sol[idx+1,1]-sol[idx,1], shape='full', lw=1.0, facecolor='black', length_includes_head=True, head_width=0.1)

ax.annotate( r'$\varphi(x,s)$', xy=(sol[idx,0]+0.05, sol[idx,1]+0.2),
        xytext=(sol[idx,0]+0.8, sol[idx,1]+0.8), textcoords='data',
        arrowprops=dict(arrowstyle="->", facecolor='black',
        connectionstyle='angle3,angleA=90,angleB=0'),fontsize=18,
        horizontalalignment='left',verticalalignment='center' )

pt2 = [ sol[-1,0], sol[-1,1] ]
scatter(pt2[0], pt2[1], s=13, color='k')


ax.annotate( r'$\Theta(x)$', xy=(pt2[0]-0.02, pt2[1]), xytext=(pt2[0]-0.8, pt2[1]-0.3), 
arrowprops=dict(arrowstyle="->", facecolor='black',
    connectionstyle='angle3,angleA=-20,angleB=-150'),
    fontsize=18,horizontalalignment='right',verticalalignment='center' )

#compute the flow map \varphi^A, by solving the ode (with AA=0.5)
nt = 101
t_list = np.linspace(0, 30, nt)
sol = odeint(f, pt1, t_list, args=(0.5,))
plot(sol[:,0], sol[:,1], linewidth=2,color='k', linestyle='--')

idx = 6 

plt.arrow(sol[idx,0], sol[idx,1], sol[idx+1,0]-sol[idx,0],
        sol[idx+1,1]-sol[idx,1], shape='full', lw=1.0, facecolor='black', length_includes_head=True, head_width=0.1)

ax.annotate( r'$\varphi^A(x,s)$', xy=(sol[idx,0]+0.05, sol[idx,1]+0.2),
        xytext=(sol[idx,0]+0.8, sol[idx,1]+0.2), textcoords='data', 
arrowprops=dict(arrowstyle="->", facecolor='black', connectionstyle='angle3,angleA=-20,angleB=0'),
        fontsize=18, horizontalalignment='left',verticalalignment='center' )


pt2 = [ sol[-1,0], sol[-1,1] ]
scatter(pt2[0], pt2[1], s=15, color='k')
ax.annotate( r'$\Theta^A(x)$', xy=(pt2[0]+0.05, pt2[1]), xytext=(pt2[0]+0.6,
    pt2[1]-1.0), 
arrowprops=dict(arrowstyle="->", facecolor='black',
    connectionstyle='angle3,angleA=120,angleB=-20'),
    fontsize=18,horizontalalignment='left',verticalalignment='center' )
#ax.text( pt2[0]+0.1, pt2[1]-0.08, r'$\Theta(x)$', fontsize=10)

plt.xlim(-xb, xb)
plt.ylim(-1.5, 2)
xticks([-4, -2, 2, 4])
yticks([2])
ax.text(0.05, -0.32, '0', fontsize=18)
ax.text(0.02, -1.15, '-1', fontsize=18)
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)

fs = 30
#ax.set_title('$G_\epsilon(x)$', fontsize=fs)
ax.set_xlabel(r'$x_1$', fontsize=18)
ax.xaxis.set_label_coords(1.0, 0.43)
ax.set_ylabel(r'$x_2$', fontsize=18, rotation=0)
ax.yaxis.set_label_coords(0.53, 0.97)
ax.set_aspect('equal')

plt.tight_layout()
plt.show()

#plt.tight_layout()
fig_file_name = '../fig/elliptic.eps' 
savefig(fig_file_name, bbox_inches='tight')

