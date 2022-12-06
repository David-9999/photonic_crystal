from cmath import cos


import numpy as np
from numpy.linalg import inv
from numpy.linalg import multi_dot
from numpy.linalg import matrix_power
import random

import matplotlib as mpl
import matplotlib.pyplot as plt

## DAVID VRBA


epsilo_a = 1
mu_a = 1

epsilo_b = 4+1j*0.3
mu_b = 1+1j*0.00000001  

n_a = np.sqrt(complex(epsilo_a*mu_a))
n_b = np.sqrt(complex(epsilo_b*mu_b))


if ((n_a.imag < 0) or (n_b.imag < 0)):
    print("Im(n) je zaporná")
    exit()


theta_a = 0 ## kolmy dopad
theta_b = np.arcsin(np.sin(theta_a)*n_a.real/n_b.real)
## quarter stack

l_a = 10.
l_b = 1.
l = l_a + l_b

omega_0 = 2*np.pi/l

N = 20 ## pocet vrstiev pred poruchou
F = 20 ## pocet vrstiev za poruchou


P = 2 ## Number of defects, can be change.
colors = ["b","g", "r", "c", "m", "y", "k", "w"]

koef = 1 ## koeficient nasobku poruhcy a->koef*a

Transmision = []
w = []
## Calculation for P lenghts of defect 
for M in range(1,P):
    Transmision.append([])
    w.append([])
    omega = 0.0001*omega_0

    while omega <= 1.5*omega_0:
        ## Calculation of wave vectors
        k_a = omega*np.sqrt(n_a**2)*np.cos(theta_a)
        k_b = omega*np.sqrt(n_b**2)*np.cos(theta_b)


        chi = (mu_b*k_a)/(mu_a*k_b)

        ## Trasnfer matrix A->B
        M_ab = 1/2*np.matrix([[1 + chi, 1 - chi], [1 - chi , 1 + chi]])

        ## Trasnfer matrix B->A
        M_ba = inv(M_ab)

        ## Matrix throght layer A
        M_aa = np.matrix([[np.exp(1j*k_a*l_a), 0], [0, np.exp(-1j*k_a*l_a)]])

        ## Matrix throght layer A
        M_bb = np.matrix([[np.exp(1j*k_b*l_b), 0], [0, np.exp(-1j*k_b*l_b)]])

        ## Transition A->B + B + B->A + A
        M_transition = multi_dot([M_aa, M_ab, M_bb, M_ba])

        ## M_transition powered to N
        M_powered = matrix_power(M_transition, N)

        ## Matrix throght defect
        M_ll = np.matrix([[np.exp(1j*k_a*l_a*M*koef), 0], [0, np.exp(-1j*k_a*l_a*M*koef)]])

        ## Total matrix of tranition
        M_celkove = multi_dot([M_powered,M_ll, M_powered])


        ## Transmiton coeficient
        T = 1/np.absolute(M_celkove[1, 1])**2

        Transmision[M-1].append(T)
        w[M-1].append(omega/omega_0)

        omega+=0.0001*omega_0
    



# PLOT
sz = 20
lw = 0.75
wd = 10
hg = 10

fig, ax = plt.subplots(figsize=(wd,hg))

ax.tick_params(which = 'major', axis = 'both', labelsize = sz, direction = 'in', length = 10, pad = 7.5,
	labelbottom = True, labeltop = False, labelleft = True, labelright = False,
	bottom = True, top = True, left = True, right = True)

ax.tick_params(which = 'minor', axis = 'both', labelsize = sz, direction = 'in', length = 5, pad = 7.5,
	labelbottom = True, labeltop = False, labelleft = True, labelright = False,
	bottom = True, top = True, left = True, right = True)

for M in range(1,P):
    ax.plot(w[M-1], Transmision[M-1], linewidth = lw, color = colors[M-1], label = "Size of defect - "+str(M*koef))
ax.legend()
txt = "N = "+str(N)+", $l_a$ = "+str(l_a)+", $l_b$ = "+str(l_b)+", $\epsilon_a$ = "+str(epsilo_a)+", $\mu_a$ = "+str(mu_a)+", $\epsilon_b$ = "+str(epsilo_b)+", $\mu_b$ = "+str(mu_b)+r", $\theta_a$ = "+str(theta_a)

fig.text(.05,.05,txt)

ax.set_xlabel('$\omega/\omega_0$', labelpad = 7.5, fontsize = sz)
ax.set_ylabel("T", labelpad = 7.5, fontsize = sz)

ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 

plt.savefig('plot_T.png', dpi = 600, format = 'png')
plt.close()
