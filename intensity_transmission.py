from cmath import cos
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import sys
from numpy.linalg import inv
from numpy.linalg import multi_dot
from numpy.linalg import matrix_power

## DAVID VRBA
## Calculation of intensity of electric field


epsilo_a = 1#+1j*0.00000001
mu_a = 1#+1j*0.00000001

epsilo_b = 4+1j*0.1
mu_b = 1+1j*0.00000001  

n_a = np.sqrt(complex(epsilo_a*mu_a))
n_b = np.sqrt(complex(epsilo_b*mu_b))

theta_a = 0 
theta_b = np.arcsin(np.sin(theta_a)*n_a.real/n_b.real)

print(n_a)
print(n_b)

if ((n_a.imag < 0) or (n_b.imag < 0)):
    print("Im(n) je zaporná")
    exit()

l_a = 1.
l_b = 1.
l = l_a + l_b

omega_0 = 2*np.pi/l

N = 10 ## lenght of A+B layers before defect 

koef = 1 ## size of defect a->koef*a
omega = 3*omega_0/8

## Calculation of wave vectors
k_a = omega*n_a*np.cos(theta_a)
k_b = omega*n_b*np.cos(theta_b)
    
chi = (mu_b*k_a)/(mu_a*k_b)

## Trasnfer matrix A->B
M_ab = 1/2*np.array([[1 + chi, 1 - chi], [1 - chi , 1 + chi]])

## Trasnfer matrix B->A
M_ba = inv(M_ab)



## Amplidute of electric field 
E_plus = 1. 
E_minus = 0.

x = []
Intensity = []

E = [E_plus, E_minus]
Intensity.append(np.absolute(E[1]+E[0])**2)

total_lenght = 2*N*l+koef*l_a
number_of_steps = 100000
lenght_of_step = total_lenght/number_of_steps ## dlzka kroku

s_a = int(l_a/lenght_of_step) ## Number of step needed to go through layer A
s_b = int(l_b/lenght_of_step) ## Number of step needed to go through layer B

s = int(s_a + s_b) ## Number of step needed to go through layer A+B

l1 = (s_a+s_b)*N ## Number of step to the defect
l2 = (s_a+s_b)*N+koef*s_a ## Number of step to and through defect 

def transmission(step,K):
    return np.array([[np.exp(1j*K*step), 0], [0, np.exp(-1j*K*step)]])
print("s: ", end ="")
print(s)
print("s_a: ", end ="")
print(s_a)
print("s_b: ", end ="")
print(s_b)


i = 1
while i <= number_of_steps:
    if (i <= l1): ## To the defect
        if 0 < (i % s) < s_a:
            E = transmission(lenght_of_step,k_a).dot(E) ## Dot product
        elif (i % s) == s_a:
            E = M_ab.dot(E)
        elif s_a < (i % s):
            E = np.dot(transmission(lenght_of_step,k_b),E)
        elif (i % s) == 0:
            E = M_ba.dot(E)
        else:
            print("Error!")
    elif ((l1 < i) and (i <= l2)): ## Defect
        E = transmission(lenght_of_step,k_a).dot(E)
    elif (l2 < i): ## After defect
        if ((i-koef*s_a) % s) < s_a:
            E = transmission(lenght_of_step,k_a).dot(E)
        elif ((i-koef*s_a) % s) == s_a:
            E = M_ab.dot(E)
        elif s_a < ((i-koef*s_a) % s):
            E = np.dot(transmission(lenght_of_step,k_b),E)
        elif ((i-koef*s_a) % s) == int(s-1):
            E = M_ba.dot(E)
        else:
            print("Error!")
    else:
        print("Error!")

    x.append(lenght_of_step*(i-1)/total_lenght)
    I = np.absolute(E[1]+E[0])**2
    Intensity.append(I)
    i+=1

Intensity.pop()
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

ax.plot(x, Intezita, linewidth = lw, color = "k", label = "Šírka poruchy - "+str(koef))
ax.legend()
txt = "Počet krokov: "+str(number_of_steps)+", N = "+str(N)+", $l_a$ = "+str(l_a)+", $l_b$ = "+str(l_b)+", $\epsilon_a$ = "+str(epsilo_a)+", $\mu_a$ = "+str(mu_a)+", $\epsilon_b$ = "+str(epsilo_b)+", $\mu_b$ = "+str(mu_b)+r", $\theta_a$ = "+str(theta_a)
plt.axvline(x = 0.5, color = 'r', linestyle = 'dotted')

fig.text(.05,.05,txt)

ax.set_xlabel("Celková dĺžka", labelpad = 10.5, fontsize = sz)
ax.set_ylabel("E$^2$", labelpad = 7.5, fontsize = sz)

ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 

plt.savefig("plot_E.png", dpi = 600, format = 'png')
plt.close()



## vypocet transmisie
print("pocitam T")
n_a = np.sqrt(complex(epsilo_a*mu_a))
n_b = np.sqrt(complex(epsilo_b*mu_b))


if ((n_a.imag < 0) or (n_b.imag < 0)):
    print("Im(n) je zaporná")
    exit()


theta_a = 0 ## kolmy dopad
theta_b = 0
## quarter stack

l_a = 2.
l_b = 1.

l = l_a + l_b

omega_0 = 2*np.pi/l

P = 4 ## cislo pre limit poctu poruch
colors = ["b","g", "r", "c", "m", "y", "k", "w"]

koef = 1 ## koeficient nasobku poruhcy a->koef*a

Transmision = []
w = []

for M in range(1,P):
    Transmision.append([])
    w.append([])
    omega = 0.0001*omega_0

    while omega < 2*omega_0:
        ## vypocet k-cok
        k_a = omega*np.sqrt(n_a**2)*np.cos(theta_a)
        k_b = omega*np.sqrt(n_b**2)*np.cos(theta_b)


        chi = (mu_b*k_a)/(mu_a*k_b)

        ## matica prechodu A->B
        M_ab = 1/2*np.matrix([[1 + chi, 1 - chi], [1 - chi , 1 + chi]])

        ## matica prechodu B->A
        M_ba = inv(M_ab)

        ## matica prechodu cez A 
        M_aa = np.matrix([[np.exp(1j*k_a*l_a), 0], [0, np.exp(-1j*k_a*l_a)]])

        ## matica prechodu cez B
        M_bb = np.matrix([[np.exp(1j*k_b*l_b), 0], [0, np.exp(-1j*k_b*l_b)]])

        ## prechod A->B + B + B->A + A
        M_prechod = multi_dot([M_aa, M_ab, M_bb, M_ba])

        ## M_prechod na pocet vrstiev
        M_powered = matrix_power(M_prechod, N)

        ## matica predstavujuca prechod cez poruchu
        M_ll = np.matrix([[np.exp(1j*k_a*l_a*M*koef), 0], [0, np.exp(-1j*k_a*l_a*M*koef)]])

        ## vypocet celkovej matice prechodu
        M_celkove = multi_dot([M_powered,M_ll, M_powered])


        ## vypocet koeficientu transmisie
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
    ax.plot(w[M-1], Transmision[M-1], linewidth = lw, color = colors[M-1], label = "Šírka poruchy - "+str(M*koef))
ax.legend()
txt = "N = "+str(N)+", $l_a$ = "+str(l_a)+", $l_b$ = "+str(l_b)+", $\epsilon_a$ = "+str(epsilo_a)+", $\mu_a$ = "+str(mu_a)+", $\epsilon_b$ = "+str(epsilo_b)+", $\mu_b$ = "+str(mu_b)+r", $\theta_a$ = "+str(theta_a)

fig.text(.05,.05,txt)

ax.set_xlabel('$\omega/\omega_0$', labelpad = 7.5, fontsize = sz)
ax.set_ylabel("T", labelpad = 7.5, fontsize = sz)

ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 

plt.savefig('plot_T.png', dpi = 600, format = 'png')
plt.close()
