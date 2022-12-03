from cmath import cos
import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
import sys

## DAVID VRBA
epsilo_a = 1
mu_a = 1

epsilo_b = 4
mu_b = 1

n_a = np.sqrt(complex(epsilo_a*mu_a))
n_b = np.sqrt(complex(epsilo_b*mu_b))

theta_a = 0 ## kolmy dopad
theta_b = 0


if ((n_a.imag < 0) or (n_b.imag < 0)):
    print("Im(n) je zaporná")
    exit()

l_a = 2.
l_b = n_a.real/n_b.real*l_a

l = l_a + l_b

omega_0 = np.pi/(2*l_a*n_a.real)

N = 20 ## pocet vrstiev

koef = 1 ## koeficient nasobku poruhcy a->koef*a

omega = float(sys.argv[1])*omega_0

    ## vypocet k-cok
k_a = omega*n_a*np.cos(theta_a)
k_b = omega*n_b*np.cos(theta_b)
    
chi = (mu_b*k_a)/(mu_a*k_b)

    ## matica prechodu A->B
M_ab = 1/2*np.array([[1 + chi, 1 - chi], [1 - chi , 1 + chi]])

    ## matica prechodu B->A
M_ba = inv(M_ab)



E_plus = 1.
E_minus = 0.
x = []
Intezita = []

E = [E_plus, E_minus]

total_lenght = 2*N*l+koef*l_a
number_of_steps = 1000000
delta = total_lenght/number_of_steps ## dlzka kroku

s_a = int(l_a/delta) ## pocet krokov potrebnych na prejdenie a
s_b = int(l_b/delta) ## pocet krokov potrebnych na prejdenie b

s = int(s_a + s_b) ## pocet krokov potrebnych pre a+b

l1 = (s_a+s_b)*N ## pocet krokov k poruche
l2 = (s_a+s_b)*N+koef*s_a ## pocet krokov za poruchu

def prechod(posunutie,K):
    return np.array([[np.exp(1j*K*posunutie), 0], [0, np.exp(-1j*K*posunutie)]])
print("s: ", end ="")
print(s)
print("s_a: ", end ="")
print(s_a)
print("s_b: ", end ="")
print(s_b)

    # logika vypoctu je takato:
    # viem kolko krokov urobim v prostredi A a B
    # takze viem kolko krokov urobim cez A+B
    # podla toho aky mam zvyskok po deleni poctu krokov cez A+B
    # viem urcite ci ide o prostredie A, rozhranie alebo B
    # ak zvysok i/s je mensi ako s_a som v A
    # ak sa rovna tak som na rozhrani AB a podobne pre ostatne limity
    # Potom si vypocitam kolko krokov l1 spravim dokym sa dostanem k poruche 
    # aplikujem tieto podmienky 
    # vynasobim maticu E z lava prechodovou maticou s l rovnce i*delta
    # podobne to spravim aj ked som v poruche a za poruchou
i = 1
while i <= number_of_steps:
    if (i <= l1): ## k poruche
        if 0 < (i % s) < s_a:
            E = prechod(delta,k_a).dot(E) ## dot product
        elif (i % s) == s_a:
            E = M_ab.dot(E)
        elif s_a < (i % s):
            E = np.dot(prechod(delta,k_b),E)
        elif (i % s) == 0:
            E = M_ba.dot(E)
        else:
            print("problem")
    elif ((l1 < i) and (i <= l2)): ## porucha
        E = prechod(delta,k_a).dot(E)
    elif (l2 < i): ## od poruchy
        if ((i-koef*s_a) % s) < s_a:
            E = prechod(delta,k_a).dot(E)
        elif ((i-koef*s_a) % s) == s_a:
            E = M_ab.dot(E)
        elif s_a < ((i-koef*s_a) % s):
            E = np.dot(prechod(delta,k_b),E)
        elif ((i-koef*s_a) % s) == int(s-1):
            E = M_ba.dot(E)
        else:
            print("problem")
    else:
        print("problem")

    x.append(delta*(i-1)/total_lenght)
    I = np.absolute(E[1]+E[0])**2
    Intezita.append(I)
    i+=1


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

ax.plot(x, Intezita, linewidth = lw, color = "k", label = "$\omega$ = "+str(sys.argv[1])+"$\omega_0$")
ax.legend()
txt = "Počet krokov: "+str(number_of_steps)+", N = "+str(N)+", $l_a$ = "+str(l_a)+", $l_b$ = "+str(l_b)+", $\epsilon_a$ = "+str(epsilo_a)+", $\mu_a$ = "+str(mu_a)+", $\epsilon_b$ = "+str(epsilo_b)+", $\mu_b$ = "+str(mu_b)+r", $\theta_a$ = "+str(theta_a)
plt.axvline(x = 0.5, color = 'r', linestyle = 'dotted')

fig.text(.05,.05,txt)

ax.set_xlabel("Celková dĺžka", labelpad = 10.5, fontsize = sz)
ax.set_ylabel("E$^2$", labelpad = 7.5, fontsize = sz)

ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 

plt.savefig("plot_E_"+str(sys.argv[1])+".png", dpi = 600, format = 'png')
plt.close()

