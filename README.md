# Photonic_crystal

This program can compute intensity of electirc field and transmition coeficient of 1D photonic crystal with defect.

## Theory


Theory is based on [this book](https://www.jstor.org/stable/j.ctt7s2wj) and any good Electromagnetism theory book.
In simplicity photonic crystal is material that disrupts translation invariance of free space and the index of refraction is periodically repeated. The simpliest photonic crystal is repeatiog of two layers of materials. We call this 1D photonic crystal.

Defect is everything that disrupts periodically repeated of index of refraction. For example photonic crystal with defect can be ilustrated by $N$ layers of $2$ materials ($A$ and $B$), then $D$ layer and another $N$ layers of $A$ and $B$.

We gonna need some physical expressions from electromagnetics and optics.
The transition of electric field can be described with transfer matrix.

$$ {\left\lbrack \matrix{E_2^+ \cr E_2^-} \right\rbrack}=M_{\text{T}}{\left\lbrack \matrix{E_1^+ \cr E_1^-} \right\rbrack} $$

The matrix of transition can be created by combination of two types of matrixes. 
First is matrix that describe transmition through layer,

$$ M_{\text{layer}}={\left\lbrack \matrix{e^{ik_z\Delta} & 0 \cr 0 & e^{-ik_z\Delta}} \right\rbrack} $$

and maxtrix of transmition throght interfrace of layers,

$$ M_{\text{TE}}={\left\lbrack \matrix{1+χ_E & 1-χ_E \cr 1-χ_E & 1+χ_E} \right\rbrack} $$

where
$$\chi_E = \frac{\mu_2 k_{1z}}{\mu_1 k_{2z}}$$
$$k_{1,2z}=\frac{\omega}{c}n_{1,2}\cos\left(\theta_{1,2}\right)$$

$\omega$ is frequency of  incident electric field, $\theta$ is angle of incidence, $\mu$ is permeability and $n=\sqrt{\epsilon \mu}$ is refraction index.
Transmition coeficient can be calculated as $T = 1/\left\|M_{22}\right\|^2$ and depend of frequency $\omega$.
We define reference frequency as $\omega_0 = 2\pi/l$, where $l$ is lenght of $A+B$ layer.

## Logic of program

You can see there are $4$ <em>.py</em> programs.


### photonic_transmition.py

The [<em>photonic_transmition.py</em>](../main/photonic_transmition.py) will generate picture of transmition coeficient. This program doesnt take any input so you must change variables by yourself. The logic is based on computing trasnfer matrix and layer matrix and then multipy them in correct order.

### photonic_intesity.py

The [<em>photonic_intesity.py</em>](../main/photonic_intesity.py) will generate picture of intesity of electric field inside photonic crystal. This program doesnt take any input so you must change variables by yourself. The only tricky part of this program is <em>while loop<\em> on line $88$. This loop function as followed:
* Calculate number of step needed to go through layer $A$ $\left(s_a\right)$ and layer $B$ $\left(s_b\right)$ and throught layer $A+B$ $\left(s\right)$
* Calculate number of step to the defect $\left(l_1\right)$ and throught defect $\left(l_2\right)$
* It is good to relize that if $s$ modulo of $i$-step is less that $s_a$ we are calculation electric field in layer $A$. If modulo of $i$ si equals to $s_a$ then we are locaded at the interface of $A$ and $B$, similary other cases.


### usr_E_T.py

The [<em>usr_E_T.py</em>](../main/usr_E_T.py) generetes the intensity of electric field inside of photonic crystal and transmition coeficient for range of $\omega$. You need to input permitivity, permeability, lenghts of layers, angle of incidence, size of defect and number of repeating of $A+B$ layer and multiplication of $\omega_0$ for frequency $\omega$. And it will generate two <em>.png</em> pictures. First <em>plot_E.png</em> shows the intensity of $E$ in photonic crystal, second <em>plot_T.png</em> represent transmission coeficient.

### intensity_transmiton.py

The [<em>intensity_transmiton.py</em>](../main/intensity_transmiton.py) is "identical" copy of [<em>usr_E_T.py</em>](../main/usr_E_T.py) but with no inputs. So you can vary any variable in program. 

