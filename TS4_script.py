
import numpy as np

alpha_max = 0.4
ws = 3

eeps = 10**(alpha_max/10)-1  # eeps = eps al cuadrado
eps = np.sqrt(eeps)
print('eps = ', eps )   
print('eeps = ', eeps)

#%
for n in range(2,6):
    alpha_min_c = 10 * np.log10(1 + eeps * np.cosh(n * np.arccosh(ws))**2)
    print('n = {:d}    ---    alpha_min_cheby {:f}'. format(n, alpha_min_c))

n = 5

#%
from pytc2.general import Chebyshev_polynomials, s, w, print_subtitle
import sympy as sp

chebn_expr = Chebyshev_polynomials(n)

print(sp.expand(chebn_expr))

#%
# obtengo T_ch(s) = T(s)*T(-s) 

# Tch_jw = (1/eeps) / (1/eeps + chebn_expr**2)
Tch_jw = -(1/(eeps*256)) / (-(1/eeps + chebn_expr**2)/256)
j = sp.I

Tch_s = Tch_jw.subs(w, s/j)

print(sp.expand(Tch_s))

#%
from pytc2.sistemas_lineales import analyze_sys, pretty_print_lti, tf2sos_analog, pretty_print_SOS
from pytc2.general import print_latex, print_subtitle
import scipy.signal as sig

num_Tch_s = np.array([(-1)/(eeps*256)])
den_Tch_s = np.array([ 1 , 0 , 5/2 , 0 , 35/16 , 0 , 25/32 , 0 , 25/256 , 0 , -1/(eeps*256) ]) 

Taux = sig.TransferFunction(num_Tch_s , den_Tch_s)
pretty_print_lti(Taux)

#%
# Obtengo las raices

roots_den_Tch_s = np.roots(den_Tch_s)
print(roots_den_Tch_s)

#%
# obtengo las raices de T(s) filtrando las de T(-s)

roots_den_T_s = roots_den_Tch_s[np.real(roots_den_Tch_s) < 0]
print(roots_den_T_s)


# verificacion con cheblap

z,p,k = sig.cheb1ap(n, alpha_max)
num_cheb, den_cheb = sig.zpk2tf(z,p,k)

sos_cheb = tf2sos_analog(num_cheb, den_cheb)

pretty_print_SOS(sos_cheb)
print(' ')
pretty_print_SOS(sos_cheb, mode='omegayq')

#%
# verifico

# %matplotlib qt
from pytc2.sistemas_lineales import analyze_sys, bodePlot, pzmap

all_sys = []

all_sys.append(sig.TransferFunction(num_cheb, den_cheb))

analyze_sys( all_sys )

















