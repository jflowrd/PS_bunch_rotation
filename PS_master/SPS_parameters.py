'''
SPS parameters

@author: jflowerdew
'''
import numpy as np
from scipy.constants import m_p, c, e

''' Machine parameters '''
Ekin = 26e9 # [eV]

charge = 1
E0 = m_p * c**2. / e
circumference_sps = 6911.5 # [m]
energy = Ekin + E0 # [eV]
momentum_sps = np.sqrt(energy**2. - E0**2.) # [eV/c]
beta = momentum_sps / energy
gamma = energy / E0

t_rev = circumference_sps / (beta*c)
f_rev = 1 / t_rev
w_0 = 2*np.pi*f_rev #revolution frequency

harmonic_number_sps = 4620
voltage_sps = 4.5e6 # [V]
f_rf = harmonic_number_sps*f_rev
t_rf = 1/ f_rf
phi_s = 0
phase = phi_s
w_r = harmonic_number_sps*w_0

gamma_transition_sps = 18
momentum_compaction_sps = 1 / gamma_transition_sps**2.
eta = momentum_compaction_sps - 1 / gamma**2. 