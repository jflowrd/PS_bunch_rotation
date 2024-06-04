'''
Fitting the sector valve impedance with resonators

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import m_p, e, c


# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))

charge = 1  # e
mass = m_p*c**2/e  # eV
BField = 1.255e4  # G
BField = BField*(
            1
            - 1.441767e-6*BField
            + 2.947290e-10*BField**2.
            - 2.357026e-14*BField**3.)
bending_radius = 70.079  # m
momentum = BField * bending_radius * charge * c * 1e-4

t_rev = (2 * np.pi * 100) / (c * momentum /
                             np.sqrt(momentum ** 2. + mass ** 2.))
f_rev = 1 / t_rev

# Saving result into text file
ImZ_over_F_save = np.array(0.96/f_rev, ndmin=2)

np.savetxt(script_path+'/ImZ_over_f.txt',
           ImZ_over_F_save, header='Impedance of step transitions \n' +
           'ImZ/f [Hz-1]')


plt.show()
