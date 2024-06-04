'''
Outputting the HOM impedance for C80

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
from scipy.constants import m_p, c, e
import matplotlib.pyplot as plt

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path+'/../../../../../impedance_toolbox'
sys.path.insert(0, toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters


# Writing cases into files and plots
impedParams = ImpedanceParameters('.')

fr = 1e6*np.array([255.9, 340.01, 436.40, 538.0, 559.13, 659.56, 702.92, 760.55, 806.65, 848.32, 912.97, 948.23, 964.31, 1056.22])
Q = np.array([75, 90, 215, 480, 80, 3400, 470, 10, 1440, 1400, 1940, 10000, 8700, 2100])
Rs = 1e3*np.array([1.5, 0.54, 1.03, 0.18, 0.36, 0.47, 2.80, 0.01, 1.07, 3.44, 3.65, 9.21, 28.85, 1.16])


freqArray = np.linspace(0, 1.2e9, 100000)

impedParams.addResonators(Rs, fr, Q, freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.clf()
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.real,
         label='Real')
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.imag,
         label='Imag')
plt.plot(impedParams.freqArray/1e6, np.abs(impedParams.impedance),
         label='Abs')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/HOMs_C80.png')

# Saving result into text file
R_S_save = np.array(Rs, ndmin=2)
fr_save = np.array(fr, ndmin=2)
Q_save = np.array(Q, ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/HOMs_C80.txt',
           saved_matrix, header='Impedance of C80 HOMs for one cavity\nBased on simulated R/Q and measured Q_loaded with 4 couplers\nf_r [Hz]\tR_s [Ohm]\tQ')

