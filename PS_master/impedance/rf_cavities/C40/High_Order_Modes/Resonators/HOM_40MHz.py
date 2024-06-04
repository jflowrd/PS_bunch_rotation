'''
Outputting the HOM impedance for C40

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

fr = 1e6*np.array([260.747, 284.0, 394.367, 444.89, 503.347, 557.903, 611.133, 655.308, 732.224, 774.713, 813.656, 862.611, 876.039])
Q = np.array([50, 10, 410, 10, 10, 6000, 410, 460, 10000, 4900, 100, 1800, 10])
Rs_over_Q = np.array([1.29, 0.01, 5.33, 0.92, 5.82, 0.01, 1.08, 0.51, 0.85, 1.36, 0.01, 2.15, 0.79])
Rs = Rs_over_Q*Q

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
plt.savefig(script_path+'/HOMs_C40.png')

# Saving result into text file
R_S_save = np.array(Rs, ndmin=2)
fr_save = np.array(fr, ndmin=2)
Q_save = np.array(Q, ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/HOMs_C40.txt',
           saved_matrix, header='Impedance of C40 HOMs for one cavity\nBased on simulated R/Q and measured Q_loaded\nf_r [Hz]\tR_s [Ohm]\tQ')

