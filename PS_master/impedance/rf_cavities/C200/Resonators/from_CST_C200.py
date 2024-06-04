'''
Outputting the design impedance for C200
'''

# General import
import os
import sys
import numpy as np
from scipy.constants import m_p, c, e
import matplotlib.pyplot as plt

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path+'/../../../../impedance_toolbox'
sys.path.insert(0, toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters


# Writing cases into files and plots
impedParams = ImpedanceParameters('.')

charge = 1  # e
mass = m_p*c**2./e  # eV
BField = 0.66681  # T
bending_radius = 70.079  # m
momentum = BField*bending_radius*charge*c
print(momentum)

t_rev = (2*np.pi*100)/(c*momentum/np.sqrt(momentum**2.+mass**2.))
f_rev = 1/t_rev

fr = 199.948e6  # 420*f_rev
Q = 134.
Rs_over_Q = 28.5
Rs = Rs_over_Q*Q

freqArray = np.linspace(0, 400e6, 100000)

impedParams.addResonators(Rs, fr, Q, freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.clf()
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.real/1e3,
         label='Real')
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.imag/1e3,
         label='Imag')
plt.plot(impedParams.freqArray/1e6, np.abs(impedParams.impedance)/1e3,
         label='Abs')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$k\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/from_CST_C200_single_damped.png')

# Saving result into text file
R_S_save = np.array(Rs, ndmin=2)
fr_save = np.array(fr, ndmin=2)
Q_save = np.array(Q, ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/from_CST_C200_single_damped.txt',
           saved_matrix, header='Impedance as measured by S. Persichelli \n' +
           'Q damped by 3 PIN lines (rf not in use)\n' +
           'Single gap, 6 in total\n' +
           'Single resonator\n' +
           'f_r [Hz]\tR_s [Ohm]\tQ')


# Writing cases into files and plots
impedParams = ImpedanceParameters('.')

Q = 971.
Rs = Rs_over_Q*Q

freqArray = np.linspace(0, 400e6, 100000)

impedParams.addResonators(Rs, fr, Q, freqArray=freqArray)

# Plotting the impedance
plt.figure('Impedance')
plt.clf()
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.real/1e3,
         label='Real')
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.imag/1e3,
         label='Imag')
plt.plot(impedParams.freqArray/1e6, np.abs(impedParams.impedance)/1e3,
         label='Abs')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$k\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/from_CST_C200_single_operational.png')

# Saving result into text file
R_S_save = np.array(Rs, ndmin=2)
fr_save = np.array(fr, ndmin=2)
Q_save = np.array(Q, ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/from_CST_C200_single_operational.txt',
           saved_matrix, header='Impedance as measured by S. Persichelli \n' +
           'Q undamped (rf in use)\n' +
           'Single resonator\n' +
           'f_r [Hz]\tR_s [Ohm]\tQ')
