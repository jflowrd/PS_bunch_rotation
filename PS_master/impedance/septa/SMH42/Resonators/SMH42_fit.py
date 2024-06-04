'''
Fitting the SMH42 septum impedance with resonators

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path+'/../../../../impedance_toolbox'
sys.path.insert(0, toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters
from impedance_toolbox.handle_impedance import handleImpedance

# Loading input impedance and defining the limit in frequency for the
# multi resonators fit
septum_name = 'SMH42'
data_path = './../Measurements/wire_long_imped_Measured_SMH42_Unit_2_bumper_power_supply_sept_short.txt'
maxfev = 10000

loaded_data = np.loadtxt(data_path)

freq_data = loaded_data[:, 0]
real_data = loaded_data[:, 1] - loaded_data[:, 1][0]
imag_data = loaded_data[:, 2]

# Importing an impedance source or a full impedance model
impedParams = ImpedanceParameters('.')

impedParams.addImpedanceTable(freq_data,
                              real_data,
                              imag_data,
                              freqArrayInterp=freq_data)

# First approximation
f_start = [17e6, 26e6, 200e6]
Q_start = [10, 10, 1]
r_start = [100, 75, 50]

# Setting the fitting parameters, first run is coarse with only one resonator
fittedParameters = impedParams.fitResonators(r_start,
                                             f_start,
                                             Q_start,
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             freqBound=None,
                                             ImZoverF=2.1e-6,
                                             ImZoverFBound=None,
                                             ImZoverFScale=1e-6,
                                             fitResidue='lin_real_and_imag',
                                             frequencyWindow=[0, 150e6])

print(fittedParameters)

# Setting the fitting parameters, first run is coarse with only one resonator
fittedParameters = impedParams.fitResonators(fittedParameters[0],
                                             fittedParameters[1],
                                             fittedParameters[2],
                                             RShuntBound=0.,
                                             freqBound=0.,
                                             QBound=0.,
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             ImZoverF=fittedParameters[3],
                                             ImZoverFBound=None,
                                             ImZoverFScale=1e-6,
                                             fitResidue='lin_real_and_imag',
                                             maxiter=maxfev)

print(fittedParameters)

# Plotting the impedance
plt.figure('Impedance real imag')
plt.clf()
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.real,
         'b', alpha=0.5,
         label='Input Real')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         impedParams.resonatorForFit.impedance.real,
         'b',
         label='2 Resonators')
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.imag,
         'g', alpha=0.5,
         label='Input Imag')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         impedParams.resonatorForFit.impedance.imag,
         'g',
         label='2 Resonators')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/fitted_broadband_realimag.png')
plt.xlim((0, 50))
plt.ylim((0, 130))
plt.tight_layout()
plt.savefig(script_path+'/fitted_resonances_realimag.png')

# Saving result into text file
R_S_save = np.array(fittedParameters[0], ndmin=2)
fr_save = np.array(fittedParameters[1], ndmin=2)
Q_save = np.array(fittedParameters[2], ndmin=2)
ImZ_over_F_save = fittedParameters[3] * np.ones(R_S_save.shape) / np.prod(R_S_save.shape)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T, ImZ_over_F_save.T))

np.savetxt(script_path+'/multi_resonator_and_ImZ_over_f.txt',
           saved_matrix, header='Impedance of %s \nMulti resonator and ImZ/f fit\nThe ImZ/n was divided over the number of resonators and should be summed\nAuthor: A. Lasheen\nf_r [Hz]\t\t\tR_s [Ohm]\t\t\tQ\t\t\tImZ/f [Hz-1]' %(septum_name))


plt.show()
