'''
Fitting the KFA71 kicker impedance with resonators

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
kicker_name = 'KFA71'
data_path = './../Tsutsui/'
maxfev = 10000

loaded_data = np.loadtxt(script_path+'/'+data_path+'/KFA71_Tsutsui.imp', skiprows=1)

freq_tsu = loaded_data[:, 0]
real_tsu = loaded_data[:, 1]
imag_tsu = loaded_data[:, 2]

freq_data = freq_tsu
real_data = real_tsu
imag_data = imag_tsu

# Importing an impedance source or a full impedance model
impedParams = ImpedanceParameters('.')

impedParams.addImpedanceTable(freq_data,
                              real_data,
                              imag_data,
                              freqArrayInterp=freq_data)

# Plotting the impedance
plt.figure('Impedance')
plt.clf()
plt.plot(impedParams.freqArray/1e6, np.abs(impedParams.impedance),
         label='Input')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

# Initial parameters
f_start = 650e6
Q_start = 1
r_start = 600


# Setting the fitting parameters, only broadband resonators on full range

fittedParameters = impedParams.fitResonators(r_start,
                                             f_start,
                                             Q_start,
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             fitResidue='lin_real_and_imag',
                                             disp=True,
                                             method='Nelder-Mead',
                                             maxfev=maxfev)

n_res_max = 3
fittedParameters = impedParams.fitMultiResonators(fittedParameters[0],
                                                  fittedParameters[1],
                                                  fittedParameters[2],
                                                  RShuntScale=1e2,
                                                  freqScale=1e6,
                                                  QScale=1,
                                                  fitResidue='lin_real_and_imag',
                                                  n_res_max=n_res_max,
                                                  disp=True,
                                                  method='Nelder-Mead',
                                                  maxfev=maxfev)

fittedParameters = impedParams.fitResonators(fittedParameters[0],
                                             fittedParameters[1],
                                             fittedParameters[2],
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             fitResidue='lin_real_and_imag',
                                             disp=True,
                                             method='Nelder-Mead',
                                             maxfev=maxfev)

print(fittedParameters)
# Plotting the result with three resonators
plt.figure('Impedance')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='3 Resonators')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/fitted_broadband.png')

plt.figure('Impedance real imag')
plt.clf()
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.real,
         'b', alpha=0.5,
         label='Input Real')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         impedParams.resonatorForFit.impedance.real,
         'b',
         label='3 Resonators')
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.imag,
         'g', alpha=0.5,
         label='Input Imag')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         impedParams.resonatorForFit.impedance.imag,
         'g',
         label='3 Resonators')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/fitted_broadband_realimag.png')

# Saving result into text file
sorted_freq = np.argsort(fittedParameters[1])

R_S_save = np.array(fittedParameters[0][sorted_freq], ndmin=2)
fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/resonators_broadband.txt',
           saved_matrix,
           header='Impedance of KFA71\nBroadband component\nThree resonators fit\nAuthor: A. Lasheen\nf_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ')

plt.show()
