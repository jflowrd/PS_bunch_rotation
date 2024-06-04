'''
Fitting the KFA04 kicker impedance with resonators

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
kicker_name = 'KFA04'
data_path = './../cst_raw_data/wake/'

loaded_imp = handleImpedance(folder=script_path)
loaded_imp.importWakeFromCST(data_path, unitFreq=1E6,
                             ZFactor=1, debug=False)

freq_data = loaded_imp.table_impedance[data_path]['fr']
real_data = loaded_imp.table_impedance[data_path]['ReZ']
imag_data = loaded_imp.table_impedance[data_path]['ImZ']

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

# Setting the fitting parameters, first run is coarse with only one resonator
# but on all the data
fittedParameters = impedParams.fitResonators(r_start,
                                             f_start,
                                             Q_start,
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             fitResidue='lin_real_and_imag',
                                             disp=True,
                                             method='Nelder-Mead')

n_res_max = 2
fittedParameters = impedParams.fitMultiResonators(fittedParameters[0],
                                                  fittedParameters[1],
                                                  fittedParameters[2],
                                                  RShuntScale=1e2,
                                                  freqScale=1e6,
                                                  QScale=1,
                                                  fitResidue='lin_real_and_imag',
                                                  n_res_max=n_res_max,
                                                  disp=True,
                                                  method='Nelder-Mead')

fittedParameters = impedParams.fitResonators(fittedParameters[0],
                                             fittedParameters[1],
                                             fittedParameters[2],
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             fitResidue='lin_real_and_imag',
                                             disp=True,
                                             method='Nelder-Mead')

print(fittedParameters)
# Plotting the result with two resonators
plt.figure('Impedance')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='2 Resonators')
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

# Saving result into text file
sorted_freq = np.argsort(fittedParameters[1])

R_S_save = np.array(fittedParameters[0][sorted_freq], ndmin=2)
fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/resonators_broadband.txt',
           saved_matrix,
           header='Impedance of KFA04\nBroadband component\nTwo resonators fit\nAuthor: A. Lasheen\nf_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ')


# Fitting extra resonances
added_shunt = [117.85871327, 197.04195815]
added_f = [6.20125891e+06, 190.065446e+06]
added_Q = [43.26241632, 389.32761719]

fittedParameters = impedParams.fitResonators(np.append(fittedParameters[0], added_shunt),
                                             np.append(fittedParameters[1], added_f),
                                             np.append(fittedParameters[2], added_Q),
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             fitResidue='lin_real_and_imag',
                                             disp=True,
                                             method='Nelder-Mead')

print(fittedParameters)
# Plotting the result with extra resonators
plt.figure('Impedance')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='4 Resonators')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/fitted_broadband_extra_resonances.png')

# Saving result into text file
sorted_freq = np.argsort(fittedParameters[1])

R_S_save = np.array(fittedParameters[0][sorted_freq], ndmin=2)
fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/resonators_broadband_extra_resonances.txt',
           saved_matrix,
           header='Impedance of KFA04\nBroadband component + resonances\nFour resonators fit\nAuthor: A. Lasheen\nf_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ')


plt.show()