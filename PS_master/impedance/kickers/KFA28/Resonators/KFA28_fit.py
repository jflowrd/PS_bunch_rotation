'''
Fitting the KFA28 kicker impedance with resonators

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
kicker_name = 'KFA28'
data_path = './../cst_raw_data/wake/'
maxfev = 10000

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

fittedParameters_BB = np.array(fittedParameters)
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
           header='Impedance of KFA28\nBroadband component\nThree resonators fit\nAuthor: A. Lasheen\nf_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ')


# Setting the fitting parameters, run only on low freq resonances
# + 1 broadband resonator

fittedParameters = impedParams.fitResonators(r_start,
                                             f_start,
                                             Q_start,
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             frequencyWindow=[0, 400e6],
                                             fitResidue='lin_real_and_imag',
                                             disp=True,
                                             method='Nelder-Mead',
                                             maxfev=maxfev)

n_res_max = 5
fittedParameters = impedParams.fitMultiResonators(fittedParameters[0],
                                                  fittedParameters[1],
                                                  fittedParameters[2],
                                                  frequencyWindow=[0, 400e6],
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
                                                  frequencyWindow=[0, 400e6],
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             fitResidue='lin_real_and_imag',
                                             disp=True,
                                             method='Nelder-Mead',
                                             maxfev=maxfev)

fittedParameters_lowfreq = np.array(fittedParameters)
print(fittedParameters)
# Plotting the result with five resonators
plt.figure('Impedance')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='4 low freq (+1 BB) Resonators')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Combining the broadband resonators with the low freq ones,
# the linear residue is taken -> better low freq fit
fittedParameters = impedParams.fitResonators(np.append(fittedParameters_lowfreq[0][1:],
                                                       fittedParameters_BB[0]),
                                             np.append(fittedParameters_lowfreq[1][1:],
                                                       fittedParameters_BB[1]),
                                             np.append(fittedParameters_lowfreq[2][1:],
                                                       fittedParameters_BB[2]),
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             fitResidue='log_real_and_imag',
                                             disp=True,
                                             method='Nelder-Mead',
                                             maxiter=maxfev)

print(fittedParameters)
# Plotting the result with extra resonators
plt.figure('Impedance')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='Combined 7 Resonators (better low freq)')
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

np.savetxt(script_path+'/resonators_broadband_extra_resonances_best_lowfreq.txt',
           saved_matrix,
           header='Impedance of KFA28\nBroadband component + resonances\nSeven resonators fit\nAuthor: A. Lasheen\nf_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ')


# Combining the broadband resonators with the low freq ones,
# the linear residue is taken -> better broadband fit
fittedParameters = impedParams.fitResonators(np.append(fittedParameters_lowfreq[0][1:],
                                                       fittedParameters_BB[0]),
                                             np.append(fittedParameters_lowfreq[1][1:],
                                                       fittedParameters_BB[1]),
                                             np.append(fittedParameters_lowfreq[2][1:],
                                                       fittedParameters_BB[2]),
                                             RShuntScale=1e2,
                                             freqScale=1e6,
                                             QScale=1,
                                             fitResidue='lin_real_and_imag',
                                             disp=True,
                                             method='Nelder-Mead',
                                             maxiter=maxfev)

print(fittedParameters)
# Plotting the result with extra resonators
plt.figure('Impedance')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='Combined 7 Resonators (better broadband)')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/fitted_broadband_extra_resonances.png')
plt.xlim((0, 300))
plt.ylim((0, 500))
plt.tight_layout()
plt.savefig(script_path+'/fitted_broadband_extra_resonances_zoom.png')

# Saving result into text file
sorted_freq = np.argsort(fittedParameters[1])

R_S_save = np.array(fittedParameters[0][sorted_freq], ndmin=2)
fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/resonators_broadband_extra_resonances_best_broadband.txt',
           saved_matrix,
           header='Impedance of KFA28\nBroadband component + resonances\nSeven resonators fit\nAuthor: A. Lasheen\nf_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ')


plt.show()
