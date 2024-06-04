'''
Fitting the stripline bpm impedance with resonators

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path + '/../../../../impedance_toolbox'
sys.path.insert(0, toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters
from impedance_toolbox.handle_impedance import handleImpedance, impedance2blond


case = 'wake'

loaded_impedance = handleImpedance(folder=script_path)

loaded_impedance.importWakeFromCST('../cst_raw_data/%s/' % (case))

exported_impedance = impedance2blond(loaded_impedance.table_impedance)
exported_impedance.impedanceList

impedParams = ImpedanceParameters('.')

# Frequency array not regular, reinterpolating!

freq_interp = np.linspace(loaded_impedance.table_impedance['../cst_raw_data/%s/' % (case)]['fr'][0],
                          loaded_impedance.table_impedance['../cst_raw_data/%s/' % (
                              case)]['fr'][-1],
                          len(loaded_impedance.table_impedance['../cst_raw_data/%s/' % (case)]['fr']) * 3)
real_interp = np.interp(freq_interp,
                        loaded_impedance.table_impedance['../cst_raw_data/%s/' % (
                            case)]['fr'],
                        loaded_impedance.table_impedance['../cst_raw_data/%s/' % (case)]['ReZ'] * 3)
imag_interp = np.interp(freq_interp,
                        loaded_impedance.table_impedance['../cst_raw_data/%s/' % (
                            case)]['fr'],
                        loaded_impedance.table_impedance['../cst_raw_data/%s/' % (case)]['ImZ'] * 3)

impedParams.addImpedanceTable(freq_interp,
                              real_interp,
                              imag_interp,
                              freqArrayInterp=freq_interp)

# Plotting the impedance
plt.figure('Impedance %s' % (case))
plt.clf()
plt.plot(impedParams.freqArray / 1e6, np.abs(impedParams.impedance),
         label='Input')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()

# Getting initial guess
r_start, f_start, Q_start = impedParams.fitInitialGuess(level=0.9)
r_start, f_start, Q_start = impedParams.fitInitialGuess(level=4500 / r_start)

# Q_start = 1e3

print(r_start, f_start, Q_start)

# Fitting with multiple resonators (first resonator changed numerical
# issue (noise))
r_start[0] = 2610
f_start[0] = 1.256e9
Q_start[0] = 300
fittedParameters = impedParams.fitResonators(r_start,
                                             f_start,
                                             Q_start,
                                             freqBound=0.,
                                             RShuntBound=0.,
                                             RShuntScale=1e4,
                                             freqScale=1e9,
                                             QScale=1e3,
                                             ImZoverF=1e-8,
                                             ImZoverFScale=1e-8,
                                             fitResidue='lin_real_and_imag',
                                             maxiter=10000)


# Second iteration fine tuning
fittedParameters = impedParams.fitResonators(fittedParameters[0],
                                             fittedParameters[1],
                                             fittedParameters[2],
                                             RShuntScale=1e4,
                                             freqScale=1e9,
                                             QScale=1e3,
                                             ImZoverF=fittedParameters[-1],
                                             ImZoverFScale=1e-8,
                                             fitResidue='lin_real_and_imag',
                                             maxiter=10000)

print(fittedParameters)

# Plotting the result with multiple resonators
plt.figure('Impedance %s' % (case))
plt.plot(impedParams.resonatorForFit.frequency_array / 1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='%d Resonators' % (len(fittedParameters[0])))
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.yscale("log", nonposy='clip')
plt.ylim((1, 1e5))
plt.tight_layout()
plt.savefig(script_path + '/fitted_impedance_%s.png' % (case))

# Saving result into text file
sorted_freq = np.argsort(fittedParameters[1])
R_S_save = np.array(fittedParameters[0][sorted_freq], ndmin=2)
fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)
ImZ_over_F_save = fittedParameters[3] * \
    np.ones(R_S_save.shape) / np.prod(R_S_save.shape)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T, ImZ_over_F_save.T))

np.savetxt(script_path + '/multi_resonator_and_ImZ_over_f.txt',
           saved_matrix, header='Impedance of Stripline BPM\n' +
           'Multi resonator and ImZ/f fit\n' +
           'The total ImZ/f=%.4e [Hz-1]\n' % (fittedParameters[3]) +
           '(The ImZ/f was divided over the number of resonators and should be summed)\n' +
           'Author: A. Lasheen\nf_r [Hz]\t\t\tR_s [Ohm]\t\t\tQ\t\t\tImZ/f [Hz-1]')

plt.show()
