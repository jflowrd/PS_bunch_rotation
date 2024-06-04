'''
Fitting the C10-36 cavity impedance with resonators

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
from impedance_toolbox.handle_impedance import handleImpedance, impedance2blond

loaded_impedance = handleImpedance(folder=script_path)

loaded_impedance.importWakeFromCST('../cst_raw_data/wake/')

exported_impedance = impedance2blond(loaded_impedance.table_impedance)
exported_impedance.impedanceList

impedParams = ImpedanceParameters('.')

impedParams.addImpedanceTable(loaded_impedance.table_impedance['../cst_raw_data/wake/']['fr'],
                              loaded_impedance.table_impedance['../cst_raw_data/wake/']['ReZ'],
                              loaded_impedance.table_impedance['../cst_raw_data/wake/']['ImZ'],
                              freqArrayInterp=loaded_impedance.table_impedance['../cst_raw_data/wake/']['fr'])


# Plotting the impedance
plt.figure('Impedance')
plt.clf()
plt.plot(impedParams.freqArray/1e6, np.abs(impedParams.impedance),
         label='Input')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Getting initial guess
r_start, f_start, Q_start = impedParams.fitInitialGuess(level=0.9)


# Fitting with multiple resonators
n_res_max = 1
fittedParameters = impedParams.fitMultiResonators(r_start,
                                                  f_start,
                                                  Q_start,
                                                  RShuntScale=1e3,
                                                  freqScale=1e6,
                                                  QScale=1,
                                                  fitResidue='lin_total',
                                                  n_res_max=n_res_max,
                                                  method='Nelder-Mead')

# Adding manually missing resonances
fittedParameters = impedParams.fitResonators(np.append(fittedParameters[0], [10, 14.5640832]),
                                             np.append(fittedParameters[1], [0.8e9, 1.67171093e+09]),
                                             np.append(fittedParameters[2], [1, 1]),
                                             RShuntScale=1e3,
                                             freqScale=1e6,
                                             QScale=1,
                                             fitResidue='lin_total',
                                             method='Nelder-Mead')

# Second iteration fine tuning
fittedParameters = impedParams.fitResonators(fittedParameters[0],
                                             fittedParameters[1],
                                             fittedParameters[2],
                                             RShuntScale=1e3,
                                             freqScale=1e6,
                                             QScale=1,
                                             fitResidue='lin_real_and_imag',
                                             method='Nelder-Mead')

# Forcing Q values >= 0.5
fittedParameters = impedParams.fitResonators(fittedParameters[0],
                                             fittedParameters[1],
                                             fittedParameters[2],
                                             RShuntScale=1e3,
                                             freqScale=1e6,
                                             QScale=1,
                                             RShuntBound=0.,
                                             freqBound=0.,
                                             QBound=[0.5, 1.0],
                                             fitResidue='lin_real_and_imag',
                                             maxiter=10000)

print(fittedParameters)

# Plotting the result with multiple resonators
plt.figure('Impedance')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='%d Resonators' % (len(fittedParameters[0])))
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.yscale("log", nonposy='clip')
plt.ylim((1, 1e3))
plt.tight_layout()
plt.savefig(script_path+'/fitted_impedance.png')
plt.xlim((0, 50))
plt.tight_layout()
plt.savefig(script_path+'/fitted_impedance_zoom.png')

plt.figure('Impedance real imag')
plt.clf()
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.real,
         'b', alpha=0.5,
         label='Input Real')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         impedParams.resonatorForFit.impedance.real,
         'b',
         label='Resonators')
plt.plot(impedParams.freqArray/1e6, impedParams.impedance.imag,
         'g', alpha=0.5,
         label='Input Imag')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         impedParams.resonatorForFit.impedance.imag,
         'g',
         label='Resonators')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/fitted_impedance_realimag.png')

# Saving result into text file
sorted_freq = np.argsort(fittedParameters[1])
R_S_save = np.array(fittedParameters[0][sorted_freq], ndmin=2)
fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)

saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

np.savetxt(script_path+'/multi_resonator.txt',
           saved_matrix, header='Impedance of Wall Current Monitor\n' +
           'Multi resonator fit\n' +
           'Author: A. Lasheen' +
           '\nf_r [Hz]\t\t\tR_s [Ohm]\t\t\tQ')

plt.show()
