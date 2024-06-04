'''
Fitting the sector valve impedance with resonators

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

case_list = ['wake', 'wake_spectrum']

for case in case_list:

    loaded_impedance = handleImpedance(folder=script_path)

    loaded_impedance.importWakeFromCST('../cst_raw_data/%s/'%(case))

    exported_impedance = impedance2blond(loaded_impedance.table_impedance)
    exported_impedance.impedanceList

    impedParams = ImpedanceParameters('.')

    impedParams.addImpedanceTable(loaded_impedance.table_impedance['../cst_raw_data/%s/'%(case)]['fr'],
                                  loaded_impedance.table_impedance['../cst_raw_data/%s/'%(case)]['ReZ'],
                                  loaded_impedance.table_impedance['../cst_raw_data/%s/'%(case)]['ImZ'],
                                  freqArrayInterp=loaded_impedance.table_impedance['../cst_raw_data/%s/'%(case)]['fr'])

    # Plotting the impedance
    plt.figure('Impedance %s'%(case))
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
    fittedParameters_main = impedParams.fitResonators(r_start,
                                                 f_start,
                                                 Q_start,
                                                 RShuntScale=1e4,
                                                 freqScale=1e9,
                                                 QScale=1e3,
                                                 fitResidue='log_real_and_imag')

    fittedParameters = impedParams.fitResonators([fittedParameters_main[0],
                                                  50.,
                                                  100.,
                                                  50.,
                                                  1e4],
                                                 [fittedParameters_main[1],
                                                  495e6,
                                                  1086e6,
                                                  1730e6,
                                                  3000e6],
                                                 [fittedParameters_main[2],
                                                  200.,
                                                  200.,
                                                  200.,
                                                  150.],
                                                 RShuntScale=1e4,
                                                 freqScale=1e9,
                                                 QScale=1e3,
                                                 fitResidue='log_total',
                                                 maxiter=10000,
                                                 method='Nelder-Mead')

    # Second iteration fine tuning
    fittedParameters = impedParams.fitResonators(fittedParameters[0],
                                                 fittedParameters[1],
                                                 fittedParameters[2],
                                                 RShuntScale=1e4,
                                                 freqScale=1e9,
                                                 QScale=1e3,
                                                 fitResidue='lin_total',
                                                 maxiter=10000,
                                                 method='Nelder-Mead')

    print(fittedParameters)

    # Plotting the result with multiple resonators
    plt.figure('Impedance %s'%(case))
    plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
             np.abs(impedParams.resonatorForFit.impedance),
             label='%d Resonators' % (len(fittedParameters[0])))
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.yscale("log", nonposy='clip')
    plt.ylim((1, 1e5))
    plt.tight_layout()
    plt.savefig(script_path+'/fitted_impedance_%s.png' % (case))

    # Saving result into text file
    sorted_freq = np.argsort(fittedParameters[1])
    R_S_save = np.array(fittedParameters[0][sorted_freq], ndmin=2)
    fr_save = np.array(fittedParameters[1][sorted_freq], ndmin=2)
    Q_save = np.array(fittedParameters[2][sorted_freq], ndmin=2)

    saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

    np.savetxt(script_path+'/multi_resonator_%s.txt' % (case),
               saved_matrix,
               header='Impedance of Sector valves \n' +
               'Multi resonator\n' +
               'Author: A. Lasheen\n' +
               'f_r [Hz]\t\t\tR_s [Ohm]\t\t\tQ')

plt.show()
