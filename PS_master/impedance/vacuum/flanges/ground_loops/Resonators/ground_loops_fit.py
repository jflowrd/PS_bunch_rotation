'''
Fitting the ground loops with rf bypass with resonators

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path + '/../../../../../impedance_toolbox'
sys.path.insert(0, toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters
from impedance_toolbox.handle_impedance import handleImpedance, impedance2blond

sys.path.insert(0, script_path + '/..')
from circuit_model.circuit_model import load_model, R, L, C


n_bypasses_list = [0, 1, 2]

for n_bypasses in n_bypasses_list:

    impedance, freq_array = load_model(None, n_bypasses=n_bypasses)

    impedParams = ImpedanceParameters('.')

    impedParams.addImpedanceTable(freq_array,
                                  impedance.real,
                                  impedance.imag,
                                  freqArrayInterp=freq_array)

    # Plotting the impedance
    plt.figure('Impedance')
    plt.clf()
    plt.plot(impedParams.freqArray / 1e6, np.abs(impedParams.impedance),
             label='Input')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()

    # Getting r,f,Q
    if n_bypasses == 0:
        r_start = np.array([R])
        f_start = np.array([1 / (2 * np.pi * np.sqrt(L * C))])
        Q_start = np.array([R * np.sqrt(C / L)])
    else:
        r_start, f_start, Q_start = impedParams.fitInitialGuess(level=0.5)
    print(r_start, f_start, Q_start)

    impedParams_res = ImpedanceParameters('.')

    impedParams_res.addResonators(r_start, f_start, Q_start,
                                  freqArray=freq_array)

    plt.plot(impedParams.freqArray / 1e6, np.abs(impedParams_res.impedance),
             label='Resonators')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()

    impedParams_res = ImpedanceParameters('.')

    impedParams_res.addResonators(r_start[0], f_start[0], Q_start[0],
                                  freqArray=freq_array)

    plt.plot(impedParams.freqArray / 1e6, np.abs(impedParams_res.impedance),
             label='Resonators - Low freq')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.xlim((0.01, 200))
    plt.xscale('log')
    plt.tight_layout()
    plt.savefig(script_path + '/fitted_impedance_%d_bypasses.png' %
                (n_bypasses))

    # Saving result into text file
    sorted_freq = np.argsort(f_start)
    if n_bypasses > 1:
        R_S_save = np.array(r_start[sorted_freq][0], ndmin=2)
        fr_save = np.array(f_start[sorted_freq][0], ndmin=2)
        Q_save = np.array(Q_start[sorted_freq][0], ndmin=2)

        saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

        np.savetxt(script_path + '/multi_resonator_%d_bypass_low_freq.txt' % (n_bypasses),
                   saved_matrix, header='Impedance of ground loops with 1 RF bypass\n' +
                   'Only low freq resonance\n' +
                   'f_r [Hz]\t\t\tR_s [Ohm]\t\t\tQ')

    # Saving result into text file
    R_S_save = np.array(r_start[sorted_freq], ndmin=2)
    fr_save = np.array(f_start[sorted_freq], ndmin=2)
    Q_save = np.array(Q_start[sorted_freq], ndmin=2)

    saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

    np.savetxt(script_path + '/multi_resonator_%d_bypass.txt' % (n_bypasses),
               saved_matrix, header='Impedance of ground loops with 1 RF bypass\n' +
               'f_r [Hz]\t\t\tR_s [Ohm]\t\t\tQ')

    plt.show()
