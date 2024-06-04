'''
Exporting the C20 impedance depending on the different scenarios in measurements

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path+'/../../../../../impedance_toolbox'
sys.path.insert(0, toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters

# Summarizing all cases into lists
case_list = ['13.3MHz_100V', '13.3MHz_20kV', '20MHz_100V', '20MHz_20kV']
V_list = [100, 100, 20e3, 20e3]
fr_list = [13.3e6, 13.3e6, 20e6, 20e6]
Rs_20kVnoFB_list = [1.7e3, 1.7e3, 2e3, 2e3]
Q_20kVnoFB_list = [82.42, 82.42, 63.6, 63.6]
Q_withFB_list = [5.0222, 4.85, 4.5801, 4.2239]

# Writing cases into files and plots
indexloop = 0
for case in case_list:

    impedParams = ImpedanceParameters('.')

    fr = fr_list[indexloop]
    Q = Q_withFB_list[indexloop]
    Rs = (Rs_20kVnoFB_list[indexloop]/Q_20kVnoFB_list[indexloop]) * Q

    freqArray = np.linspace(0, 40e6, 1000)

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
    plt.savefig(script_path+'/single_resonator_'+case+'.png')

    # Saving result into text file
    R_S_save = np.array(Rs, ndmin=2)
    fr_save = np.array(fr, ndmin=2)
    Q_save = np.array(Q, ndmin=2)

    saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))

    np.savetxt(script_path+'/single_resonator_'+case+'.txt',
               saved_matrix, header='Impedance of C20 \nSingle resonator at f=%.1f MHz and V=%.1e V\nAuthor: M. Morvillo\nf_r [Hz]\tR_s [Ohm]\tQ' %(fr/1e6, V_list[indexloop]))

    indexloop += 1

