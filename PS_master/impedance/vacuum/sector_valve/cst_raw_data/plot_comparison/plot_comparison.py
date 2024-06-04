'''
Comparing the wake/eigen impedance

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
from impedance_toolbox.handle_impedance import handleImpedance, impedance2blond

case_list = ['wake', 'wake_spectrum', 'eigen']

for case in case_list:

    loaded_impedance = handleImpedance(folder=script_path)

    if case == 'eigen':
        loaded_impedance.importEigenFromCST('../%s/'%(case))

        exported_impedance = impedance2blond(loaded_impedance.table_impedance)
        exported_impedance.impedanceList

        impedParams = ImpedanceParameters('.')

        freArray = np.linspace(0, 3e9, 10000)

        impedParams.addResonators(loaded_impedance.table_impedance['../%s/'%(case)]['Rsh'],
                                  loaded_impedance.table_impedance['../%s/'%(case)]['fr'],
                                  loaded_impedance.table_impedance['../%s/'%(case)]['Q'],
                                  freqArray=freArray)
    else:
        loaded_impedance.importWakeFromCST('../%s/'%(case))

        exported_impedance = impedance2blond(loaded_impedance.table_impedance)
        exported_impedance.impedanceList

        impedParams = ImpedanceParameters('.')

        impedParams.addImpedanceTable(loaded_impedance.table_impedance['../%s/'%(case)]['fr'],
                                      loaded_impedance.table_impedance['../%s/'%(case)]['ReZ'],
                                      loaded_impedance.table_impedance['../%s/'%(case)]['ImZ'],
                                      freqArrayInterp=loaded_impedance.table_impedance['../%s/'%(case)]['fr'])

    # Plotting the impedance
    plt.figure('Impedance')
#     plt.clf()
    plt.plot(impedParams.freqArray/1e6, np.abs(impedParams.impedance),
             label=case)
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Impedance [$\\Omega$]')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.legend(loc='best')
    plt.yscale("log", nonposy='clip')
    plt.ylim((1, 1e5))
    plt.tight_layout()
    plt.savefig(script_path+'/cst_comparison.png')


plt.show()
