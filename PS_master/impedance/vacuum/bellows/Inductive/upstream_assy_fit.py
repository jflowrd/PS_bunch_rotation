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

# ImZ/f fit
fittedParameters = impedParams.fitResonators([],
                                             [],
                                             [],
                                             RShuntScale=1e4,
                                             freqScale=1e9,
                                             QScale=1e3,
                                             ImZoverF=3e-8,
                                             ImZoverFScale=1e-8,
                                             fitResidue='lin_real_and_imag',
                                             method='Nelder-Mead')

print(fittedParameters)

# Plotting the result with ImZ/f
plt.figure('Impedance')
plt.plot(impedParams.resonatorForFit.frequency_array/1e6,
         np.abs(impedParams.resonatorForFit.impedance),
         label='%d Resonators' % (len(fittedParameters[0])))
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/fitted_impedance.png')

# Saving result into text file
ImZ_over_F_save = np.array(fittedParameters[3], ndmin=2)

np.savetxt(script_path+'/ImZ_over_f.txt',
           ImZ_over_F_save, header='Impedance of upstream assembly (bellow) \n' +
           'ImZ/f fit\n' +
           'Author: A. Lasheen\n' +
           'ImZ/f [Hz-1]')


plt.show()
