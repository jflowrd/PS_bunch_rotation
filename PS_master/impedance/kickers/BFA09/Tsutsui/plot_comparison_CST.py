'''
Plotting the Tsutsui model as simulated in
https://gitlab.cern.ch/IRIS/PS_IW_model/-/tree/master/Impedances/Longitudinal
compared to the CST simulation

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

# Loading the Tsutsui model

loaded_data_1 = np.loadtxt(script_path+'/BFA09P_Tsutsui.imp', skiprows=1)

freq_tsu_1 = loaded_data_1[:, 0]
real_tsu_1 = loaded_data_1[:, 1]
imag_tsu_1 = loaded_data_1[:, 2]
imp_tsu_1 = real_tsu_1 + 1j*imag_tsu_1

loaded_data_2 = np.loadtxt(script_path+'/BFA09S_Tsutsui.imp', skiprows=1)

freq_tsu_2 = loaded_data_2[:, 0]
real_tsu_2 = loaded_data_2[:, 1]
imag_tsu_2 = loaded_data_2[:, 2]
imp_tsu_2 = real_tsu_2 + 1j*imag_tsu_2

freq_tsu = freq_tsu_1
imp_tsu = imp_tsu_1 + imp_tsu_2

# Plotting the impedance
plt.figure('Impedance')
plt.clf()
plt.plot(freq_tsu/1e6, np.abs(imp_tsu_1), '--',
         label='Tsutsui P')
plt.plot(freq_tsu/1e6, np.abs(imp_tsu_2), '--',
         label='Tsutsui S')
plt.plot(freq_tsu/1e6, np.abs(imp_tsu),
         label='Tsutsui P+S')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()


# Loading input impedance and defining the limit in frequency for the
# multi resonators fit
kicker_name = 'BFA09'
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
plt.plot(impedParams.freqArray/1e6, np.abs(impedParams.impedance),
         label='CST P+S')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.xlim((0, impedParams.freqArray[-1]/1e6))
plt.tight_layout()
plt.savefig(script_path+'/comparison_tsu_cst.png')

plt.show()
