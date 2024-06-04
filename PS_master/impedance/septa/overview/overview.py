'''
Overview of all the septa impedance

@author: alasheen
'''

# General import
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Impedance tools imports
script_path = os.path.dirname(os.path.realpath(__file__))
toolbox_path = script_path+'/../../../impedance_toolbox'
sys.path.insert(0, toolbox_path)
from impedance_toolbox.impedance_params import ImpedanceParameters
from impedance_toolbox.handle_impedance import handleImpedance


freqArray = np.linspace(0, 2e9, 10000)
f_rev = 477e3

# Load SMH42
filename = '../SMH42/Resonators/multi_resonator_and_ImZ_over_f.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

smh42_loader = handleImpedance(folder=script_path)
smh42_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

smh42 = ImpedanceParameters('.')
smh42.addResonators(smh42_loader.table_impedance[filename]['Rsh'],
                    smh42_loader.table_impedance[filename]['fr'],
                    smh42_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

smh42.impedance.imag += im_z_over_f * freqArray

# Plotting the impedance
plt.figure('Impedance')
plt.clf()
plt.plot(smh42.freqArray/1e6, np.abs(smh42.impedance),
         label='SMH42')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_indiv_septa.png')
plt.xlim((0, 50))
plt.ylim((0, 200))
plt.tight_layout()
plt.savefig(script_path+'/overview_indiv_septa_zoom.png')

impedance_real_summed = smh42.impedance.real

plt.figure('Impedance real summed')
plt.clf()
plt.plot(smh42.freqArray/1e6, impedance_real_summed,
         label='SMH42')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_real_summed.png')
plt.xlim((0, 50))
plt.ylim((0, 200))
plt.tight_layout()
plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed = smh42.impedance.imag

plt.figure('Impedance imag summed')
plt.clf()
plt.plot(smh42.freqArray/1e6, impedance_imag_summed/(smh42.freqArray/f_rev),
         label='SMH42')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_imag_summed.png')
plt.xlim((0, 50))
plt.ylim((-1, 5))
plt.tight_layout()
plt.savefig(script_path+'/overview_imag_summed_zoom.png')

plt.show()
