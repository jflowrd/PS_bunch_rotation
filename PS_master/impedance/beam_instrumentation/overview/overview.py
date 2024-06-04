'''
Overview of all the BI equipment impedance

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


freqArray = np.linspace(0, 3e9, 1000000)
f_rev = 477e3

# Load WCM_SD03

n_WCM_SD03 = 2

sys.path.insert(0, script_path+'/../')
from WCM_SD03.circuit_model.circuit_model import load_model

WCM_SD03 = ImpedanceParameters('.')
WCM_SD03.addResonators(0,1,1,
                      freqArray=freqArray)

WCM_SD03.impedance += load_model(freqArray)[0] * n_WCM_SD03


# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(WCM_SD03.freqArray/1e6, np.abs(WCM_SD03.impedance),
         label='WCM_SD03')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_BI.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_BI_zoom.png')

impedance_real_summed = WCM_SD03.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(WCM_SD03.freqArray/1e6, impedance_real_summed,
         label='WCM_SD03')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed = WCM_SD03.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(WCM_SD03.freqArray/1e6, impedance_imag_summed/(WCM_SD03.freqArray/f_rev),
         label='WCM_SD03')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')


# Load Stripline_BPM_SD72
filename = '../Stripline_BPM_SD72/Resonators/multi_resonator.txt'
# loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

n_Stripline_BPM_SD72 = 1
Stripline_BPM_SD72_loader = handleImpedance(folder=script_path)
Stripline_BPM_SD72_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

Stripline_BPM_SD72 = ImpedanceParameters('.')
Stripline_BPM_SD72.addResonators(Stripline_BPM_SD72_loader.table_impedance[filename]['Rsh'] * n_Stripline_BPM_SD72,
                    Stripline_BPM_SD72_loader.table_impedance[filename]['fr'],
                    Stripline_BPM_SD72_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# Stripline_BPM_SD72.impedance.imag += im_z_over_f * freqArray * n_Stripline_BPM_SD72

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(Stripline_BPM_SD72.freqArray/1e6, np.abs(Stripline_BPM_SD72.impedance),
         label='Stripline_BPM_SD72')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_BI.png')
# plt.xscale('log')
# plt.xlim((0, 100))
# plt.ylim((0, 5000))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_BI_zoom.png')

impedance_real_summed += Stripline_BPM_SD72.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(Stripline_BPM_SD72.freqArray/1e6, impedance_real_summed,
         label='Stripline_BPM_SD72')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed.png')
# plt.xscale('log')
# plt.xlim((0, 100))
# plt.ylim((0, 5000))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed += Stripline_BPM_SD72.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(Stripline_BPM_SD72.freqArray/1e6, impedance_imag_summed/(Stripline_BPM_SD72.freqArray/f_rev),
         label='Stripline_BPM_SD72')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xscale('log')
# plt.xlim((0, 100))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom_x.png')
# plt.ylim((-20, 20))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom_xy.png')



# Load BGI_Vertical
filename = '../BGI/Vertical/Resonators/multi_resonator.txt'
# loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

n_BGI_Vertical = 1
BGI_Vertical_loader = handleImpedance(folder=script_path)
BGI_Vertical_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

BGI_Vertical = ImpedanceParameters('.')
BGI_Vertical.addResonators(BGI_Vertical_loader.table_impedance[filename]['Rsh'] * n_BGI_Vertical,
                    BGI_Vertical_loader.table_impedance[filename]['fr'],
                    BGI_Vertical_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# BGI_Vertical.impedance.imag += im_z_over_f * freqArray * n_BGI_Vertical

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(BGI_Vertical.freqArray/1e6, np.abs(BGI_Vertical.impedance),
         label='BGI_Vertical')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_BI.png')
# plt.xscale('log')
# plt.xlim((0, 100))
# plt.ylim((0, 5000))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_BI_zoom.png')

impedance_real_summed += BGI_Vertical.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(BGI_Vertical.freqArray/1e6, impedance_real_summed,
         label='BGI_Vertical')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed.png')
# plt.xscale('log')
# plt.xlim((0, 100))
# plt.ylim((0, 5000))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed += BGI_Vertical.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(BGI_Vertical.freqArray/1e6, impedance_imag_summed/(BGI_Vertical.freqArray/f_rev),
         label='BGI_Vertical')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xscale('log')
# plt.xlim((0, 100))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom_x.png')
# plt.ylim((-20, 20))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom_xy.png')


# Load WS
filename = '../WS/cst_raw_data/PS_SPS_wire_scanner.dat'

n_WS = 4
WS_loader = handleImpedance(folder=script_path)
WS_loader.importResonatorFromFile(filename, unitFreq=1e9, unitRsh=1e3)

WS = ImpedanceParameters('.')
WS.addResonators(WS_loader.table_impedance[filename]['Rsh'] * n_WS,
                    WS_loader.table_impedance[filename]['fr'],
                    WS_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# WS.impedance.imag += im_z_over_f * freqArray * n_WS

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(WS.freqArray/1e6, np.abs(WS.impedance),
         label='WS')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_indiv_BI.png')

impedance_real_summed += WS.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(WS.freqArray/1e6, impedance_real_summed,
         label='WS')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_real_summed.png')

impedance_imag_summed += WS.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(WS.freqArray/1e6, impedance_imag_summed/(WS.freqArray/f_rev),
         label='WS')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_imag_summed.png')
plt.xlim((0, 50))
plt.ylim((0, 1))
plt.tight_layout()
plt.savefig(script_path+'/overview_imag_summed_zoom.png')

plt.show()
