'''
Overview of all the misc impedance

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

# Load internal_dump

n_internal_dump = 2

filename = '../internal_dump/cst_raw_data/eigen'
# loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

internal_dump_loader = handleImpedance(folder=script_path)
internal_dump_loader.importFromCST(filename)

internal_dump = ImpedanceParameters('.')
internal_dump.addResonators(internal_dump_loader.table_impedance[filename]['Rsh'] * n_internal_dump,
                    internal_dump_loader.table_impedance[filename]['fr'],
                    internal_dump_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# internal_dump.impedance.imag += im_z_over_f * freqArray * n_internal_dump


# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(internal_dump.freqArray/1e6, np.abs(internal_dump.impedance),
         label='internal_dump')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_misc.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_misc_zoom.png')

impedance_real_summed = internal_dump.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(internal_dump.freqArray/1e6, impedance_real_summed,
         label='internal_dump')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed = internal_dump.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(internal_dump.freqArray/1e6, impedance_imag_summed/(internal_dump.freqArray/f_rev),
         label='internal_dump')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')


# Load ralentisseur
filename = '../ralentisseur/Resonators/cst_eig.txt'
# loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

n_ralentisseur = 1
ralentisseur_loader = handleImpedance(folder=script_path)
ralentisseur_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1,
                                            delimiter=',')

ralentisseur = ImpedanceParameters('.')
ralentisseur.addResonators(ralentisseur_loader.table_impedance[filename]['Rsh'] * n_ralentisseur,
                    ralentisseur_loader.table_impedance[filename]['fr'],
                    ralentisseur_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# ralentisseur.impedance.imag += im_z_over_f * freqArray * n_ralentisseur

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(ralentisseur.freqArray/1e6, np.abs(ralentisseur.impedance),
         label='ralentisseur')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_misc.png')
# plt.xscale('log')
# plt.xlim((0, 100))
# plt.ylim((0, 5000))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_misc_zoom.png')

impedance_real_summed += ralentisseur.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(ralentisseur.freqArray/1e6, impedance_real_summed,
         label='ralentisseur')
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

impedance_imag_summed += ralentisseur.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(ralentisseur.freqArray/1e6, impedance_imag_summed/(ralentisseur.freqArray/f_rev),
         label='ralentisseur')
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


# Load TFB_Kicker_SD97
filename = '../TFB_Kicker_SD97/Resonators/multi_resonator_and_ImZ_over_f.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

n_TFB_Kicker_SD97 = 1
TFB_Kicker_SD97_loader = handleImpedance(folder=script_path)
TFB_Kicker_SD97_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

TFB_Kicker_SD97 = ImpedanceParameters('.')
TFB_Kicker_SD97.addResonators(TFB_Kicker_SD97_loader.table_impedance[filename]['Rsh'] * n_TFB_Kicker_SD97,
                    TFB_Kicker_SD97_loader.table_impedance[filename]['fr'],
                    TFB_Kicker_SD97_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

TFB_Kicker_SD97.impedance.imag += im_z_over_f * freqArray * n_TFB_Kicker_SD97
  
# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(TFB_Kicker_SD97.freqArray/1e6, np.abs(TFB_Kicker_SD97.impedance),
         label='TFB_Kicker_SD97')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_indiv_misc.png')
  
impedance_real_summed += TFB_Kicker_SD97.impedance.real
  
plt.figure('Impedance real summed')
# plt.clf()
plt.plot(TFB_Kicker_SD97.freqArray/1e6, impedance_real_summed,
         label='TFB_Kicker_SD97')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_real_summed.png')
  
impedance_imag_summed += TFB_Kicker_SD97.impedance.imag
  
plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(TFB_Kicker_SD97.freqArray/1e6, impedance_imag_summed/(TFB_Kicker_SD97.freqArray/f_rev),
         label='TFB_Kicker_SD97')
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
