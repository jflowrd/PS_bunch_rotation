'''
Overview of all the vacuum equipment impedance

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

# Load flanges_ground

n_flanges_ground = 200

# script_path = os.path.dirname(os.path.realpath(__file__))
# toolbox_path = script_path+'/../../../impedance_toolbox'
sys.path.insert(0, script_path+'/../')
from flanges.ground_loops.circuit_model import load_model

flanges_ground = ImpedanceParameters('.')
flanges_ground.addResonators(0,1,1,
                      freqArray=freqArray)

flanges_ground.impedance += load_model(freqArray) * n_flanges_ground

# Plotting the impedance
plt.figure('Impedance')
plt.clf()
plt.plot(flanges_ground.freqArray/1e6, np.abs(flanges_ground.impedance),
         label='flanges_ground')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_vacuum.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_vacuum_zoom.png')

impedance_real_summed = flanges_ground.impedance.real

plt.figure('Impedance real summed')
plt.clf()
plt.plot(flanges_ground.freqArray/1e6, impedance_real_summed,
         label='flanges_ground')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed = flanges_ground.impedance.imag

plt.figure('Impedance imag summed')
plt.clf()
plt.plot(flanges_ground.freqArray/1e6, impedance_imag_summed/(flanges_ground.freqArray/f_rev),
         label='flanges_ground')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')

impedance_imag_summed_highfreq_only = 0


# Load flanges_gap_ps195

n_flanges_gap_ps195 = 259 - 100

filename = '../flanges/gap_ps195/cst_raw_data/eigen'
# loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

flanges_gap_ps195_loader = handleImpedance(folder=script_path)
flanges_gap_ps195_loader.importFromCST(filename, unitFreq=1e9)

flanges_gap_ps195 = ImpedanceParameters('.')
flanges_gap_ps195.addResonators(flanges_gap_ps195_loader.table_impedance[filename]['Rsh'] * n_flanges_gap_ps195,
                    flanges_gap_ps195_loader.table_impedance[filename]['fr'],
                    flanges_gap_ps195_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# flanges_gap_ps195.impedance.imag += im_z_over_f * freqArray * n_flanges_gap_ps195

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(flanges_gap_ps195.freqArray/1e6, np.abs(flanges_gap_ps195.impedance),
         label='flanges_gap_ps195')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_vacuum.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_vacuum_zoom.png')

impedance_real_summed += flanges_gap_ps195.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(flanges_gap_ps195.freqArray/1e6, impedance_real_summed,
         label='flanges_gap_ps195')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed += flanges_gap_ps195.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(flanges_gap_ps195.freqArray/1e6, impedance_imag_summed/(flanges_gap_ps195.freqArray/f_rev),
         label='flanges_gap_ps195')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')

impedance_imag_summed_highfreq_only += flanges_gap_ps195.impedance.imag

plt.figure('Impedance imag summed high freq')
# plt.clf()
plt.plot(flanges_gap_ps195.freqArray/1e6, impedance_imag_summed_highfreq_only/(flanges_gap_ps195.freqArray/f_rev),
         label='flanges_gap_ps195')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')


# Load pumping_manifold_CODD
filename = '../pumping_manifold_CODD/Resonators/multi_resonator_and_ImZ_over_f.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

n_pumping_manifold_CODD = 100-36
pumping_manifold_CODD_loader = handleImpedance(folder=script_path)
pumping_manifold_CODD_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

pumping_manifold_CODD = ImpedanceParameters('.')
pumping_manifold_CODD.addResonators(pumping_manifold_CODD_loader.table_impedance[filename]['Rsh'] * n_pumping_manifold_CODD,
                    pumping_manifold_CODD_loader.table_impedance[filename]['fr'],
                    pumping_manifold_CODD_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

pumping_manifold_CODD.impedance.imag += im_z_over_f * freqArray * n_pumping_manifold_CODD

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(pumping_manifold_CODD.freqArray/1e6, np.abs(pumping_manifold_CODD.impedance),
         label='pumping_manifold_CODD')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_vacuum.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_vacuum_zoom.png')

impedance_real_summed += pumping_manifold_CODD.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(pumping_manifold_CODD.freqArray/1e6, impedance_real_summed,
         label='pumping_manifold_CODD')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed += pumping_manifold_CODD.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(pumping_manifold_CODD.freqArray/1e6, impedance_imag_summed/(pumping_manifold_CODD.freqArray/f_rev),
         label='pumping_manifold_CODD')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')

impedance_imag_summed_highfreq_only += pumping_manifold_CODD.impedance.imag

plt.figure('Impedance imag summed high freq')
# plt.clf()
plt.plot(pumping_manifold_CODD.freqArray/1e6, impedance_imag_summed_highfreq_only/(pumping_manifold_CODD.freqArray/f_rev),
         label='pumping_manifold_CODD')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')


# Load pumping_manifold_empty
filename = '../pumping_manifold_empty/Resonators/multi_resonator_and_ImZ_over_f.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

n_pumping_manifold_empty = 36
pumping_manifold_empty_loader = handleImpedance(folder=script_path)
pumping_manifold_empty_loader.importResonatorFromFile(filename, unitFreq=1, unitRsh=1)

pumping_manifold_empty = ImpedanceParameters('.')
pumping_manifold_empty.addResonators(pumping_manifold_empty_loader.table_impedance[filename]['Rsh'] * n_pumping_manifold_empty,
                    pumping_manifold_empty_loader.table_impedance[filename]['fr'],
                    pumping_manifold_empty_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

pumping_manifold_empty.impedance.imag += im_z_over_f * freqArray * n_pumping_manifold_empty

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(pumping_manifold_empty.freqArray/1e6, np.abs(pumping_manifold_empty.impedance),
         label='pumping_manifold_empty')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_vacuum.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_vacuum_zoom.png')

impedance_real_summed += pumping_manifold_empty.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(pumping_manifold_empty.freqArray/1e6, impedance_real_summed,
         label='pumping_manifold_empty')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed += pumping_manifold_empty.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(pumping_manifold_empty.freqArray/1e6, impedance_imag_summed/(pumping_manifold_empty.freqArray/f_rev),
         label='pumping_manifold_empty')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')

impedance_imag_summed_highfreq_only += pumping_manifold_empty.impedance.imag

plt.figure('Impedance imag summed high freq')
# plt.clf()
plt.plot(pumping_manifold_empty.freqArray/1e6, impedance_imag_summed_highfreq_only/(pumping_manifold_empty.freqArray/f_rev),
         label='pumping_manifold_empty')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')

# Load bellows

n_bellows = 99

filename = '../bellows/Inductive/ImZ_over_f.txt'
loaded_impedance = np.loadtxt(script_path+'/'+filename)
im_z_over_f = loaded_impedance

bellows = ImpedanceParameters('.')
bellows.addResonators(0,1,1,
                      freqArray=freqArray)

bellows.impedance.imag += im_z_over_f * freqArray * n_bellows

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(bellows.freqArray/1e6, np.abs(bellows.impedance),
         label='bellows')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_vacuum.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_indiv_vacuum_zoom.png')

impedance_real_summed += bellows.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(bellows.freqArray/1e6, impedance_real_summed,
         label='bellows')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed.png')
# plt.xlim((0, 50))
# plt.ylim((0, 200))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed += bellows.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(bellows.freqArray/1e6, impedance_imag_summed/(bellows.freqArray/f_rev),
         label='bellows')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')

impedance_imag_summed_highfreq_only += bellows.impedance.imag

plt.figure('Impedance imag summed high freq')
# plt.clf()
plt.plot(bellows.freqArray/1e6, impedance_imag_summed_highfreq_only/(bellows.freqArray/f_rev),
         label='bellows')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
# plt.xlim((0, 50))
# plt.ylim((-1, 5))
# plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed_zoom.png')


# Load sector_valve
filename = '../sector_valve/cst_raw_data/eigen'
# loaded_impedance = np.loadtxt(script_path+'/'+filename)
# im_z_over_f = loaded_impedance[-1, -1] * loaded_impedance.shape[0]

n_sector_valve = 10
sector_valve_loader = handleImpedance(folder=script_path)
sector_valve_loader.importFromCST(filename)

sector_valve = ImpedanceParameters('.')
sector_valve.addResonators(sector_valve_loader.table_impedance[filename]['Rsh'] * n_sector_valve,
                    sector_valve_loader.table_impedance[filename]['fr'],
                    sector_valve_loader.table_impedance[filename]['Q'],
                    freqArray=freqArray)

# sector_valve.impedance.imag += im_z_over_f * freqArray * n_sector_valve

# Plotting the impedance
plt.figure('Impedance')
# plt.clf()
plt.plot(sector_valve.freqArray/1e6, np.abs(sector_valve.impedance),
         label='sector_valve')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_indiv_vacuum.png')
plt.xscale('log')
plt.xlim((0, 100))
plt.ylim((0, 5000))
plt.tight_layout()
plt.savefig(script_path+'/overview_indiv_vacuum_zoom.png')

impedance_real_summed += sector_valve.impedance.real

plt.figure('Impedance real summed')
# plt.clf()
plt.plot(sector_valve.freqArray/1e6, impedance_real_summed,
         label='sector_valve')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ReZ [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_real_summed.png')
plt.xscale('log')
plt.xlim((0, 100))
plt.ylim((0, 5000))
plt.tight_layout()
plt.savefig(script_path+'/overview_real_summed_zoom.png')

impedance_imag_summed += sector_valve.impedance.imag

plt.figure('Impedance imag summed')
# plt.clf()
plt.plot(sector_valve.freqArray/1e6, impedance_imag_summed/(sector_valve.freqArray/f_rev),
         label='sector_valve')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig(script_path+'/overview_imag_summed.png')
plt.xscale('log')
plt.xlim((0, 100))
plt.tight_layout()
plt.savefig(script_path+'/overview_imag_summed_zoom_x.png')
plt.ylim((-20, 20))
plt.tight_layout()
plt.savefig(script_path+'/overview_imag_summed_zoom_xy.png')

impedance_imag_summed_highfreq_only += sector_valve.impedance.imag

plt.figure('Impedance imag summed high freq')
# plt.clf()
plt.plot(sector_valve.freqArray/1e6, impedance_imag_summed_highfreq_only/(sector_valve.freqArray/f_rev),
         label='sector_valve')
plt.xlabel('Frequency [MHz]')
plt.ylabel('Impedance ImZ/n [$\\Omega$]')
plt.legend(loc='best')
plt.tight_layout()
# plt.savefig(script_path+'/overview_imag_summed.png')
plt.xlim((0, 50))
plt.ylim((-1, 5))
plt.tight_layout()
plt.savefig(script_path+'/overview_imag_summed_highfreqonly_zoom.png')

plt.show()
