'''
Example loagin resonators

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


case_list = ['PS_SPS_ws_0deg_fork_rot_fRQRoQ_EIG.txt',
             'PS_SPS_ws_135deg_fork_rot_fRQRoQ_EIG.txt',
             'PS_SPS_wire_scanner.dat']

for case in case_list:

    loaded_impedance = handleImpedance(folder=script_path)
    
    if case == 'PS_SPS_wire_scanner.dat':
        delimiter='\t'
        unitRsh = 1e3
    else:
        delimiter=','
        unitRsh = 1
    loaded_impedance.importResonatorFromFile(case,
                                             unitFreq=1e9, unitRsh=unitRsh,
                                             delimiter=delimiter)
    
    
    exported_impedance = impedance2blond(loaded_impedance.table_impedance)
    exported_impedance.impedanceList
     
    impedParams = ImpedanceParameters('.')
     
    freArray = np.linspace(0, 2e9, 10000)
     
    impedParams.addResonators(loaded_impedance.table_impedance['%s'%(case)]['Rsh'],
                              loaded_impedance.table_impedance['%s'%(case)]['fr'],
                              loaded_impedance.table_impedance['%s'%(case)]['Q'],
                              freqArray=freArray)
    
    
    # Plotting the impedance
    plt.figure('%s'%(case))
    plt.clf()
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
    plt.savefig(script_path+'/%s.png'%(case))


plt.show()


