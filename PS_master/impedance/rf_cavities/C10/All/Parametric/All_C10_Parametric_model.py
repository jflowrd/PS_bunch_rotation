'''
Created on 14 feb. 2018

@author: alasheen
'''

# General imports
import os
import numpy as np
import matplotlib.pyplot as plt

# Toolbox import
script_path = os.path.dirname(os.path.realpath(__file__))


def load_model(frequency):
    '''
    Method to load the parametric model
    Requires the rf frequency as an argument
    '''

    loaded_impedance = np.loadtxt(
        script_path+'/single_resonator_param_model.txt')

    R_S_interp = np.interp(frequency,
                           loaded_impedance[:, 0],
                           loaded_impedance[:, 1])
    Q_interp = np.interp(frequency,
                         loaded_impedance[:, 0],
                         loaded_impedance[:, 2])
    f_r_interp = frequency

    return f_r_interp, R_S_interp, Q_interp


if __name__ == '__main__':

    # Loading resonator fits and fitting evolution of parameters as a function
    # of the cavity harmonic

    harmonic_list = [8, 16, 21]

    R_S = np.zeros(len(harmonic_list))
    f_r = np.zeros(len(harmonic_list))
    Q = np.zeros(len(harmonic_list))

    indexloop = 0
    for harmonic in harmonic_list:

        loaded_file = np.loadtxt(script_path +
                                 '/../Resonators/single_resonator_h' +
                                 str(harmonic) +
                                 '.txt')

        R_S[indexloop] = loaded_file[1]
        f_r[indexloop] = loaded_file[0]
        Q[indexloop] = loaded_file[2]

        indexloop += 1

    sorted_freq = np.argsort(f_r)

    R_S_save = np.array(R_S[sorted_freq], ndmin=2)
    fr_save = np.array(f_r[sorted_freq], ndmin=2)
    Q_save = np.array(Q[sorted_freq], ndmin=2)

    saved_matrix = np.hstack((fr_save.T, R_S_save.T, Q_save.T))
    np.savetxt(script_path+'/single_resonator_param_model.txt',
               saved_matrix, header='Impedance of all cavities except C10-11' +
               '\nTotal impedance divided by 10 to have "single" cavity impedance' +
               '\nParametric model based on single resonators fit' +
               '\nAuthor: A. Lasheen\nf_r [Hz]\t\t\t\tR_s [Ohm]\t\t\t\tQ')

    plt.figure('Impedance all')
    plt.clf()
    plt.plot(f_r/1e6, R_S/1e3, 'bo--')
    plt.ylim((0, 1))
    plt.xlabel('Resonant frequency $f_r$ [MHz]')
    plt.ylabel('Shunt impedance $R_s$ [$\\mathrm{k\\Omega}$]')
    plt.twinx()
    plt.plot(f_r/1e6, Q, 'go--')
    plt.ylabel('Q')
    plt.ylim((0, 10))
    plt.grid('on')
    plt.tight_layout()
    plt.savefig(script_path+'/all_cavities.png')

    # Taking all the cavities one by one
    cavity_list = ['C10-36', 'C10-46', 'C10-51', 'C10-56', 'C10-66', 'C10-76',
                   'C10-81', 'C10-86', 'C10-91', 'C10-96']
    harmonic_list = [8, 16, 21]

    R_S = np.zeros((len(cavity_list), len(harmonic_list)))
    f_r = np.zeros((len(cavity_list), len(harmonic_list)))
    Q = np.zeros((len(cavity_list), len(harmonic_list)))

    plt.figure('Impedance single Rs')
    plt.clf()
    plt.figure('Impedance single Q')
    plt.clf()

    indexloop_h = 0
    for harmonic in harmonic_list:

        indexloop_c = 0
        for cavity_name in cavity_list:

            loaded_file = np.loadtxt(script_path +
                                     '/../../Individual/' + cavity_name +
                                     '/Resonators/single_resonator_h' +
                                     str(harmonic) +
                                     '.txt')

            R_S[indexloop_c, indexloop_h] = loaded_file[1]
            f_r[indexloop_c, indexloop_h] = loaded_file[0]
            Q[indexloop_c, indexloop_h] = loaded_file[2]

            indexloop_c += 1

        plt.figure('Impedance single Rs')
        plt.plot(f_r[:, indexloop_h]/1e6, R_S[:, indexloop_h]/1e3, '.',
                 label='h=%d' % (harmonic))

        plt.figure('Impedance single Q')
        plt.plot(f_r[:, indexloop_h]/1e6, Q[:, indexloop_h], '.',
                 label='h=%d' % (harmonic))

        indexloop_h += 1

    plt.figure('Impedance single Rs')
    plt.legend(loc='best')
    plt.xlabel('Resonant frequency $f_r$ [MHz]')
    plt.ylabel('Shunt impedance $R_s$ [$\\mathrm{k\\Omega}$]')
    plt.errorbar(np.mean(f_r, axis=0)/1e6,
                 np.mean(R_S, axis=0)/1e3,
                 xerr=np.std(f_r, axis=0)/1e6,
                 yerr=np.std(R_S, axis=0)/1e3,
                 fmt='bo--')
    plt.title('Average $R_s = %.2f~\\mathrm{k\\Omega}$' % (np.mean(R_S)/1e3))
    plt.tight_layout()
    plt.savefig(script_path+'/Rs_single_cavities_summed.png')

    plt.figure('Impedance single Q')
    plt.legend(loc='best')
    plt.xlabel('Resonant frequency $f_r$ [MHz]')
    plt.ylabel('Q factor $Q$')
    plt.errorbar(np.mean(f_r, axis=0)/1e6,
                 np.mean(Q, axis=0),
                 xerr=np.std(f_r, axis=0)/1e6,
                 yerr=np.std(Q, axis=0),
                 fmt='bo--')
    plt.tight_layout()
    plt.savefig(script_path+'/Q_single_cavities_summed.png')

