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

n_resonators = 2
harmonic_list = [8, 16, 21]


def load_model(harmonic):
    '''
    Method to load the parametric model
    Requires the rf frequency as an argument
    '''

    loaded_impedance = np.loadtxt(
        script_path+'/multi_resonator_and_ImZ_over_f_param_model.txt')

    f_r_interp = []
    R_S_interp = []
    Q_interp = []
    ImZ_over_f_interp = []

    for indexRes in range(n_resonators):

        f_r_interp.append(
            np.interp(harmonic,
                      loaded_impedance[indexRes::n_resonators, 0],
                      loaded_impedance[indexRes::n_resonators, 1]))

        R_S_interp.append(
            np.interp(
                harmonic,
                loaded_impedance[indexRes::n_resonators, 0],
                loaded_impedance[indexRes::n_resonators, 2]))

        Q_interp.append(
            np.interp(
                harmonic,
                loaded_impedance[indexRes::n_resonators, 0],
                loaded_impedance[indexRes::n_resonators, 3]))

        ImZ_over_f_interp.append(
            np.interp(
                harmonic,
                loaded_impedance[indexRes::n_resonators, 0],
                loaded_impedance[indexRes::n_resonators, 4]))

    f_r_interp = np.asarray(f_r_interp)
    R_S_interp = np.asarray(R_S_interp)
    Q_interp = np.asarray(Q_interp)
    ImZ_over_f_interp = np.sum(ImZ_over_f_interp)

    return f_r_interp, R_S_interp, Q_interp, ImZ_over_f_interp


if __name__ == '__main__':

    # Loading resonator fits and fitting evolution of parameters as a function
    # of the cavity harmonic

    R_S = np.zeros((len(harmonic_list), n_resonators))
    f_r = np.zeros((len(harmonic_list), n_resonators))
    Q = np.zeros((len(harmonic_list), n_resonators))
    ImZ_over_f = np.zeros((len(harmonic_list), n_resonators))

    indexloop = 0
    for harmonic in harmonic_list:

        loaded_file = np.loadtxt(script_path +
                                 '/../Resonators/multi_resonator_and_ImZ_over_f_h' +
                                 str(harmonic) +
                                 '.txt')

        R_S[indexloop] = loaded_file[:, 1]
        f_r[indexloop] = loaded_file[:, 0]
        Q[indexloop] = loaded_file[:, 2]
        ImZ_over_f[indexloop] = loaded_file[:, 3]

        indexloop += 1

    harmonic_save = np.array(np.sort(harmonic_list+harmonic_list), ndmin=2)
    R_S_save = np.array(R_S.reshape(n_resonators*len(harmonic_list)), ndmin=2)
    fr_save = np.array(f_r.reshape(n_resonators*len(harmonic_list)), ndmin=2)
    Q_save = np.array(Q.reshape(n_resonators*len(harmonic_list)), ndmin=2)
    ImZ_over_f_save = np.array(ImZ_over_f.reshape(n_resonators*len(harmonic_list)), ndmin=2)

    saved_matrix = np.hstack((harmonic_save.T, fr_save.T, R_S_save.T, Q_save.T, ImZ_over_f_save.T))
    np.savetxt(script_path+'/multi_resonator_and_ImZ_over_f_param_model.txt',
               saved_matrix, header='Impedance for 1xC10 cavities with closed gap relay based on C10-11 test' +
               '\nParametric model based on multi resonators fit with extra ImZ/f' +
               '\nAll contributions for a giver harmonic should be summed' +
               '\nAuthor: A. Lasheen\nh\t\t\tf_r [Hz]\t\t\tR_s [Ohm]\t\t\tQ\t\t\tImZ/f [Hz-1]')

    plt.figure('Impedance all')
    plt.clf()
    for indexRes in range(n_resonators):
        plt.plot(harmonic_list, R_S[:, indexRes]/1e3, 'bo--')
    plt.xlabel('Harmonic $h$')
    plt.ylabel('Shunt impedance $R_s$ [$\\mathrm{k\\Omega}$]')
    plt.grid('on')
    plt.tight_layout()
    plt.savefig(script_path+'/all_cavities_shunt.png')

    plt.figure('Freq all')
    plt.clf()
    for indexRes in range(n_resonators):
        plt.plot(harmonic_list, f_r[:, indexRes]/1e6, 'bo--')
    plt.xlabel('Harmonic $h$')
    plt.ylabel('Resonant frequency $f_r$ [MHz]')
    plt.grid('on')
    plt.tight_layout()
    plt.savefig(script_path+'/all_cavities_fr.png')

    plt.figure('Q all')
    plt.clf()
    for indexRes in range(n_resonators):
        plt.plot(harmonic_list, Q[:, indexRes], 'bo--')
    plt.xlabel('Harmonic $h$')
    plt.ylabel('Quality factor Q')
    plt.grid('on')
    plt.tight_layout()
    plt.savefig(script_path+'/all_cavities_Q.png')

    plt.figure('ImZ_over_f all')
    plt.clf()
    plt.plot(harmonic_list, np.sum(ImZ_over_f, axis=1), 'bo--')
    plt.xlabel('Harmonic $h$')
    plt.ylabel('ImZ/f [Hz-1]')
    plt.grid('on')
    plt.tight_layout()
    plt.savefig(script_path+'/all_cavities_ImZ_over_f.png')
