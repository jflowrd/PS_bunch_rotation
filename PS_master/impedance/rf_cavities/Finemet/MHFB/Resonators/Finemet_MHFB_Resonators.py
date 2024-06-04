'''
Created on 24 juin 2019

@author: alasheen
'''

# General imports
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import RegularGridInterpolator

# Toolbox import
sys.path.append('../../../impedance_toolbox')
from impedance_toolbox.impedance_toolbox.miniblond.impedances.impedance_sources import Resonators

# Script path
script_path = os.path.dirname(os.path.realpath(__file__))

# Loading LUT
Finemet_2013_loaded_filter_data = np.load(
    script_path+'/2013_LUT_gains_and_delays.npz')
Finemet_2018_loaded_filter_data = np.load(
    script_path+'/2018_LUT_gains_and_delays.npz')


def load_model(filename_finemet, revolution_frequency,
               impedance_reduction_target, bandwidth_target,
               frequency_array=None,
               manual_params_finemet=None,
               RshFactor=1., QFactor=1.):
    '''
    - Method to load the parametric model
    Requires the main harmonic, the revolution frequency, the impedance
    reduction target on the central harmonic.
    The feedback can be applied or removed using central_harmonic_FB=True/False
    '''

    if frequency_array is None:
        frequency_array = np.linspace(0, 1000e6, int(1e6))

    final_impedance = 0

    harmonic_list = np.arange(1, 13)

    # Check used cavities
    if 'Resonators_2013' in filename_finemet.split('/'):
        filepath = '/../../Resonators_2013/multi_resonators_impedance.txt'
        Finemet_loaded_filter_data = Finemet_2013_loaded_filter_data
    elif 'Resonators_2018' in filename_finemet.split('/'):
        filepath = '/../../Resonators_2018/multi_resonators_impedance.txt'
        Finemet_loaded_filter_data = Finemet_2018_loaded_filter_data

    loaded_finemet = np.loadtxt(
        script_path+filepath,
        skiprows=4)
    Finemet_Resonator = Resonators(loaded_finemet[:, 1]*RshFactor,
                                   loaded_finemet[:, 0],
                                   loaded_finemet[:, 2]*QFactor)

    Finemet_Resonator.imped_calc(frequency_array)

    if impedance_reduction_target < 0:
        if manual_params_finemet is None:

            finemet_fittedParameters = np.zeros(
                Finemet_loaded_filter_data[
                    'all_filter_parameters_save'].shape[-1])
            for index_param in range(
                    Finemet_loaded_filter_data[
                        'all_filter_parameters_save'].shape[-1]):

                    finemet_fittedParameters[index_param] = np.interp(
                        revolution_frequency,
                        Finemet_loaded_filter_data[
                            'revolution_frequency_range'],
                        Finemet_loaded_filter_data[
                            'all_filter_parameters_save'][
                                :, 0, 0, index_param])

#                 interp_func = RegularGridInterpolator(
#                      (Finemet_loaded_filter_data[
#                          'revolution_frequency_range'],
#                       -Finemet_loaded_filter_data[
#                           'impedance_reduction_target_range'],
#                       Finemet_loaded_filter_data[
#                           'bandwidth_target_range']),
#                      Finemet_loaded_filter_data[
#                          'all_filter_parameters_save'][:, :, :, index_param])
# 
#                 finemet_fittedParameters[index_param] = interp_func(
#                     np.array([revolution_frequency,
#                               -impedance_reduction_target,
#                               bandwidth_target]))

        else:
            finemet_fittedParameters = list(manual_params_finemet)

        finemet_final_filter = final_filter_calc(
                frequency_array, revolution_frequency,
                harmonic_list, *finemet_fittedParameters)

        final_impedance = final_impedance_calc(
                Finemet_Resonator.impedance, finemet_final_filter)

    else:

        final_impedance = Finemet_Resonator.impedance
        finemet_fittedParameters = 0

    return final_impedance, frequency_array, finemet_fittedParameters


def final_filter_calc(frequency_array, revolution_frequency, harmonic_list,
                      *params):

    final_filter = 0
    for index_h in range(len(harmonic_list)):

        resonator_filter = Resonators(
            1,
            harmonic_list[index_h]*revolution_frequency, params[index_h*3+2])
        resonator_filter.imped_calc(frequency_array)
        resonator_filter = resonator_filter.impedance

        delay = params[index_h*3]
        gain = params[index_h*3+1]

        final_filter += gain*resonator_filter*np.exp(
            -1j*2*np.pi*frequency_array*delay)

    return final_filter


def final_impedance_calc(initial_impedance, final_filter):

    final_impedance = initial_impedance / (1+initial_impedance*final_filter)

    return final_impedance


def impedance_to_dB(impedance, reference=None, ratio=False):

    if reference is None:
        return 20*np.log10(np.abs(impedance)/np.nanmax(np.abs(impedance)))
    elif isinstance(reference, float) or (
            isinstance(reference, np.ndarray) and not ratio):
        return 20*np.log10(np.abs(impedance)/np.nanmax(np.abs(reference)))
    elif isinstance(reference, float) or (
            isinstance(reference, np.ndarray) and ratio):
        return 20*np.log10(np.abs(impedance)/np.abs(reference))
    else:
        return None


def residue_function_delay_single(params, *args):

    (harmonic_list, selected_index_h,
     revolution_frequency, cavity_impedance, delay_ratio, gain_ratio,
     Q_ratio, fixed_delay, fixed_gain, fixed_Q) = args

    filter_params = []
    for index_h in range(len(harmonic_list)):
        if index_h == selected_index_h:
            filter_params += [1/revolution_frequency+params*delay_ratio]
        else:
            filter_params += [1/revolution_frequency+fixed_delay[index_h]*delay_ratio]
        filter_params += [fixed_gain*gain_ratio]
        filter_params += [fixed_Q*Q_ratio]

    final_filter = final_filter_calc(frequency_array,
                                     revolution_frequency,
                                     harmonic_list,
                                     *filter_params)

    final_impedance = final_impedance_calc(cavity_impedance, final_filter)

    impedance_ratio_dB = impedance_to_dB(final_impedance,
                                         cavity_impedance,
                                         True)

    impedance_in_window_center = impedance_ratio_dB[
            (frequency_array >
             (harmonic_list[selected_index_h]*revolution_frequency -
              revolution_frequency/2)) *
            (frequency_array <
             (harmonic_list[selected_index_h]*revolution_frequency +
              revolution_frequency/2))]

    min_center = np.nanmin(impedance_in_window_center)

    residue_center = np.exp(min_center)

    residue = residue_center

    residue = np.sqrt(residue)

    return residue


def residue_function_gain_Q(params, *args):

    (impedance_reduction_target, bandwdith_target, main_harmonic,
     harmonic_list, revolution_frequency, cavity_impedance, delay_ratio,
     gain_ratio, Q_ratio, fixed_delay) = args

    filter_params = []

    for index_h in range(len(harmonic_list)):
        filter_params += [fixed_delay[index_h]*delay_ratio]
        filter_params += [params[(index_h*2) % 2]*gain_ratio]
        filter_params += [np.abs(params[(index_h*2+1) % 2]*Q_ratio)]

    final_filter = final_filter_calc(frequency_array,
                                     revolution_frequency,
                                     harmonic_list,
                                     *filter_params)

    final_impedance = final_impedance_calc(cavity_impedance, final_filter)

    impedance_ratio_dB = impedance_to_dB(final_impedance,
                                         cavity_impedance,
                                         True)

    impedance_in_window_center = impedance_ratio_dB[
            (frequency_array >
             (main_harmonic*revolution_frequency-revolution_frequency/4)) *
            (frequency_array <
             (main_harmonic*revolution_frequency+revolution_frequency/4))]

    min_center = np.nanmin(impedance_in_window_center)

    residue_center = (min_center-impedance_reduction_target)**2.

    impedance_ratio_dB -= min_center + 3

    frequency_array_in_window_left = frequency_array[
            (frequency_array >
             (main_harmonic*revolution_frequency-revolution_frequency/4)) *
            (frequency_array <
             (main_harmonic*revolution_frequency))]
    frequency_array_in_window_right = frequency_array[
            (frequency_array >
             (main_harmonic*revolution_frequency)) *
            (frequency_array <
             (main_harmonic*revolution_frequency+revolution_frequency/4))]

    impedance_in_window_left = impedance_ratio_dB[
            (frequency_array >
             (main_harmonic*revolution_frequency-revolution_frequency/4)) *
            (frequency_array <
             (main_harmonic*revolution_frequency))]
    impedance_in_window_right = impedance_ratio_dB[
            (frequency_array >
             (main_harmonic*revolution_frequency)) *
            (frequency_array <
             (main_harmonic*revolution_frequency+revolution_frequency/4))]

    left_freq_BW = np.interp(0, -impedance_in_window_left,
                             frequency_array_in_window_left)
    right_freq_BW = np.interp(0, impedance_in_window_right,
                              frequency_array_in_window_right)
    bandwidth_3dB = right_freq_BW-left_freq_BW

    residue_BW = ((bandwidth_3dB-bandwdith_target)/1e3)**2.

    residue = np.sqrt(residue_center + residue_BW)

    return residue


def residue_function_central_flat(params, *args):

    (main_harmonic, harmonic_list,
     revolution_frequency, cavity_impedance, delay_ratio, gain_ratio,
     Q_ratio, fixed_filter_params) = args

    filter_params = []

    for index_h in range(len(harmonic_list)):
        filter_params += [fixed_filter_params[index_h*3]*delay_ratio]
        filter_params += [fixed_filter_params[index_h*3+1]*gain_ratio]
        filter_params += [fixed_filter_params[index_h*3+2]*Q_ratio]

    filter_params += [filter_params[-3]]
    filter_params += [-filter_params[-3]]
    filter_params += [filter_params[-3]]

    filter_params += [1/revolution_frequency+params[0]*delay_ratio] #
    filter_params += [params[1]*gain_ratio]
    filter_params += [np.abs(params[2])]

    new_harmonic_list = np.append(harmonic_list, np.array([main_harmonic,
                                                           main_harmonic]))

    final_filter = final_filter_calc(frequency_array,
                                     revolution_frequency,
                                     new_harmonic_list,
                                     *filter_params)

    final_impedance = final_impedance_calc(cavity_impedance, final_filter)

    impedance_ratio_dB = impedance_to_dB(final_impedance,
                                         cavity_impedance,
                                         True)

    impedance_in_window_center = impedance_ratio_dB[
            (frequency_array >
             (main_harmonic*revolution_frequency-revolution_frequency/2)) *
            (frequency_array <
             (main_harmonic*revolution_frequency+revolution_frequency/2))]

    residue = np.sqrt(np.sum(impedance_in_window_center**2.))

    return residue


if __name__ == '__main__':

    # Save result
    save_result = True

    # Basic machine parameters
    revolution_frequency_range = [458e3, 459e3]  # 3e8/(2*np.pi*100)
    impedance_reduction_target_range = [-15]
    bandwidth_target_range = [8e3]

#    revolution_frequency_range = np.arange(430e3, 481e3, 5e3)
#    impedance_reduction_target_range = -np.arange(10, 30, 2)
#    bandwidth_target_range = [5e3]

    main_harmonic = 7
    harmonic_list_init = np.arange(1, 13)
    case_selector = ['2013', '2018']

    all_filter_parameters_save = np.zeros(
        (len(case_selector),
         len(revolution_frequency_range),
         len(impedance_reduction_target_range),
         len(bandwidth_target_range),
         3*(len(harmonic_list_init))))

    # Generating main impedances
    frequency_array = np.linspace(0, 20e6, 100000)

    loaded_Finemet_2013 = np.loadtxt(
        script_path +
        '/../../Resonators_2013/multi_resonators_impedance.txt',
        skiprows=4)
    Finemet_2013_Resonator = Resonators(loaded_Finemet_2013[:, 1],
                                        loaded_Finemet_2013[:, 0],
                                        loaded_Finemet_2013[:, 2])
    Finemet_2013_Resonator.imped_calc(frequency_array)

    loaded_Finemet_2018 = np.loadtxt(
        script_path +
        '/../../Resonators_2018/multi_resonators_impedance.txt',
        skiprows=4)
    Finemet_2018_Resonator = Resonators(loaded_Finemet_2018[:, 1],
                                        loaded_Finemet_2018[:, 0],
                                        loaded_Finemet_2018[:, 2])
    Finemet_2018_Resonator.imped_calc(frequency_array)

    impedance_selector = [Finemet_2013_Resonator.impedance,
                          Finemet_2018_Resonator.impedance]

    for index_freq in range(len(revolution_frequency_range)):
        for index_red in range(len(impedance_reduction_target_range)):
            for index_bw in range(len(bandwidth_target_range)):

                print('%d/%d-%d/%d-%d/%d' % (
                    index_freq+1,
                    len(revolution_frequency_range),
                    index_red+1,
                    len(impedance_reduction_target_range),
                    index_bw+1,
                    len(bandwidth_target_range)))

                revolution_frequency = revolution_frequency_range[index_freq]
                impedance_reduction_target = impedance_reduction_target_range[index_red]
                bandwidth_target = bandwidth_target_range[index_bw]

                all_revolution_freq = harmonic_list_init*revolution_frequency
                main_harmonic_freq = main_harmonic*revolution_frequency

                # Generating filter
                for index_cavity in range(len(impedance_selector)):

                    harmonic_list = np.array(harmonic_list_init)

                    # Fit ratio
                    delay_ratio = 1e-9
                    gain_ratio = 1e-4
                    Q_ratio = 1e3

                    # Adjusting delay
                    fixed_gain = 1e-4 / gain_ratio
                    fixed_Q = 1000 / Q_ratio

                    startingParams_delay = []
                    for index_h in range(len(harmonic_list)):
                        startingParams_delay += [0.1e-9 / delay_ratio]
                    fittedParameters_delay = np.array(startingParams_delay)

                    for reiteration in range(2):
                        for selected_index_h in range(len(harmonic_list)):
                            fittedParameters_delay[
                                selected_index_h] = minimize(
                                    residue_function_delay_single,
                                    fittedParameters_delay[selected_index_h],
                                    args=(harmonic_list,
                                          selected_index_h,
                                          revolution_frequency,
                                          impedance_selector[index_cavity],
                                          delay_ratio,
                                          gain_ratio,
                                          Q_ratio,
                                          fittedParameters_delay,
                                          fixed_gain,
                                          fixed_Q),
                                    options={'disp': True,
                                             'maxiter': 1000,
                                             'maxfev': 10000},
                                    method='Nelder-Mead')['x']

                    # Adjusting gain and Q
                    fixed_delay = 1/revolution_frequency/delay_ratio+fittedParameters_delay
                    startingParams_gain_Q = [fixed_gain,
                                             fixed_Q]

                    fittedParameters_gain_delay = minimize(
                                residue_function_gain_Q,
                                startingParams_gain_Q,
                                args=(impedance_reduction_target,
                                      bandwidth_target,
                                      main_harmonic,
                                      harmonic_list,
                                      revolution_frequency,
                                      impedance_selector[index_cavity],
                                      delay_ratio,
                                      gain_ratio,
                                      Q_ratio,
                                      fixed_delay),
                                options={'disp': True,
                                         'maxiter': 1000,
                                         'maxfev': 10000},
                                method='Nelder-Mead')['x']

                    # Final results
                    fittedParameters = []
                    for index_h in range(len(harmonic_list)):
                        fittedParameters += [
                            1/revolution_frequency +
                            fittedParameters_delay[index_h]*delay_ratio]
                        fittedParameters += [
                            fittedParameters_gain_delay[
                                (index_h*2) % 2]*gain_ratio]
                        fittedParameters += [
                            np.abs(fittedParameters_gain_delay[
                                (index_h*2+1) % 2]*Q_ratio)]

                    print(fittedParameters)
                    all_filter_parameters_save[index_cavity,
                                               index_freq,
                                               index_red,
                                               index_bw,
                                               :] = np.array(fittedParameters)

                    if save_result:
                        np.savez(script_path+'/%s_LUT_gains_and_delays.npz' % (
                            case_selector[index_cavity]),
                            revolution_frequency_range=revolution_frequency_range,
                            impedance_reduction_target_range=impedance_reduction_target_range,
                            bandwidth_target_range=bandwidth_target_range,
                            harmonic_list=harmonic_list,
                            all_filter_parameters_save=all_filter_parameters_save[index_cavity, :, :, :])

                    if np.isclose(revolution_frequency, 458e3) and \
                            np.isclose(bandwidth_target_range, 8e3) and \
                            np.isclose(impedance_reduction_target, -15):

                        # Final impedance
                        final_filter = final_filter_calc(
                            frequency_array, revolution_frequency,
                            harmonic_list, *fittedParameters)

                        final_impedance = final_impedance_calc(
                            impedance_selector[index_cavity], final_filter)

                        main_impedance_dB = impedance_to_dB(
                            impedance_selector[index_cavity])

                        final_impedance_dB = impedance_to_dB(
                            final_impedance,
                            np.abs(impedance_selector[index_cavity]))

                        impedance_ratio_dB = impedance_to_dB(
                            final_impedance,
                            np.abs(impedance_selector[index_cavity]),
                            True)

                        plt.figure('Filter')
                        plt.clf()
                        plt.plot(frequency_array, np.abs(final_filter))

                        plt.figure('Filter phase')
                        plt.clf()
                        plt.plot(frequency_array, np.angle(final_filter))

                        plt.figure('Impedance %s' % (
                            case_selector[index_cavity]))
                        plt.clf()
                        plt.plot(frequency_array/1e6,
                                 np.abs(impedance_selector[index_cavity])/1e3)
                        plt.plot(frequency_array/1e6,
                                 np.abs(final_impedance)/1e3)
                        plt.xlim(
                            (main_harmonic*revolution_frequency/1e6-6/2,
                             main_harmonic*revolution_frequency/1e6+6/2))
                        plt.hlines(0,
                                   main_harmonic*revolution_frequency/1e6-6/2,
                                   main_harmonic*revolution_frequency/1e6+6/2,
                                   linewidth=0.5, alpha=0.5)
                        plt.vlines(all_revolution_freq/1e6, 0,
                                   np.nanmax(np.abs(impedance_selector[index_cavity])/1e3),
                                   linewidth=0.5, alpha=0.5)
                        plt.xlabel('Frequency [MHz]')
                        plt.ylabel('Magnitude [$k\\Omega$]')
                        plt.tight_layout()
                        plt.savefig(script_path+'/%s_impedance_with_MHFB_-15dB_458kHz.png' % (
                            case_selector[index_cavity]))

                        min_final_impedance_dB = np.nanmin(final_impedance_dB[
                                (frequency_array > (main_harmonic*revolution_frequency-5e6/2)) *
                                (frequency_array < (main_harmonic*revolution_frequency+5e6/2))])
                        max_final_impedance_dB = np.nanmax(final_impedance_dB[
                                (frequency_array > (main_harmonic*revolution_frequency-5e6/2)) *
                                (frequency_array < (main_harmonic*revolution_frequency+5e6/2))])
                        plt.figure('Impedance dB %s' % (case_selector[index_cavity]))
                        plt.clf()
                        plt.plot(frequency_array/1e6, main_impedance_dB)
                        plt.plot(frequency_array/1e6, final_impedance_dB)
                        plt.xlim(
                            (main_harmonic*revolution_frequency/1e6-6/2,
                             main_harmonic*revolution_frequency/1e6+6/2))
                        plt.ylim((1.1*min_final_impedance_dB,
                                  max_final_impedance_dB+1))
                        plt.vlines(all_revolution_freq/1e6,
                                   min_final_impedance_dB,
                                   max_final_impedance_dB,
                                   linewidth=0.5,
                                   alpha=0.5)
                        plt.xlabel('Frequency [MHz]')
                        plt.ylabel('Magnitude [dB]')
                        plt.tight_layout()
                        plt.savefig(script_path+'/%s_impedance_dB_with_MHFB_-15dB_458kHz.png' % (
                            case_selector[index_cavity]))

                        min_impedance_ratio_dB = np.nanmin(impedance_ratio_dB[
                                (frequency_array > (main_harmonic*revolution_frequency-5e6/2)) *
                                (frequency_array < (main_harmonic*revolution_frequency+5e6/2))])
                        max_impedance_ratio_dB = np.nanmax(impedance_ratio_dB[
                                (frequency_array > (main_harmonic*revolution_frequency-5e6/2)) *
                                (frequency_array < (main_harmonic*revolution_frequency+5e6/2))])
                        plt.figure('Impedance ratio %s' % (
                            case_selector[index_cavity]))
                        plt.clf()
                        plt.plot(frequency_array/1e6, impedance_ratio_dB)
                        plt.xlim(
                            (main_harmonic*revolution_frequency/1e6-6/2,
                             main_harmonic*revolution_frequency/1e6+6/2))
                        plt.vlines(all_revolution_freq/1e6,
                                   min_impedance_ratio_dB,
                                   max_impedance_ratio_dB,
                                   linewidth=0.5, alpha=0.5)
                        plt.hlines(0,
                                   main_harmonic*revolution_frequency/1e6-6/2,
                                   main_harmonic*revolution_frequency/1e6+6/2,
                                   linewidth=0.5, alpha=0.5)
                        plt.xlabel('Frequency [MHz]')
                        plt.ylabel('Magnitude [dB]')
                        plt.tight_layout()
                        plt.savefig(script_path+'/%s_impedance_ratio_dB_with_MHFB_-15dB_458kHz.png'%(
                            case_selector[index_cavity]))

#                            plt.figure('Impedance phase')
#                            plt.clf()
#                            plt.plot(frequency_array/1e6, np.angle(C10_Resonator.impedance))
#                            plt.plot(frequency_array/1e6, np.angle(final_impedance))
#                            plt.xlabel('Frequency [MHz]')
#                            plt.ylabel('Phase')
#                            plt.tight_layout()
