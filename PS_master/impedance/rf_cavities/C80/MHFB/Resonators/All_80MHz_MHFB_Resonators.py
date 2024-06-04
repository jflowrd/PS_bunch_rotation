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
C80_08_loaded_filter_data = np.load(
    script_path+'/C80-08_LUT_gains_and_delays.npz')
C80_88_loaded_filter_data = np.load(
    script_path+'/C80-88_LUT_gains_and_delays.npz')
C80_89_loaded_filter_data = np.load(
    script_path+'/C80-89_LUT_gains_and_delays.npz')


def load_model(filename_80MHz, revolution_frequency,
               impedance_reduction_target, bandwidth_target,
               main_harmonic_FB=None, frequency_array=None,
               manual_params_C80_08=None, manual_params_C80_88=None,
               manual_params_C80_89=None,
               RshFactor=1., QFactor=1.):
    '''
    - Method to load the parametric model
    Requires the main harmonic, the revolution frequency, the impedance
    reduction target on sidebands +-1 f_rev of the central harmonic.
    The feedback can be applied or removed using central_harmonic_FB=True/False
    '''

    if frequency_array is None:
        frequency_array = np.linspace(0, 500e6, int(1e6))

    final_impedance = 0

    harmonic_list = np.array(
        [162, 172, 163, 171, 164, 170, 165, 169, 166, 167]+[168]+[168]+[168])

    # Check used cavities
    C80_08_enabled = False
    C80_88_enabled = False
    C80_89_enabled = False
    if isinstance(filename_80MHz, str):
        if 'C08' in filename_80MHz.split('/'):
            C80_08_enabled = True
        if 'C88' in filename_80MHz.split('/'):
            C80_88_enabled = True
        if 'C89' in filename_80MHz.split('/'):
            C80_89_enabled = True
    elif isinstance(filename_80MHz, list) or \
            isinstance(filename_80MHz, np.ndarray):
        for filename in filename_80MHz:
            if 'C08' in filename.split('/'):
                C80_08_enabled = True
            if 'C88' in filename.split('/'):
                C80_88_enabled = True
            if 'C89' in filename.split('/'):
                C80_89_enabled = True

    # Check if central h should be damped
    C80_08_central_h_FB = False
    C80_88_central_h_FB = False
    C80_89_central_h_FB = False
    if isinstance(main_harmonic_FB, str):
        if 'C08' in main_harmonic_FB.split('/'):
            C80_08_central_h_FB = True
        if 'C88' in main_harmonic_FB.split('/'):
            C80_88_central_h_FB = True
        if 'C89' in main_harmonic_FB.split('/'):
            C80_89_central_h_FB = True
    elif isinstance(main_harmonic_FB, list) or \
            isinstance(main_harmonic_FB, np.ndarray):
        for main_harmonic_FB_elem in main_harmonic_FB:
            if 'C08' in main_harmonic_FB_elem.split('/'):
                C80_08_central_h_FB = True
            if 'C88' in main_harmonic_FB_elem.split('/'):
                C80_88_central_h_FB = True
            if 'C89' in main_harmonic_FB_elem.split('/'):
                C80_89_central_h_FB = True

    # C80-08
    if C80_08_enabled:
        loaded_C80_08 = np.loadtxt(
            script_path+'/../../Individual/C08/Resonators/single_resonator_b72.txt',
            skiprows=4)
        C80_08_Resonator = Resonators(loaded_C80_08[1]*RshFactor,
                                      loaded_C80_08[0],
                                      loaded_C80_08[2]*QFactor)

        C80_08_Resonator.imped_calc(frequency_array)

        if impedance_reduction_target < 0:
            if manual_params_C80_08 is None:

                C80_08_fittedParameters = np.zeros(
                    C80_08_loaded_filter_data[
                        'all_filter_parameters_save'].shape[-1])
                for index_param in range(
                        C80_08_loaded_filter_data[
                            'all_filter_parameters_save'].shape[-1]):

                    interp_func = RegularGridInterpolator(
                         (C80_08_loaded_filter_data[
                             'revolution_frequency_range'],
                          -C80_08_loaded_filter_data[
                              'impedance_reduction_target_range'],
                          C80_08_loaded_filter_data[
                              'bandwidth_target_range']),
                         C80_08_loaded_filter_data[
                             'all_filter_parameters_save'][:, :, :, index_param])

                    C80_08_fittedParameters[index_param] = interp_func(
                        np.array([revolution_frequency,
                                  -impedance_reduction_target,
                                  bandwidth_target]))

            else:
                C80_08_fittedParameters = list(manual_params_C80_08)

            if C80_08_central_h_FB:
                C80_08_fittedParameters[-2] = 0
                C80_08_fittedParameters[-5] = 0

            C80_08_final_filter = final_filter_calc(
                    frequency_array, revolution_frequency,
                    harmonic_list, *C80_08_fittedParameters)

            C80_08_final_impedance = final_impedance_calc(
                    C80_08_Resonator.impedance, C80_08_final_filter)

            final_impedance += C80_08_final_impedance

        else:

            final_impedance += C80_08_Resonator.impedance
            C80_08_fittedParameters = 0

    else:

        C80_08_fittedParameters = 0

    # C80-88
    if C80_88_enabled:
        loaded_C80_88 = np.loadtxt(
            script_path+'/../../Individual/C88/Resonators/single_resonator_b72.txt',
            skiprows=4)
        C80_88_Resonator = Resonators(loaded_C80_88[1]*RshFactor,
                                      loaded_C80_88[0],
                                      loaded_C80_88[2]*QFactor)

        C80_88_Resonator.imped_calc(frequency_array)

        if impedance_reduction_target < 0:
            if manual_params_C80_88 is None:

                C80_88_fittedParameters = np.zeros(
                    C80_88_loaded_filter_data[
                        'all_filter_parameters_save'].shape[-1])
                for index_param in range(
                    C80_88_loaded_filter_data[
                        'all_filter_parameters_save'].shape[-1]):

                    interp_func = RegularGridInterpolator(
                        (C80_88_loaded_filter_data[
                            'revolution_frequency_range'],
                         -C80_88_loaded_filter_data[
                             'impedance_reduction_target_range'],
                         C80_88_loaded_filter_data[
                             'bandwidth_target_range']),
                        C80_88_loaded_filter_data[
                            'all_filter_parameters_save'][:, :, :, index_param])

                    C80_88_fittedParameters[index_param] = interp_func(
                        np.array([revolution_frequency,
                                  -impedance_reduction_target,
                                  bandwidth_target]))

            else:
                C80_88_fittedParameters = list(manual_params_C80_88)

            if C80_88_central_h_FB:
                C80_88_fittedParameters[-2] = 0
                C80_88_fittedParameters[-5] = 0

            C80_88_final_filter = final_filter_calc(
                    frequency_array, revolution_frequency,
                    harmonic_list, *C80_88_fittedParameters)

            C80_88_final_impedance = final_impedance_calc(
                    C80_88_Resonator.impedance, C80_88_final_filter)

            final_impedance += C80_88_final_impedance

        else:

            final_impedance += C80_88_Resonator.impedance
            C80_88_fittedParameters = 0

    else:

        C80_88_fittedParameters = 0

    # C80-89
    if C80_89_enabled:
        loaded_C80_89 = np.loadtxt(
            script_path+'/../../Individual/C89/Resonators/single_resonator_b72.txt',
            skiprows=4)
        C80_89_Resonator = Resonators(loaded_C80_89[1]*RshFactor,
                                      loaded_C80_89[0],
                                      loaded_C80_89[2]*QFactor)

        C80_89_Resonator.imped_calc(frequency_array)

        if impedance_reduction_target < 0:
            if manual_params_C80_89 is None:

                C80_89_fittedParameters = np.zeros(
                    C80_89_loaded_filter_data[
                        'all_filter_parameters_save'].shape[-1])
                for index_param in range(
                    C80_89_loaded_filter_data[
                        'all_filter_parameters_save'].shape[-1]):

                    interp_func = RegularGridInterpolator(
                        (C80_89_loaded_filter_data[
                            'revolution_frequency_range'],
                         -C80_89_loaded_filter_data[
                             'impedance_reduction_target_range'],
                         C80_89_loaded_filter_data[
                             'bandwidth_target_range']),
                        C80_89_loaded_filter_data[
                            'all_filter_parameters_save'][:, :, :, index_param])

                    C80_89_fittedParameters[index_param] = interp_func(
                        np.array([revolution_frequency,
                                  -impedance_reduction_target,
                                  bandwidth_target]))

            else:
                C80_89_fittedParameters = list(manual_params_C80_89)

            if C80_89_central_h_FB:
                C80_89_fittedParameters[-2] = 0
                C80_89_fittedParameters[-5] = 0

            C80_89_final_filter = final_filter_calc(
                    frequency_array, revolution_frequency,
                    harmonic_list, *C80_89_fittedParameters)

            C80_89_final_impedance = final_impedance_calc(
                    C80_89_Resonator.impedance, C80_89_final_filter)

            final_impedance += C80_89_final_impedance

        else:

            final_impedance += C80_89_Resonator.impedance
            C80_89_fittedParameters = 0

    else:

        C80_89_fittedParameters = 0

    return final_impedance, frequency_array, (C80_08_fittedParameters,
                                              C80_88_fittedParameters,
                                              C80_89_fittedParameters)


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
     revolution_frequency, C80_impedance, delay_ratio, gain_ratio,
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

    final_impedance = final_impedance_calc(C80_impedance, final_filter)

    impedance_ratio_dB = impedance_to_dB(final_impedance,
                                         C80_impedance,
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
     harmonic_list, revolution_frequency, C80_impedance, delay_ratio,
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

    final_impedance = final_impedance_calc(C80_impedance, final_filter)

    impedance_ratio_dB = impedance_to_dB(final_impedance,
                                         C80_impedance,
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
     revolution_frequency, C80_impedance, delay_ratio, gain_ratio,
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

    final_impedance = final_impedance_calc(C80_impedance, final_filter)

    impedance_ratio_dB = impedance_to_dB(final_impedance,
                                         C80_impedance,
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
    save_result = False

    # Basic machine parameters
    revolution_frequency_range = [477e3]  # 3e8/(2*np.pi*100)
    impedance_reduction_target_range = [-20]
    bandwidth_target_range = [5e3]

#    revolution_frequency_range = np.arange(430e3, 481e3, 5e3)
#    impedance_reduction_target_range = -np.arange(10, 30, 2)
#    bandwidth_target_range = [5e3]

    main_harmonic = 168
    harmonic_list_init = np.array(
        [162, 172, 163, 171, 164, 170, 165, 169, 166, 167]+[168])
    cavity_name_selector = ['C80-08', 'C80-88', 'C80-89']

    all_filter_parameters_save = np.zeros(
        (len(cavity_name_selector),
         len(revolution_frequency_range),
         len(impedance_reduction_target_range),
         len(bandwidth_target_range),
         3*(len(harmonic_list_init)+2)))

    # Generating main impedances
    frequency_array = np.linspace(0, 160e6, 100000)

    loaded_C80_08 = np.loadtxt(
        script_path+'/../../Individual/C08/Resonators/single_resonator_b72.txt',
        skiprows=4)
    C80_08_Resonator = Resonators(loaded_C80_08[1],
                                  loaded_C80_08[0],
                                  loaded_C80_08[2])

    C80_08_Resonator.imped_calc(frequency_array)

    loaded_C80_88 = np.loadtxt(
        script_path+'/../../Individual/C88/Resonators/single_resonator_b72.txt',
        skiprows=4)
    C80_88_Resonator = Resonators(loaded_C80_88[1],
                                  loaded_C80_88[0],
                                  loaded_C80_88[2])

    C80_88_Resonator.imped_calc(frequency_array)

    loaded_C80_89 = np.loadtxt(
        script_path+'/../../Individual/C89/Resonators/single_resonator_b72.txt',
        skiprows=4)
    C80_89_Resonator = Resonators(loaded_C80_89[1],
                                  loaded_C80_89[0],
                                  loaded_C80_89[2])

    C80_89_Resonator.imped_calc(frequency_array)

    impedance_selector = [C80_08_Resonator.impedance,
                          C80_88_Resonator.impedance,
                          C80_89_Resonator.impedance]

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
                    fixed_gain = 5e-4 / gain_ratio
                    fixed_Q = 5000 / Q_ratio

                    startingParams_delay = []
                    for index_h in range(len(harmonic_list)-1):
                        startingParams_delay += [
                            (-index_h*13/10.+6)*1e-9 / delay_ratio]
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
                                    options={'disp': False,
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
                                options={'disp': False,
                                         'maxiter': 1000,
                                         'maxfev': 10000},
                                method='Nelder-Mead')['x']

                    # Central harmonic compensation
                    fixed_filter_params = []
                    for index_h in range(len(harmonic_list)):
                        fixed_filter_params += [
                            1/revolution_frequency/delay_ratio +
                            fittedParameters_delay[index_h]]
                        fixed_filter_params += [
                            fittedParameters_gain_delay[(index_h*2) % 2]]
                        fixed_filter_params += [
                            np.abs(fittedParameters_gain_delay[
                                (index_h*2+1) % 2])]

                    startingParams_central_compensation = [
                            1e-10/delay_ratio,
                            0.5e-4/gain_ratio,
                            10/1] # Removed Q_ratio

                    fittedParameters_central_compensation = minimize(
                                residue_function_central_flat,
                                startingParams_central_compensation,
                                args=(main_harmonic,
                                      harmonic_list,
                                      revolution_frequency,
                                      impedance_selector[index_cavity],
                                      delay_ratio,
                                      gain_ratio,
                                      Q_ratio,
                                      fixed_filter_params),
                                options={'disp': False,
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

                    fittedParameters += [fittedParameters[-3]]
                    fittedParameters += [-fittedParameters[-3]]
                    fittedParameters += [fittedParameters[-3]]

                    fittedParameters += [
                            1/revolution_frequency +
                            fittedParameters_central_compensation[0]*delay_ratio] # 
                    fittedParameters += [
                            fittedParameters_central_compensation[1]*gain_ratio]
                    fittedParameters += [
                        np.abs(fittedParameters_central_compensation[2])]

                    harmonic_list = np.array(
                        [162, 172, 163, 171, 164, 170, 165, 169, 166, 167]+[168]+[168]+[168])

                    print(fittedParameters)
                    all_filter_parameters_save[index_cavity,
                                               index_freq,
                                               index_red,
                                               index_bw,
                                               :] = np.array(fittedParameters)

                    if save_result:
                        np.savez(script_path+'/%s_LUT_gains_and_delays.npz' % (
                            cavity_name_selector[index_cavity]),
                            revolution_frequency_range=revolution_frequency_range,
                            impedance_reduction_target_range=impedance_reduction_target_range,
                            bandwidth_target_range=bandwidth_target_range,
                            harmonic_list=harmonic_list,
                            all_filter_parameters_save=all_filter_parameters_save[index_cavity, :, :, :])

                    if np.isclose(revolution_frequency, 477e3) and \
                            np.isclose(bandwidth_target_range, 5e3) and \
                            np.isclose(impedance_reduction_target, -20):

                        for plot_iteration in range(2):

                            central_bool_plot = True
                            if plot_iteration == 1:
                                central_bool_plot = False
                                fittedParameters[-5] = 0
                                fittedParameters[-2] = 0

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
                                cavity_name_selector[index_cavity]))
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
                            plt.savefig(script_path+'/%s_impedance_with_MHFB_-20dB_477kHz_central_%s.png' % (
                                cavity_name_selector[index_cavity],
                                central_bool_plot))

                            min_final_impedance_dB = np.nanmin(final_impedance_dB[
                                    (frequency_array > (main_harmonic*revolution_frequency-5e6/2)) *
                                    (frequency_array < (main_harmonic*revolution_frequency+5e6/2))])
                            max_final_impedance_dB = np.nanmax(final_impedance_dB[
                                    (frequency_array > (main_harmonic*revolution_frequency-5e6/2)) *
                                    (frequency_array < (main_harmonic*revolution_frequency+5e6/2))])
                            plt.figure('Impedance dB %s' % (cavity_name_selector[index_cavity]))
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
                            plt.savefig(script_path+'/%s_impedance_dB_with_MHFB_-20dB_477kHz_central_%s.png' % (
                                cavity_name_selector[index_cavity],
                                central_bool_plot))

                            min_impedance_ratio_dB = np.nanmin(impedance_ratio_dB[
                                    (frequency_array > (main_harmonic*revolution_frequency-5e6/2)) *
                                    (frequency_array < (main_harmonic*revolution_frequency+5e6/2))])
                            max_impedance_ratio_dB = np.nanmax(impedance_ratio_dB[
                                    (frequency_array > (main_harmonic*revolution_frequency-5e6/2)) *
                                    (frequency_array < (main_harmonic*revolution_frequency+5e6/2))])
                            plt.figure('Impedance ratio %s' % (
                                cavity_name_selector[index_cavity]))
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
                            plt.savefig(script_path+'/%s_impedance_ratio_dB_with_MHFB_-20dB_477kHz_central_%s.png'%(
                                cavity_name_selector[index_cavity], central_bool_plot))

#                            plt.figure('Impedance phase')
#                            plt.clf()
#                            plt.plot(frequency_array/1e6, np.angle(C10_Resonator.impedance))
#                            plt.plot(frequency_array/1e6, np.angle(final_impedance))
#                            plt.xlabel('Frequency [MHz]')
#                            plt.ylabel('Phase')
#                            plt.tight_layout()
