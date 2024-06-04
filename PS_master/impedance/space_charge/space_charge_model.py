'''
Created on 14 feb. 2018

@author: alasheen
'''

# General imports
import os
import numpy as np
from scipy.constants import c, epsilon_0

# Toolbox import
script_path = os.path.dirname(os.path.realpath(__file__))

# Localization of BI equipment
BI_devices_loc = {'BWSH54': 334.372,
                  'BWSV64': 397.204,
                  'BWSH65': 403.207,
                  'BWSH68': 422.616,
                  'BWSV85': 528.871,
                  'BGI82': 510.833}


def load_model(emittance_x_norm, emittance_y_norm, mass, momentum,
               momentum_spread_program, method='rectangle',
               aperture_X=None, aperture_Y=None,
               optics_S=None, beta_X=None,
               beta_Y=None, disp_X=None,
               reinterpolate=False,
               return_diagnose=False,
               BI_device_x='BWSH65',
               BI_device_y='BWSV85'):
    '''
    Method to compute space charge impedance over the momentum array

    The emittances should be passed in [m rad]
    The mass should be passed in [eV/c**2]
    The momentum should be passed in [eV/c]
    The momentum spread is dp/p [unitless]

    The apertures functions can be given directly by passing a
    2D matrix, the first dimension ([0,:]) being the longitudinal position
    and the second ([1,:]) being the actual function.
    The aperture correspond to the full aperture in [m]

    The aperture is interpolated linearly on the longitudinal position provided
    for the optics: optics_S.
    The optics are given as a 1d array, for which the longitudinal
    position corresponds to optics_S (all optics functions should
    share the same). 
    The optics functions are in [m]

    The dispersion in the vertical plane is neglected and set to 0

    '''

    if isinstance(momentum, (int, float)):
        momentum = [momentum]
    if isinstance(momentum_spread_program, (int, float)):
        momentum_spread_program = [momentum_spread_program]

    momentum = np.asarray(momentum)
    momentum_spread_program = np.asarray(momentum_spread_program)

    if len(momentum_spread_program) == 1:
        momentum_spread_program *= np.ones(len(momentum))

    final_ImZ_over_n = np.zeros(len(momentum))

    if return_diagnose:
        geometrical_factor_average = np.zeros(len(momentum))
        sigma_x_diag = np.zeros(len(momentum))
        sigma_y_diag = np.zeros(len(momentum))
        beta_x_diag = np.zeros(len(momentum))
        disp_x_diag = np.zeros(len(momentum))
        beta_y_diag = np.zeros(len(momentum))

    beta_rel_program = momentum / (np.sqrt(momentum**2. + mass**2.))
    gamma_rel_program = np.sqrt(momentum**2. + mass**2.) / mass

    if aperture_X is None:
        loaded_aperture_X = np.loadtxt(script_path + '/model/aperture_X.txt')
        sX = loaded_aperture_X[:, 0]
        aperExtX = loaded_aperture_X[:, 1]
        aperIntX = loaded_aperture_X[:, 2]
        aperture_X = aperExtX - aperIntX
        aperture_X = aperture_X[sX >= 0]
        sX = sX[sX >= 0]
    else:
        sX = aperture_X[0]
        aperture_X = aperture_X[1]

    if aperture_Y is None:
        loaded_aperture_Y = np.loadtxt(script_path + '/model/aperture_Y.txt')
        sY = loaded_aperture_Y[:, 0]
        aperExtY = loaded_aperture_Y[:, 1]
        aperIntY = loaded_aperture_Y[:, 2]
        aperture_Y = aperExtY - aperIntY
        aperture_Y = aperture_Y[sY >= 0]
        sY = sY[sY >= 0]
    else:
        sY = aperture_Y[0]
        aperture_Y = aperture_Y[1]

    if (beta_X is None) or (beta_Y is None) or (disp_X is None) \
            or (optics_S is None):
        loaded_optics = np.loadtxt(script_path + '/model/optics.txt')

    if optics_S is None:
        optics_S = loaded_optics[:, 0]
    else:
        optics_S = optics_S

    if beta_X is None:
        beta_X = loaded_optics[:, 2]

    if beta_Y is None:
        beta_Y = loaded_optics[:, 3]

    if disp_X is None:
        disp_X = loaded_optics[:, 4]

    disp_Y = 0.

    if reinterpolate:
        semiSizeX = np.interp(optics_S, sX, aperture_X / 2)
        semiSizeY = np.interp(optics_S, sY, aperture_Y / 2)
    else:
        semiSizeX = aperture_X / 2
        semiSizeY = aperture_Y / 2

    for index_momentum in range(len(momentum)):

        beta_rel = beta_rel_program[index_momentum]
        gamma_rel = gamma_rel_program[index_momentum]
        momentum_spread = momentum_spread_program[index_momentum]

        sigma_x = np.sqrt(emittance_x_norm / (beta_rel * gamma_rel) * beta_X +
                          (disp_X * momentum_spread)**2.)
        sigma_y = np.sqrt(emittance_y_norm / (beta_rel * gamma_rel) * beta_Y +
                          (disp_Y * momentum_spread)**2.)

        # Calculating the geometrical factor
        sigma_XY = (sigma_x + sigma_y) / 2

        # Rectangular
        if method == 'rectangle':
            geoFac = 1 / 4 + np.log(
                2 * 2 * semiSizeY / (np.pi * np.sqrt(2) * sigma_XY) *
                np.tanh(np.pi * 2 * semiSizeX / (2 * 2 * semiSizeY)))

        # Circular
        elif method == 'circle':
            geoFac = 1 / 4 + np.log(
                (semiSizeX + semiSizeY) / 2 / (np.sqrt(2) * sigma_XY))

        local_ImZ_over_n = - 1 / \
            (epsilon_0 * c * beta_rel * gamma_rel**2.) * geoFac

        final_ImZ_over_n[index_momentum] = np.trapz(
            local_ImZ_over_n, x=optics_S) / (optics_S[-1] - optics_S[0])

        if return_diagnose:
            geometrical_factor_average[index_momentum] = np.trapz(
                geoFac, x=optics_S) / (optics_S[-1] - optics_S[0])
            if BI_device_x in list(BI_devices_loc.keys()):
                sigma_x_diag[index_momentum] = np.interp(BI_devices_loc[BI_device_x],
                                                         optics_S, sigma_x)
                beta_x_diag[index_momentum] = np.interp(BI_devices_loc[BI_device_x],
                                                        optics_S, beta_X)
                disp_x_diag[index_momentum] = np.interp(BI_devices_loc[BI_device_x],
                                                        optics_S, disp_X)
            else:
                sigma_x_diag[index_momentum] = np.trapz(
                    sigma_x, x=optics_S) / (optics_S[-1] - optics_S[0])
                beta_x_diag[index_momentum] = np.trapz(
                    beta_X, x=optics_S) / (optics_S[-1] - optics_S[0])
                disp_x_diag[index_momentum] = np.trapz(
                    disp_X, x=optics_S) / (optics_S[-1] - optics_S[0])
            if BI_device_y in list(BI_devices_loc.keys()):
                sigma_y_diag[index_momentum] = np.interp(BI_devices_loc[BI_device_y],
                                                         optics_S, sigma_y)
                beta_y_diag[index_momentum] = np.interp(BI_devices_loc[BI_device_y],
                                                        optics_S, beta_Y)
            else:
                sigma_y_diag[index_momentum] = np.trapz(
                    sigma_y, x=optics_S) / (optics_S[-1] - optics_S[0])
                beta_y_diag[index_momentum] = np.trapz(
                    beta_Y, x=optics_S) / (optics_S[-1] - optics_S[0])

    if return_diagnose:
        return (final_ImZ_over_n, geometrical_factor_average,
                sigma_x_diag, sigma_y_diag, beta_x_diag, beta_y_diag,
                disp_x_diag)
    else:
        return final_ImZ_over_n


if __name__ == '__main__':

    # Check the ipynb on how the model was generated
    # Gathering all the required info from the aperture info in
    # '/eos/project/l/liu/' and the optics information from
    # https://gitlab.cern.ch/injector-optics/PS
    # The ipynb should be ran using python 2

    pass
