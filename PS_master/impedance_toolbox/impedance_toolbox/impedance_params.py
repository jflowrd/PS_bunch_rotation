'''
Created on 15 May 2017

@author: alasheen
'''

import warnings
import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import interp1d, sproot
from scipy.constants import e
from .miniblond.impedances.impedance_sources import Resonators, InputTable

from .residual_types import log_real_imag, log_total, lin_real_imag, lin_total


import matplotlib.pyplot as plt


class ImpedanceParameters(object):
    '''
    Object that gathers all impedance related functions
    '''

    def __init__(self, directory, MachineParameters=None):
        '''
        Constructor
        '''

        self.directory = directory

        self.machineParams = MachineParameters

        self.fitResiduals = {'log_real_and_imag': log_real_imag,
                             'log_total': log_total,
                             'lin_real_and_imag': lin_real_imag,
                             'lin_total': lin_total
                             }

    def addResonators(self, R_S, frequency_R, Q, freqArray=None):

        if hasattr(self, 'freqArray'):

            pass

        elif hasattr(self.machineParams, 'freqArray'):

            self.freqArray = self.machineParams.freqArray

        elif freqArray is not None:

            self.freqArray = freqArray

        elif freqArray is None and not hasattr(self.machineParams, 'freqArray'):

            raise RuntimeError("The addResonators function requires to have " +
                               "a freqArray defined or a MachineParameters " +
                               "objects with a calculated beam spectrum")

        resonators = Resonators(R_S, frequency_R, Q)

        resonators.imped_calc(self.freqArray)

        if hasattr(self, 'impedance'):

            self.impedance += resonators.impedance

        else:

            self.impedance = resonators.impedance

    def addImpedanceTable(self, frequency, realPart, imagPart,
                          freqArrayInterp=None):

        if hasattr(self, 'freqArray'):

            pass

        elif hasattr(self.machineParams, 'freqArray'):

            self.freqArray = self.machineParams.freqArray

        elif freqArrayInterp is not None:

            self.freqArray = freqArrayInterp

        elif freqArrayInterp is None and not hasattr(self.machineParams,
                                                     'freqArray'):

            raise RuntimeError("The addImpedanceTable function requires to have " +
                               "a freqArray defined or a MachineParameters " +
                               "objects with a calculated beam spectrum")

        impedanceTable = InputTable(frequency, realPart, imagPart)

        impedanceTable.imped_calc(self.freqArray)

        if hasattr(self, 'impedance'):

            self.impedance += impedanceTable.impedance

        else:

            self.impedance = impedanceTable.impedance

    def importResonatorsList(self, ResonatorsList, freqArray=None):

        if hasattr(self, 'freqArray'):

            pass

        elif hasattr(self.machineParams, 'freqArray'):

            self.freqArray = self.machineParams.freqArray

        elif freqArray is not None:

            self.freqArray = freqArray

        elif freqArray is None and not hasattr(self.machineParams, 'freqArray'):

            raise RuntimeError("The importResonatorsList function requires to have " +
                               "a freqArray defined or a MachineParameters " +
                               "objects with a calculated beam spectrum")

        for resonators in ResonatorsList:

            resonators.imped_calc(self.freqArray)

            if hasattr(self, 'impedance'):

                self.impedance += resonators.impedance

            else:

                self.impedance = resonators.impedance

    def importImpedanceTableList(self, ImpedanceTableList, freqArray=None):

        if hasattr(self, 'freqArray'):

            pass

        elif hasattr(self.machineParams, 'freqArray'):

            self.freqArray = self.machineParams.freqArray

        elif freqArray is not None:

            self.freqArray = freqArray

        elif freqArray is None and not hasattr(self.machineParams, 'freqArray'):

            raise RuntimeError("The importImpedanceTableList function requires to have " +
                               "a freqArray defined or a MachineParameters " +
                               "objects with a calculated beam spectrum")

        for impedance_table in ImpedanceTableList:

            impedance_table.imped_calc(self.freqArray)

            if hasattr(self, 'impedance'):

                self.impedance += impedance_table.impedance

            else:

                self.impedance = impedance_table.impedance

    def importImZ_over_f_List(self, ImZ_over_f_List, freqArray=None):

        if hasattr(self, 'freqArray'):

            pass

        elif hasattr(self.machineParams, 'freqArray'):

            self.freqArray = self.machineParams.freqArray

        elif freqArray is not None:

            self.freqArray = freqArray

        elif freqArray is None and not hasattr(self.machineParams, 'freqArray'):

            raise RuntimeError("The importImpedanceTableList function requires to have " +
                               "a freqArray defined or a MachineParameters " +
                               "objects with a calculated beam spectrum")

        for ImZ_over_f in ImZ_over_f_List:

            impedance_table = np.zeros(len(self.freqArray), dtype=complex)
            impedance_table[:].imag = ImZ_over_f * self.freqArray

            if hasattr(self, 'impedance'):

                self.impedance += impedance_table

            else:

                self.impedance = impedance_table

    def inducedVoltageGeneration(self):

        if hasattr(self, 'beamSpectrum'):

            pass

        elif hasattr(self.machineParams, 'beamSpectrum'):

            self.beamSpectrum = self.machineParams.beamSpectrum
            self.timeArray = self.machineParams.timeArray

        else:

            raise RuntimeError('A beam spectrum needs to be calculated in ' +
                               'the MachineParameters and passed as an ' +
                               'argument to the ImpedanceParameters object')

        self.inducedVoltage = - e * (
            self.machineParams.nTotalBunches *
            self.machineParams.intensityPerBunch /
            np.sum(self.machineParams.beamProfileNormalized)) * \
            np.fft.irfft(self.impedance *
                         self.machineParams.beamSpectrumNormalized) * \
            (self.freqArray[1] - self.freqArray[0]) * \
            (len(self.beamSpectrum) - 1.) * len(self.timeArray)

        self.inducedVoltage = self.inducedVoltage

    def getRFLosses(self):
        '''
        Losses per turn
        '''

        len_array_time = np.min([len(self.machineParams.beamProfileNormalized),len(self.inducedVoltage)])
        self.energyLossTime = np.trapz(self.machineParams.beamProfileNormalized[:len_array_time]
                                       * self.machineParams.intensityPerBunch
                                       * self.machineParams.nTotalBunches *
                                       self.inducedVoltage[:len_array_time],
                                       dx=self.timeArray[1] - self.timeArray[0]) \
            / self.machineParams.nTurns

        self.energyLossFreq = -(e * self.machineParams.intensityPerBunch
                                * self.machineParams.nTotalBunches / self.machineParams.nTurns)**2. \
            * self.machineParams.generalParams.omega_rev[0] / np.pi \
            * np.sum(self.impedance.real * (np.abs(self.machineParams.beamSpectrumNormalized) /
                                            (np.abs(self.machineParams.beamSpectrumNormalized)[0]))**2.) \
            / e

        self.powerLossTime = self.energyLossTime * e / self.machineParams.nTurns \
            / self.machineParams.generalParams.t_rev[0]

        self.powerLossFreq = self.energyLossFreq * e / \
            self.machineParams.generalParams.t_rev[0]

    def scanFit(self, R_S, frequency_R, Q, frequencyWindow=None,
                RShuntBound=None, freqBound=None, QBound=None,
                RShuntScale=1, freqScale=1, QScale=1,
                method=None, fitResidue='log_total',
                disp=False, itt_max=5, nDivide=False,
                fit_target='R/Q'):

        if frequencyWindow is not None:
            self.freqIndexes = np.where(
                (self.freqArray > frequencyWindow[0]) * (self.freqArray < frequencyWindow[1]))[0]
        else:
            self.freqIndexes = np.arange(len(self.freqArray))

        impedanceTable = np.zeros([3, self.freqIndexes.shape[0]])
        impedanceTable[0] = self.freqArray[self.freqIndexes]
        impedanceTable[1] = self.impedance.real[self.freqIndexes]
        impedanceTable[2] = self.impedance.imag[self.freqIndexes]

        # Sort resonators by R_Shunt
        sort_index = np.argsort(frequency_R)

        R_S = np.array(R_S, dtype=float, ndmin=1)[
            list(reversed(sort_index.tolist()))]
        frequency_R = np.array(frequency_R, dtype=float, ndmin=1)[
            list(reversed(sort_index.tolist()))]
        Q = np.array(Q, dtype=float, ndmin=1)[
            list(reversed(sort_index.tolist()))]

        n_resonators = R_S.shape[0]

        self.resonatorForFit = Resonators(R_S * RShuntScale,
                                          frequency_R * freqScale,
                                          Q * QScale)

        self.resonatorForFit.imped_calc(impedanceTable[0])

        # Used for optional division by frequency
        if nDivide:
            divisor = impedanceTable[0]
        else:
            divisor = 1

        for j in range(itt_max):

            for i in range(n_resonators):

                if fit_target == 'R/Q':
                    fittedParameters = minimize(self.scanFitFunction, [R_S[i], Q[i]],
                                                args=(
                                                    impedanceTable, i, R_S, frequency_R, Q, 'R/Q', divisor),
                                                method='Powell')['x']

                    R_S[i] = fittedParameters[0]
                    Q[i] = fittedParameters[1]

                elif fit_target == 'Freq':

                    fittedParameters = minimize(self.scanFitFunction, [frequency_R[i]],
                                                args=(
                                                    impedanceTable, i, R_S, frequency_R, Q, 'Freq', divisor),
                                                method='Powell')['x']

                    frequency_R[i] = fittedParameters[0]

                elif fit_target == 'All':

                    fittedParameters = minimize(self.scanFitFunction, [frequency_R[i]],
                                                args=(
                                                    impedanceTable, i, R_S, frequency_R, Q, 'All', divisor),
                                                method='Powell')['x']

                    R_S[i] = fittedParameters[0]
                    frequency_R[i] = fittedParameters[1]
                    Q[i] = fittedParameters[2]

                else:

                    raise RuntimeError("Fit target not recognised")

        return [R_S, frequency_R, Q]

    def fitInitialGuess(self, level=None, plot=False):
        '''
        Generates an initial guess for startingParams in fitResonators
        '''

        absoluteImpedance = np.abs(self.impedance)

        [min_x_position,
         max_x_position], [min_values,
                           max_values] = minmax_location(
                               self.freqArray,
                               absoluteImpedance)

        if level is None:
            R_S_start = max_values
            frequency_R_start = max_x_position
        else:
            max_RS = np.max(max_values)
            R_S_start = max_values[max_values > level * max_RS]
            frequency_R_start = max_x_position[max_values > level * max_RS]

        Q_start = np.zeros(len(R_S_start))

        for indexRes in range(len(R_S_start)):

            freqRange = self.freqArray[(self.freqArray > frequency_R_start[indexRes] * 0.8) * (
                self.freqArray < frequency_R_start[indexRes] * 1.2)]
            impedanceRange = absoluteImpedance[(self.freqArray > frequency_R_start[indexRes] * 0.8) * (
                self.freqArray < frequency_R_start[indexRes] * 1.2)]
            deltaF = freqRange[impedanceRange > 1 / 10**(3. / 20) * np.max(
                impedanceRange)][-1] - freqRange[impedanceRange > 1 / 10**(3. / 20) * np.max(impedanceRange)][0]

            Q_start[indexRes] = frequency_R_start[indexRes] / deltaF

        if plot:
            plt.figure('fitInitialGuess')
            plt.clf()
            plt.plot(self.freqArray, absoluteImpedance)
            plt.plot(frequency_R_start, R_S_start, 'o')
            plt.show()

        return R_S_start, frequency_R_start, Q_start

    def fitResonators(self, R_S, frequency_R, Q, frequencyWindow=None,
                      RShuntBound=None, freqBound=None, QBound=None,
                      RShuntScale=1, freqScale=1, QScale=1,
                      ImZoverF=None, ImZoverFBound=None, ImZoverFScale=1,
                      method=None, fitResidue='log_total', constraints=None,
                      disp=False, maxiter=None, maxfev=None, nDivide=False):
        '''
        Fit impedance table with resonators
        '''

        if frequencyWindow is not None:
            self.freqIndexes = np.where(
                (self.freqArray > frequencyWindow[0]) * (self.freqArray < frequencyWindow[1]))[0]
        else:
            self.freqIndexes = np.arange(len(self.freqArray))

        impedanceTable = np.hstack((np.hstack((self.freqArrayFitted.reshape((1, len(self.freqArrayFitted))).T,
                                               self.impedanceFitted.real.reshape((1, len(self.impedanceFitted.real))).T)),
                                    self.impedanceFitted.imag.reshape((1, len(self.impedanceFitted.imag))).T))

        R_S = np.array(R_S, dtype=float, ndmin=1)
        frequency_R = np.array(frequency_R, dtype=float, ndmin=1)
        Q = np.array(Q, dtype=float, ndmin=1)

        if (RShuntBound is not None) or (freqBound is not None) or (QBound is
                                                                    not None):
            bounds_min = np.zeros((len(R_S) * 3, 1))
            bounds_max = np.zeros((len(R_S) * 3, 1))

            for resIndex in range(0, int(R_S.shape[0])):
                if isinstance(RShuntBound, float):
                    bounds_min[3 * resIndex, 0] = (1 - np.sign(
                        R_S[resIndex]) * RShuntBound) * R_S[resIndex] / RShuntScale
                    bounds_max[3 * resIndex, 0] = (1 + np.sign(
                        R_S[resIndex]) * RShuntBound) * R_S[resIndex] / RShuntScale
                elif isinstance(RShuntBound, list):
                    bounds_min[3 * resIndex, 0] = RShuntBound[0] / RShuntScale
                    bounds_max[3 * resIndex, 0] = RShuntBound[1] / RShuntScale
                else:
                    bounds_min[3 * resIndex, 0] = None
                    bounds_max[3 * resIndex, 0] = None

                if isinstance(freqBound, float):
                    bounds_min[3 * resIndex + 1,
                               0] = (1 - freqBound) * frequency_R[resIndex] / freqScale
                    bounds_max[3 * resIndex + 1,
                               0] = (1 + freqBound) * frequency_R[resIndex] / freqScale
                elif isinstance(freqBound, list):
                    bounds_min[3 * resIndex + 1, 0] = freqBound[0] / freqScale
                    bounds_max[3 * resIndex + 1, 0] = freqBound[1] / freqScale
                else:
                    bounds_min[3 * resIndex + 1, 0] = None
                    bounds_max[3 * resIndex + 1, 0] = None

                if isinstance(QBound, float):
                    bounds_min[3 * resIndex + 2,
                               0] = (1 - QBound) * Q[resIndex] / QScale
                    if bounds_min[3 * resIndex + 2, 0] * QScale < 0.5:
                        bounds_min[3 * resIndex + 2, 0] = 0.5 / QScale
                    bounds_max[3 * resIndex + 2,
                               0] = (1 + QBound) * Q[resIndex] / QScale
                elif isinstance(QBound, list):
                    bounds_min[3 * resIndex + 2, 0] = QBound[0] / QScale
                    bounds_max[3 * resIndex + 2, 0] = QBound[1] / QScale
                else:
                    bounds_min[3 * resIndex + 2, 0] = None
                    bounds_max[3 * resIndex + 2, 0] = None

            bounds_array = np.hstack((bounds_min, bounds_max))
            bounds_array = np.where(np.isnan(bounds_array), None, bounds_array)
            bounds = tuple(map(tuple, bounds_array))

            if ImZoverF is not None:
                if isinstance(ImZoverFBound, float):
                    bounds_min = (1 - np.sign(ImZoverF) *
                                  ImZoverFBound) * ImZoverF / ImZoverFScale
                    bounds_max = (1 + np.sign(ImZoverF) *
                                  ImZoverFBound) * ImZoverF / ImZoverFScale
                elif isinstance(ImZoverFBound, list):
                    bounds_min = ImZoverFBound[0] / ImZoverFScale
                    bounds_max = ImZoverFBound[1] / ImZoverFScale
                else:
                    bounds_min = None
                    bounds_max = None

                bounds += ((bounds_min, bounds_max),)

        else:
            bounds = None

        R_S /= RShuntScale
        frequency_R /= freqScale
        Q /= QScale

        startingParams = np.hstack((np.hstack((R_S.reshape(1, len(R_S)).T,
                                               frequency_R.reshape(1, len(frequency_R)).T)),
                                    Q.reshape(1, len(Q)).T))

        self.resonatorForFit = Resonators(startingParams[:, 0] * RShuntScale,
                                          startingParams[:, 1] * freqScale,
                                          startingParams[:, 2] * QScale)

        startingParams = startingParams.reshape(
            (int(startingParams.shape[0] * startingParams.shape[1]), ))

        if ImZoverF is not None:
            ImZoverF /= ImZoverFScale
            startingParams = np.append(startingParams, ImZoverF)

        if constraints is not None:
            method = None
            if constraints == 'positive_real':
                constraints = {'type': 'ineq',
                               'fun': self.fitResidualFunction,
                               'args': (impedanceTable,
                                        self.resonatorForFit,
                                        RShuntScale,
                                        freqScale,
                                        QScale,
                                        self.freqIndexes,
                                        positive_real_constraint,
                                        nDivide,
                                        ImZoverF,
                                        ImZoverFScale)}

        fittedParameters = minimize(self.fitResidualFunction, startingParams,
                                    args=(impedanceTable,
                                          self.resonatorForFit,
                                          RShuntScale,
                                          freqScale,
                                          QScale,
                                          self.freqIndexes,
                                          fitResidue,
                                          nDivide,
                                          ImZoverF,
                                          ImZoverFScale),
                                    bounds=bounds,
                                    method=method,
                                    constraints=constraints,
                                    options={'disp': disp,
                                             'maxfev': maxfev,
                                             'maxiter': maxiter})['x']

        if ImZoverF is None:
            fittedParameters = fittedParameters.reshape(
                (int(fittedParameters.shape[0] / 3), 3))
        else:
            ImZoverF_fit = fittedParameters[-1] * ImZoverFScale
            fittedParameters = fittedParameters[0:-1].reshape(
                (int(fittedParameters[0:-1].shape[0] / 3), 3))

        R_S_fit = fittedParameters[:, 0] * RShuntScale
        frequency_R_fit = np.abs(fittedParameters[:, 1] * freqScale)
        Q_fit = np.abs(fittedParameters[:, 2] * QScale)

        if ImZoverF is None:
            return R_S_fit, frequency_R_fit, Q_fit
        else:
            return R_S_fit, frequency_R_fit, Q_fit, ImZoverF_fit

    def fitMultiResonators(self, R_S, frequency_R, Q, frequencyWindow=None,
                           RShuntBound=None, freqBound=None, QBound=None,
                           RShuntScale=1, freqScale=1, QScale=1,
                           ImZoverF=None, ImZoverFBound=None, ImZoverFScale=1,
                           method=None, fitResidue='log_total',
                           constraints=None,
                           n_res_max=None, disp=False,
                           maxiter=None, maxfev=None):
        '''
        Fit impedance table with resonators
        '''

        if n_res_max is None:
            n_res_max = len(R_S)

        R_S_loop = np.array(R_S)
        frequency_R_loop = np.array(frequency_R)
        Q_loop = np.array(Q)

#         plt.ioff()

        for index_fit in range(n_res_max):

            fitted_parameters = self.fitResonators(
                R_S_loop, frequency_R_loop, Q_loop,
                frequencyWindow=frequencyWindow,
                RShuntBound=RShuntBound, freqBound=freqBound, QBound=QBound,
                RShuntScale=RShuntScale, freqScale=freqScale, QScale=QScale,
                ImZoverF=ImZoverF, ImZoverFBound=ImZoverFBound,
                ImZoverFScale=ImZoverFScale,
                method=method, constraints=constraints,
                fitResidue=fitResidue,
                disp=disp,
                maxiter=maxiter,
                maxfev=maxfev)

            if ImZoverF is None:
                (R_S_fit, frequency_R_fit, Q_fit) = fitted_parameters
            else:
                (R_S_fit, frequency_R_fit, Q_fit, ImZoverF_fit) = fitted_parameters

            if (index_fit < n_res_max - 1) and (method is not 'Nelder-Mead'):

                if fitResidue[-5:] == 'total':
                    discrepancy = np.abs(
                        self.impedanceFitted) - np.abs(self.resonatorForFit.impedance)
                else:
                    discrepancy = self.impedanceFitted - self.resonatorForFit.impedance

                max_peak = np.sign(discrepancy.real[np.abs(discrepancy) == np.max(
                    np.abs(discrepancy))]) * np.max(np.abs(discrepancy))
                max_peak_freq = self.freqArrayFitted[np.abs(
                    discrepancy) == np.max(np.abs(discrepancy))]

                R_S_loop = np.append(R_S_fit, max_peak)
                frequency_R_loop = np.append(frequency_R_fit, max_peak_freq)
                Q_loop = np.append(Q_fit, np.mean(Q))
                if ImZoverF is None:
                    ImZoverF_loop = None
                else:
                    ImZoverF_loop = ImZoverF_fit

                fit_parameters_test = self.fitResonators(
                    R_S_loop, frequency_R_loop, Q_loop,
                    frequencyWindow=frequencyWindow,
                    RShuntBound=0., freqBound=0., QBound=QBound,
                    RShuntScale=RShuntScale, freqScale=freqScale,
                    QScale=QScale, method=method, constraints=constraints,
                    fitResidue=fitResidue,
                    ImZoverF=ImZoverF_loop, ImZoverFBound=ImZoverFBound,
                    ImZoverFScale=ImZoverFScale,
                    disp=disp,
                    maxiter=maxiter,
                    maxfev=maxfev)

                if ImZoverF is None:
                    (R_S_fit_test, frequency_R_fit_test,
                     Q_fit_test) = fit_parameters_test
                else:
                    (R_S_fit_test, frequency_R_fit_test, Q_fit_test,
                     ImZoverF_fit_test) = fit_parameters_test

                R_S_loop = R_S_fit_test
                frequency_R_loop = frequency_R_fit_test
                Q_loop = Q_fit_test

#                 print(max_peak_freq, max_peak)
#                 plt.figure('test')
#                 plt.clf()
#                 plt.plot(self.freqArrayFitted, discrepancy.real)
#                 plt.plot(self.freqArrayFitted, discrepancy.imag)
#                 plt.show()

        if ImZoverF is None:
            return R_S_fit, frequency_R_fit, Q_fit
        else:
            return R_S_fit, frequency_R_fit, Q_fit, ImZoverF_fit

    def fitResidualFunction(self, params, *inputImpedance):
        '''
        The residual function for fits
        '''

        (impedanceTable,
         resonatorForFit,
         RShuntScale,
         freqScale,
         QScale,
         freqIndexes,
         fitResidue,
         nDivide,
         ImZoverF,
         ImZoverFScale) = inputImpedance

        if ImZoverF is None:
            params = params.reshape((int(params.shape[0] / 3), 3))
        else:
            params_ImZoverF = params[-1] * ImZoverFScale
            params = params[0:-1].reshape((int(params[0:-1].shape[0] / 3), 3))

        resonatorForFit.R_S = params[:, 0] * RShuntScale
        resonatorForFit.frequency_R = np.abs(params[:, 1] * freqScale)
        resonatorForFit.Q = np.abs(params[:, 2] * QScale)
        resonatorForFit.Q[resonatorForFit.Q < 0.5] = 0.5

        resonatorForFit.imped_calc(self.freqArray[freqIndexes])

        if ImZoverF is not None:
            resonatorForFit.impedance.imag += params_ImZoverF * \
                self.freqArray[freqIndexes]

        if isinstance(fitResidue, str):
            residue = self.fitResiduals[fitResidue](
                impedanceTable, resonatorForFit)
        else:
            residue = fitResidue(impedanceTable, resonatorForFit)

        residue /= RShuntScale
        return residue

    def scanFitFunction(params, *InputParams):

        # divisor allows division by frequency to emphasise lower frequency
        # resonances
        impedanceTarget, resonatorForFit, resonator_number, allR, allF, allQ, fitTarget, divisor = InputParams

        # Option to scan only R&Q, only Frequency or all 3 parameters
        if fitTarget == 'R/Q':
            allR[resonator_number] = params[0]
            allQ[resonator_number] = params[1]
        elif fitTarget == 'Freq':
            allF[resonator_number] = params[0]
        elif fitTarget == 'All':
            allR[resonator_number] = params[0]
            allF[resonator_number] = params[1]
            allQ[resonator_number] = params[1]
        else:
            raise RuntimeError("Fit target not recognised")

        resonatorForFit.R_S = allR * RShuntScale
        resonatorForFit.frequency_R = allF * freqScale
        resonatorForFit.Q = allQ * QScale
        resonatorForFit.Q[resonatorForFit.Q < 0.5] = 0.5
        resonatorForFit.imped_calc(impedanceTable[0])

        sourceImag = impedanceTarget[2] / divisor
        fitImag = resonatorForFit.impedance.imag / divisor

        sourceReal = impedanceTarget[1] / divisor
        fitReal = resonatorForFit.impedance.real / divisor

        return np.sum(np.abs(sourceImag - fitImag) + np.abs(sourceReal - fitReal))

    @property
    def freqArrayFitted(self):

        if hasattr(self, 'freqIndexes'):
            return self.freqArray[self.freqIndexes]
        else:
            return self.freqArray

    @property
    def impedanceFitted(self):

        if hasattr(self, 'freqIndexes'):
            return self.impedance[self.freqIndexes]
        else:
            return self.impedance


def positive_real_constraint(impedanceTable, resonatorForFit):
    '''
    Constraint function checking the resistive impedance is positive
    '''

    return(np.sum(resonatorForFit.impedance.real >= 0) - len(resonatorForFit.impedance.real) + 1)


def minmax_location(x, f):
    '''
    *Function to locate the minima and maxima of the f(x) numerical function.*
    '''

    f_derivative = np.diff(f)

    f_derivative_second = np.diff(f_derivative)

    warnings.filterwarnings("ignore")
    f_derivative_zeros = np.unique(np.append(np.where(
        f_derivative == 0), np.where(f_derivative[1:] / f_derivative[0:-1] < 0)))

    min_x_position = x[f_derivative_zeros[f_derivative_second[f_derivative_zeros] > 0] + 1]
    max_x_position = x[f_derivative_zeros[f_derivative_second[f_derivative_zeros] < 0] + 1]

    min_values = f[f_derivative_zeros[f_derivative_second[f_derivative_zeros] > 0] + 1]
    max_values = f[f_derivative_zeros[f_derivative_second[f_derivative_zeros] < 0] + 1]

    warnings.filterwarnings("default")

    return [min_x_position, max_x_position], [min_values, max_values]
