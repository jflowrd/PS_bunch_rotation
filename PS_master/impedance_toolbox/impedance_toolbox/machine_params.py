'''
Created on 9 May 2017

@author: alasheen
'''

import yaml
import numpy as np
from .miniblond.input_parameters.general_parameters import GeneralParameters
from .miniblond.input_parameters.rf_parameters import RFSectionParameters
from .miniblond.beams.beams import Beam
from .miniblond.beams.slices import Slices
from scipy.special import gamma, hyp0f1
from scipy.constants import e
from bisect import bisect_left
# from scipy.fftpack import next_fast_len


class MachineParameters(object):
    '''
    The MachineParameters object is based on BLonD objects and gathers all
    information about the machine, beam and impedance
    '''

    def __init__(self, inputFile=None, generalParameters=None, rfParameters=None,
                 bunchParameters=None, fillingScheme=None):
        '''
        Constructor
        '''

        if inputFile is not None:

            loadedFile = open(inputFile, 'r')
            loadedYAML = yaml.safe_load(loadedFile)
            loadedFile.close()

            generalParameters = loadedYAML['General parameters']
            rfParameters = loadedYAML['RF parameters']
            self.bunchParameters = loadedYAML['Bunch parameters']
            self.fillingScheme = loadedYAML['Filling scheme']

        elif (generalParameters is not None) and \
            (rfParameters is not None) and \
            (bunchParameters is not None) and \
                (fillingScheme is not None):

            self.bunchParameters = bunchParameters
            self.fillingScheme = fillingScheme

        else:
            raise IOError('The input for the MachineParameters object should either be an input file, or manually specifying general_parameters, rf_parameters, bunch_parameters and filling_scheme')

        self.generalParams = GeneralParameters(1, generalParameters['Circumference [m]'],
                                               1 /
                                               generalParameters['Transition gamma']**2.,
                                               generalParameters['Momentum [1e9 eV/c]'] * 1e9,
                                               generalParameters['Particle type'])

        self.rfParams = RFSectionParameters(self.generalParams,
                                            len(rfParameters['Harmonic numbers']),
                                            rfParameters['Harmonic numbers'],
                                            rfParameters['RF voltages'],
                                            rfParameters['RF phases'])

        self.beamCurrent = 0

        self.beamSpectrum = 0

    def generateBeamCurrent(self, nTurns, nPoints=None, resolutionTime=None, maxFreq=None, forceSingleTurn=False):

        self.nTurns = nTurns

        if (nPoints == None) and (resolutionTime == None) and (maxFreq == None):
            raise RuntimeError(
                'You should include at least one of those parameters in generateBeamCurrent function: nPoints or resolutionTime or maxFreq')

        if resolutionTime != None:
            nPoints = next_fast_len(
                int(self.generalParams.t_rev[0] * nTurns / resolutionTime) + 1)

        if maxFreq != None:
            resolutionTime = 1 / (2 * maxFreq)
            nPoints = next_fast_len(
                int(self.generalParams.t_rev[0] * nTurns / resolutionTime) + 1)

        self.exponentDistrib = self.bunchParameters['Binomial exponent']
        self.bunchLength4Sig = self.bunchParameters['Bunch length [ns]'] * 1e-9
        self.intensityPerBunch = self.bunchParameters['Intensity per bunch [1e10 ppb]'] * 1e10

        self.nBunches = self.fillingScheme['Number of bunches']
        self.nBatches = self.fillingScheme['Number of batches']
        self.nFills = self.fillingScheme['Number of fills']

        if isinstance(self.nBunches, int):
            self.nBunches = self.nBunches * np.ones(self.nFills, dtype='int')
        else:
            if len(self.nBunches) != self.nFills:
                raise RuntimeError('The length of the list of number of bunches ' +
                                   'does not match the number of fills !')
            self.nBunches = np.array(self.nBunches, dtype='int')

        if isinstance(self.nBatches, int):
            self.nBatches = self.nBatches * np.ones(self.nFills, dtype='int')
        else:
            if len(self.nBatches) != self.nFills:
                raise RuntimeError('The length of the list of number of batches ' +
                                   'does not match the number of fills !')
            self.nBatches = np.array(self.nBatches, dtype='int')

        self.bunchesSpacing = self.fillingScheme['Buckets between bunches']
        self.batchesSpacing = self.fillingScheme['Buckets between batches']
        self.fillsSpacing = self.fillingScheme['Buckets between fills']

        self.timeArray = np.linspace(
            0, self.generalParams.t_rev[0] * nTurns, nPoints)
        fullBunchLength = np.sqrt(
            3 + 2 * self.exponentDistrib) / 2. * self.bunchLength4Sig

        self.beamProfileNormalized = np.zeros(nPoints)

        bunch_position = self.generalParams.t_rev[0] / \
            self.rfParams.harmonic[0, 0] / 2

        for indexTurn in range(nTurns):
            bunch_position = indexTurn * \
                self.generalParams.t_rev[0] + self.generalParams.t_rev[0] / \
                self.rfParams.harmonic[0, 0] / 2
            for indexFill in range(self.nFills):
                for indexBatch in range(self.nBatches[indexFill]):
                    for indexBunch in range(self.nBunches[indexFill]):
                        self.beamProfileNormalized[(self.timeArray >= bunch_position - fullBunchLength / 2) * (self.timeArray <= bunch_position + fullBunchLength / 2)] += 2 * gamma(1.5 + self.exponentDistrib) / (fullBunchLength * np.sqrt(np.pi) * gamma(
                            1 + self.exponentDistrib)) * (1 - 4 * ((self.timeArray[(self.timeArray >= bunch_position - fullBunchLength / 2) * (self.timeArray <= bunch_position + fullBunchLength / 2)] - bunch_position) / fullBunchLength)**2.)**self.exponentDistrib
                        bunch_position += self.bunchesSpacing * \
                            self.generalParams.t_rev[0] / \
                            self.rfParams.harmonic[0, 0]
                    bunch_position += self.batchesSpacing * \
                        self.generalParams.t_rev[0] / self.rfParams.harmonic[0, 0] - \
                        self.bunchesSpacing * \
                        self.generalParams.t_rev[0] / \
                        self.rfParams.harmonic[0, 0]
                bunch_position += self.fillsSpacing * \
                    self.generalParams.t_rev[0] / self.rfParams.harmonic[0, 0] - \
                    self.batchesSpacing * \
                    self.generalParams.t_rev[0] / self.rfParams.harmonic[0, 0]
            if forceSingleTurn:
                break

        if not forceSingleTurn:
            self.nTotalBunches = np.sum(self.nBatches * self.nBunches * nTurns)
        else:
            self.nTotalBunches = np.sum(self.nBatches * self.nBunches)

        self.beamProfileNormalized /= (self.nTotalBunches)

        self.beamCurrent = e * self.intensityPerBunch * \
            self.nTotalBunches * self.beamProfileNormalized

#         self.beamObject = Beam(self.generalParams, 1, self.nTotalBunches*self.intensityPerBunch)
#
#         self.slicesObject = Slices(self.rfParams, self.beamObject, 1)
#
#         self.slicesObject.n_macroparticles = beamProfileNormalized
#         self.slicesObject.bin_centers = self.timeArray
#         self.slicesObject.n_slices = len(self.timeArray)
#         self.slicesObject.bin_size = self.timeArray[1]-self.timeArray[0]
#         self.slicesObject.cut_left = self.timeArray[0] - 0.5*self.slicesObject.bin_size
#         self.slicesObject.cut_right = self.timeArray[-1] + 0.5*self.slicesObject.bin_size
#         self.slicesObject.cuts_unit = 's'

    def generateBeamSpectrum(self):

        self.freqArray = np.fft.rfftfreq(
            len(self.beamCurrent), d=self.timeArray[1] - self.timeArray[0])

        self.beamSpectrum = np.fft.rfft(
            self.beamCurrent) * 2 / len(self.timeArray)
        self.beamSpectrumNormalized = np.fft.rfft(
            self.beamProfileNormalized) * 2 / len(self.timeArray)

        self.dcCurrent = self.beamSpectrum[0] / 2.

#         self.slicesObject.beam_spectrum = self.beamSpectrum
#         self.slicesObject.beam_spectrum_freq = self.freqArray

    @property
    def analyticalSpectrum(self, freqArray=None, amplitude=None):

        if hasattr(self, 'freqArray') and freqArray == None:
            freqArray = self.freqArray

        if not hasattr(self, 'freqArray') and freqArray == None:
            raise RuntimeError(
                'A frequency array is required to evaluate the analyticalSpectrum')

        if hasattr(self, 'dcCurrent'):
            zeroAmplitude = 2 * self.dcCurrent.real
        else:
            zeroAmplitude = 1.

        exponentDistrib = self.bunchParameters['Binomial exponent']
        bunchLength4Sig = self.bunchParameters['Bunch length [ns]'] * 1e-9

        fullBunchLength = np.sqrt(
            3 + 2 * exponentDistrib) / 2. * bunchLength4Sig

        return zeroAmplitude * hyp0f1(1.5 + exponentDistrib, -(np.pi * fullBunchLength * freqArray)**2. / 4.)


def next_fast_len(target):
    """
    Find the next fast size of input data to `fft`, for zero-padding, etc.

    SciPy's FFTPACK has efficient functions for radix {2, 3, 4, 5}, so this
    returns the next composite of the prime factors 2, 3, and 5 which is
    greater than or equal to `target`. (These are also known as 5-smooth
    numbers, regular numbers, or Hamming numbers.)

    Parameters
    ----------
    target : int
        Length to start searching from.  Must be a positive integer.

    Returns
    -------
    out : int
        The first 5-smooth number greater than or equal to `target`.

    Examples
    --------
    On a particular machine, an FFT of prime length takes 133 ms:

    >>> from scipy import fftpack
    >>> min_len = 10007  # prime length is worst case for speed
    >>> a = np.random.randn(min_len)
    >>> b = fftpack.fft(a)

    Zero-padding to the next 5-smooth length reduces computation time to
    211 us, a speedup of 630 times:

    >>> fftpack.helper.next_fast_len(min_len)
    10125
    >>> b = fftpack.fft(a, 10125)

    Rounding up to the next power of 2 is not optimal, taking 367 us to
    compute, 1.7 times as long as the 5-smooth size:

    >>> b = fftpack.fft(a, 16384)

    """
    hams = (8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 27, 30, 32, 36, 40, 45, 48,
            50, 54, 60, 64, 72, 75, 80, 81, 90, 96, 100, 108, 120, 125, 128,
            135, 144, 150, 160, 162, 180, 192, 200, 216, 225, 240, 243, 250,
            256, 270, 288, 300, 320, 324, 360, 375, 384, 400, 405, 432, 450,
            480, 486, 500, 512, 540, 576, 600, 625, 640, 648, 675, 720, 729,
            750, 768, 800, 810, 864, 900, 960, 972, 1000, 1024, 1080, 1125,
            1152, 1200, 1215, 1250, 1280, 1296, 1350, 1440, 1458, 1500, 1536,
            1600, 1620, 1728, 1800, 1875, 1920, 1944, 2000, 2025, 2048, 2160,
            2187, 2250, 2304, 2400, 2430, 2500, 2560, 2592, 2700, 2880, 2916,
            3000, 3072, 3125, 3200, 3240, 3375, 3456, 3600, 3645, 3750, 3840,
            3888, 4000, 4050, 4096, 4320, 4374, 4500, 4608, 4800, 4860, 5000,
            5120, 5184, 5400, 5625, 5760, 5832, 6000, 6075, 6144, 6250, 6400,
            6480, 6561, 6750, 6912, 7200, 7290, 7500, 7680, 7776, 8000, 8100,
            8192, 8640, 8748, 9000, 9216, 9375, 9600, 9720, 10000)

    if target <= 6:
        return target

    # Quickly check if it's already a power of 2
    if not (target & (target - 1)):
        return target

    # Get result quickly for small sizes, since FFT itself is similarly fast.
    if target <= hams[-1]:
        return hams[bisect_left(hams, target)]

    match = float('inf')  # Anything found will be smaller
    p5 = 1
    while p5 < target:
        p35 = p5
        while p35 < target:
            # Ceiling integer division, avoiding conversion to float
            # (quotient = ceil(target / p35))
            quotient = -(-target // p35)

            # Quickly find next power of 2 >= quotient
            p2 = 2**((quotient - 1).bit_length())

            N = p2 * p35
            if N == target:
                return N
            elif N < match:
                match = N
            p35 *= 3
            if p35 == target:
                return p35
        if p35 < match:
            match = p35
        p5 *= 5
        if p5 == target:
            return p5
    if p5 < match:
        match = p5
    return match
