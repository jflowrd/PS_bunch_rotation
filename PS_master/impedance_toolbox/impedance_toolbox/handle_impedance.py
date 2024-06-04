"""
**Module used to import the SPS impedance scenario and apply wanted modifications**

:Authors: **Joel Repond**, **Alexandre Lasheen**
"""

# Imports
import os
import inspect
import numpy as np
import matplotlib.pyplot as plt
from .miniblond.impedances.impedance_sources import Resonators, InputTable, TravelingWaveCavity


class handleImpedance(object):
    """Fundamental class defining how to interact with the impedance tables

    This class contains the main methods to import the data from the impedance
    folder and to act on it (damping, frequency shifting)

    It contains a dictionary called self.table_impedance which contains:
        -- the .keys() of the dictionary are based on the file name and arborescence
        -- the ['type'] of impedance ('resonator', 'twc', 'inputtable')
        -- the data depending in the type of impedance

    The three types of impedance used are as followed:
        -- 'resonator': resonating impedance defined by fr, Rsh and Q.
        -- 'inputtable': raw data table (from CST for example), three column
            containing fr, the frequency, ReZ the real part of Z and ImZ the
            imaginary part of the impedance.
        -- 'twc': theoretical travelling wave cavity impedance, containing fr,
            Rsh and alpha the time constant [1]_.

    Attributes
    ----------
    table_impedance : dict
        contains the impedance data.
    impedanceFolder : str
        root of the impedance folder.

    References
    ----------
    .. [1]  "The SPS acceleration system travelling wave drif-tube structure for
            the CERN SPS", G. D�me, CERN-SPS-ARF-77-11, 1977

    Examples
    --------
    >>> from blabla import blibli
    >>>
    """

# Init method
# ------------

    def __init__(self, **kwargs):
        self.table_impedance = dict()

        if 'folder' in kwargs:
            self.impedanceFolder = kwargs.get('folder')
        else:
            fold = os.path.dirname(
                os.path.abspath(
                    inspect.getfile(inspect.currentframe())))
            self.impedanceFolder = fold + '/../../impedance'


# Importation method
# -------------------

    def check_Q(self, tableimport):
        """
        This method checks (for a resonator) if Q < 0.5 and rescals everything
        accordingly to 0.5

        Parameters
        ----------
        tableimport : np.array
            imported 3D array containing [fr, Rsh, Q]

        Returns
        --------
        tableimport : np.array
            rescaled array if necessary
        """

        if len(tableimport[:, 2][tableimport[:, 2] <= 0.5]) > 0:
            tableimport[:, 1][tableimport[:, 2] <= 0.5] /= tableimport[:, 2][
                tableimport[:, 2] <= 0.5] / 0.5
            tableimport[:, 2][tableimport[:, 2] <= 0.5] = 0.5
            return tableimport
        else:
            return tableimport

    def importWakeFromCST(self, filestr, unitFreq=1E6, ZFactor=1.,
                          debug=False):
        """
        This method import the data from the CST wakefield data and puts it
        into the table_impedance.

        NOTE: ASSUMING UNITS FOR CST ARE mm MHz ns K

        Parameters
        ----------
        filestr : str
            path of the corresponding file, relative to impedanceFolder, this
            str is used a key in the table_impedance dictionary to easily find
            your impedance if necessary.
        unitFreq (1E6): float
            the frequency must be in Hz, are assuming CST is giving frequencies 
            in MHz and then we import it as Hz (CST User Defined)
        ZFactor (1.): float
            used to modify the real and imaginary part of the impedance of the
            entire file loaded
        """
        def debug_print(*args, debug=False):
            if debug is True:
                for it, elem in enumerate(args):
                    print(elem)

        def find_files_loc(directory, fileEnding):
            CSVFileLoc = []
            # Search for all .csv files in directory's sub directorys
            for root, dirs, files in os.walk(directory):
                for file in files:
                    if file.endswith(fileEnding):
                        CSVFileLoc.append(os.path.join(root, file))
                        debug_print('Meas. File \'{:s}\'Imported'.format(file))
            return CSVFileLoc

        def import_re_im(file, header):
            if np.size(header) != 4:
                print('ERROR: This new header should contain two columns: Freq | complexNum along with their data types')
                return
            tempHeader = ['freq', 're', 'im']
            importedReIM = np.genfromtxt(file,
                                         dtype=None,
                                         skip_header=0,
                                         names=tempHeader)
        #    print(np.size(importedReIM),header)
            importedComplex = np.zeros(np.size(importedReIM), dtype=header)
            importedComplex[header[0][0]] = importedReIM['freq']
            importedComplex[header[1][0]] = importedReIM['re']+1j*importedReIM['im']
            return importedComplex
        sim_direc = self.impedanceFolder + '/' + filestr
        txt_files = find_files_loc(sim_direc, ".txt")

        numEfiles = 0
        for it, elem in enumerate(txt_files):
            fName = elem[elem.rfind('\\')+1:elem.find('.txt')]

    #        debug_print("fName: {:}".format(fName))
            debug_print("FileLoc", elem)
            wakeHeader = [('freq', np.float32),
                          ('impedX', np.complex64),
                          ('impedY', np.complex64),
                          ('impedZ', np.complex64)]
            # These keys are the exact file names that are output from CST Post Processing
            cst_file_keys = ['Particle Beams_ParticleBeam1_Wake impedance [WakeSpectrum]_X',
                             'Particle Beams_ParticleBeam1_Wake impedance [WakeSpectrum]_Y',
                             'Particle Beams_ParticleBeam1_Wake impedance [WakeSpectrum]_Z']
            cst_file_keys_bis = ['Particle Beams_ParticleBeam1_Wake impedance (WakeSpectrum)_X',
                                 'Particle Beams_ParticleBeam1_Wake impedance (WakeSpectrum)_Y',
                                 'Particle Beams_ParticleBeam1_Wake impedance (WakeSpectrum)_Z']
            if ('impedance_X' in fName) or (cst_file_keys[0] in fName) or (cst_file_keys_bis[0] in fName):
                debug_print('Importing Impedance X from: {:}'.format(fName),
                            debug=debug)

                wakeImpedX = import_re_im(elem, [wakeHeader[0], wakeHeader[1]])
                numEfiles = numEfiles+1
            elif ('impedance_Y' in fName) or (cst_file_keys[1] in fName) or (cst_file_keys_bis[1] in fName):
                debug_print('Importing Impedance Y from: {:}'.format(fName),
                            debug=debug)

                wakeImpedY = import_re_im(elem, [wakeHeader[0], wakeHeader[2]])
                numEfiles = numEfiles + 1
            elif ('impedance_Z' in fName) or (cst_file_keys[2] in fName) or (cst_file_keys_bis[2] in fName):
                debug_print('Importing Impedance Z from: {:}'.format(fName),
                            debug=debug)

                wakeImpedZ = import_re_im(elem, [wakeHeader[0], wakeHeader[3]])
                numEfiles = numEfiles + 1
            else:
                debug_print('ERROR: Not an Wakefield related file: {:}'.format(fName))
        if numEfiles != 3:
            print('ERROR: {:} Files Found, should have 3 Wakefield files in the directory...'.format(numEfiles))
            print('Check the CST Post Processing to Confirm')
            return

        wakeData = np.zeros(np.size(wakeImpedX), dtype=wakeHeader)
        wakeData['freq'] = wakeImpedX['freq']*unitFreq
        wakeData['impedX'] = wakeImpedX['impedX']
        wakeData['impedY'] = wakeImpedY['impedY']
        wakeData['impedZ'] = wakeImpedZ['impedZ']

        frequency = wakeData['freq']
        RealZ = np.real(wakeData['impedZ'])*ZFactor
        ImagZ = np.imag(wakeData['impedZ'])*ZFactor

        self.table_impedance[filestr] = dict()

        self.table_impedance[filestr]['fr'] = frequency
        self.table_impedance[filestr]['ReZ'] = RealZ
        self.table_impedance[filestr]['ImZ'] = ImagZ
        self.table_impedance[filestr]['type'] = 'inputtable'

    def importEigenFromCST(self, filestr, unitFreq=1E6, unitRsh=1,
                           RshFactor=1., QFactor=1.,
                           debug=False):
        """
        This method import the data for resonators from a the CST Eigenmode
        post-processing folder (export)

        NOTE: ASSUMING UNITS FOR CST ARE mm MHz ns K

        0.5 factor is used to account for CST's different definition of Rshunt

        Parameters
        ----------
        filestr : str
            path of the corresponding file, relative to impedanceFolder, this
            str is used a key in the table_impedance dictionary to easily find
            your impedance if necessary.
        unitFreq (1E6): float
            the frequency must be in Hz
        unitRsh (1.): float
            the shunt impedance must be in Ohm
        RshFactor (1.): float
            used to modify the shunt impedance of the entire file loaded
        QFactor (1.): float
            used to modify the quality factor Q of the entire file loaded

        """
        def debug_print(*args, debug=False):
            if debug is True:
                for it, elem in enumerate(args):
                    print(elem)

        def find_files_loc(directory, fileEnding):
            CSVFileLoc = []
            # Search for all .csv files in directory's sub directorys
            for root, dirs, files in os.walk(directory):
                for file in files:
                    if file.endswith(fileEnding):
                        CSVFileLoc.append(os.path.join(root, file))
                        debug_print('Meas. File \'{:s}\'Imported'.format(file))
            return CSVFileLoc
        sim_direc = self.impedanceFolder + '/' + filestr
        ASCIIFileLoc = find_files_loc(sim_direc, ".txt")

        numEfiles = 0
        for it, elem in enumerate(ASCIIFileLoc):
            fName = elem[elem.rfind('\\') + 1:elem.find('.txt')]
            eigenHeader = [('mode', np.int32),
                           ('resFreq', np.float32),
                           ('Q', np.float32),
                           ('RShunt', np.float32),
                           ('RoQ', np.float32)]

            if 'Frequency' in fName:
                debug_print('Importing Res Freq from file: {:}'.format(fName),
                            debug=debug)
                eigenFreq = np.genfromtxt(elem,
                                          dtype=None,
                                          skip_header=0,
                                          names=[eigenHeader[0][0],
                                                 eigenHeader[1][0]])
                if eigenFreq.size == 1:
                    eigenFreq = {eigenHeader[0][0]: 1,
                                 eigenHeader[1][0]: np.array(eigenFreq, dtype=float)}
                numEfiles = numEfiles + 1
            elif 'Q-Factor (Perturbation)' in fName:
                debug_print('Importing Q from file: {:}'.format(fName),
                            debug=debug)
                eigenQ = np.genfromtxt(elem,
                                       dtype=None,
                                       skip_header=0,
                                       names=[eigenHeader[0][0],
                                              eigenHeader[2][0]])
                if eigenQ.size == 1:
                    eigenQ = {eigenHeader[0][0]: 1,
                              eigenHeader[2][0]: np.array(eigenQ, dtype=float)}
                numEfiles = numEfiles+1
            elif 'R over Q' in fName:
                debug_print('Importing R over Q from file: {:}'.format(fName),
                            debug=debug)
                eigenRoQ = np.genfromtxt(elem,
                                         dtype=None,
                                         skip_header=0,
                                         names=[eigenHeader[0][0],
                                                eigenHeader[4][0]])
                if eigenRoQ.size == 1:
                    eigenRoQ = {eigenHeader[0][0]: 1,
                                eigenHeader[4][0]: np.array(eigenRoQ, dtype=float)}
                numEfiles = numEfiles+1
            elif 'Shunt Impedance' in fName:
                debug_print('Importing R Shunt from file: {:}'.format(fName),
                            debug=debug)
                eigenRShunt = np.genfromtxt(elem,
                                            dtype=None,
                                            skip_header=0,
                                            names=[eigenHeader[0][0],
                                                   eigenHeader[3][0]])
                if eigenRShunt.size == 1:
                    eigenRShunt = {eigenHeader[0][0]: 1,
                                   eigenHeader[3][0]: np.array(eigenRShunt, dtype=float)}
                numEfiles = numEfiles+1

        if numEfiles < 3:
            print('ERROR: not all eigenmode files are present in the directory')
            print('Check the CST Post Processing to Confirm 3 Files At Least')
            print(ASCIIFileLoc)
            return

        if 'eigenRShunt' not in locals():
            eigenRShunt = {}
            eigenRShunt['RShunt'] = eigenRoQ['RoQ'] * eigenQ['Q']

        self.table_impedance[filestr] = dict()

        self.table_impedance[filestr]['fr'] = eigenFreq['resFreq']*unitFreq
        self.table_impedance[filestr]['Rsh'] = 0.5*eigenRShunt['RShunt']*unitRsh * RshFactor
        self.table_impedance[filestr]['Q'] = eigenQ['Q']*QFactor
        self.table_impedance[filestr]['type'] = 'resonator'

    def importFromCST(self, filename_list, ZFactor=1., RshFactor=1.,
                      QFactor=1., unitRsh=1., unitFreq=1E6, debug=False):
        '''
        Generic function to import impedance from CST (wake or eigen).
        Can treat list of impedance, and identifies wether it comes
        from wakefield or eigenmode simulations based on the filename
        '''

        if isinstance(filename_list, str):
            if filename_list[-4:] == 'wake':
                self.importWakeFromCST(filename_list, unitFreq=unitFreq,
                                       ZFactor=ZFactor, debug=debug)
            elif filename_list[-5:] == 'eigen':
                self.importEigenFromCST(filename_list, unitFreq=unitFreq,
                                        unitRsh=unitRsh, RshFactor=RshFactor,
                                        QFactor=QFactor, debug=debug)

        elif isinstance(filename_list, list):
            for filename in filename_list:
                if filename[-4:] == 'wake':
                    self.importWakeFromCST(filename, unitFreq=unitFreq,
                                           ZFactor=ZFactor, debug=debug)
                elif filename[-5:] == 'eigen':
                    self.importEigenFromCST(filename, unitFreq=unitFreq,
                                            unitRsh=unitRsh,
                                            RshFactor=RshFactor,
                                            QFactor=QFactor, debug=debug)

    def importResonatorFromFile(self, filestr, unitFreq=1e9, unitRsh=1e3,
                                RshFactor=1., QFactor=1., delimiter=None):

        """
        This method import the data for resonators from a file and store it in
        table_impedance

        Parameters
        ----------
        filestr : str
            path of the corresponding file, relative to impedanceFolder, this
            str is used a key in the table_impedance dictionary to easily find
            your impedance if necessary.
        unitFreq : float
            the frequency must be in Hz
        unitRsh : float
            the shunt impedance must be in Ohm
        RshFactor : float
            used to modify the shunt impedance of the entire file loaded
        QFactor : float
            used to modify the quality factor Q of the entire file loaded

        """

        if not isinstance(filestr, str):
            raise SystemError('The first argument of import_resonator_fromfile must be a string containing the name of the file')

        importedFile = np.atleast_2d(np.loadtxt(self.impedanceFolder + '/' +
                                                filestr, comments=['!',
                                                                   '%',
                                                                   '#'],
                                                delimiter=delimiter))
        importedFile[:, 0] *= unitFreq
        importedFile[:, 1] *= unitRsh * RshFactor
        importedFile[:, 2] *= QFactor
        self.check_Q(importedFile)

        self.table_impedance[filestr] = dict()

        self.table_impedance[filestr]['fr'] = importedFile[:, 0]
        self.table_impedance[filestr]['Rsh'] = importedFile[:, 1]
        self.table_impedance[filestr]['Q'] = importedFile[:, 2]
        self.table_impedance[filestr]['type'] = 'resonator'

    def importResonatorsFromList(self, thelist, thestr, unitFreq=1e9,
                                 unitRsh=1e3, RshFactor=1., QFactor=1.):
        """
        This method import an inputtable from a 1D list ([f_r, R_S, Q]) or
        2D list ([f_r_list, R_S_list, Q_list]) and store it in table_impedance

        Parameters
        ----------
        thelist : list
            The list containing the impedance table in the form [f_r, R_S, Q]
            or ([f_r_list, R_S_list, Q_list])
        thestr : str
            this str is used a key in the table_impedance dictionary to easily
            find your impedance if necessary.
        unitFreq : float
            the frequency must be in Hz
        RshFactor : float
            used to modify the shunt impedance of the entire file loaded
        QFactor : float
            used to modify the quality factor Q of the entire file loaded
        """

        importedList = np.array(thelist, ndmin=2, dtype=float)

        # Transpose if 2D list as input
        if importedList.shape[0] != 1:
            importedList = importedList.T

        importedList[:, 0] *= unitFreq
        importedList[:, 1] *= unitRsh * RshFactor
        importedList[:, 2] *= QFactor

        self.check_Q(importedList)

        self.table_impedance[thestr] = dict()

        self.table_impedance[thestr]['fr'] = importedList[:, 0]
        self.table_impedance[thestr]['Rsh'] = importedList[:, 1]
        self.table_impedance[thestr]['Q'] = importedList[:, 2]
        self.table_impedance[thestr]['type'] = 'resonator'

    def importInputTableFromFile(self, filestr, unitFreq=1., ZFactor=1.):

        """
        This method import the data from an inputtable in a file and store it in
        table_impedance

        Parameters
        ----------
        filestr : str
            path of the corresponding file, relative to impedanceFolder, this
            str is used a key in the table_impedance dictionary to easily find
            your impedance if necessary.
        unitFreq : float
            the frequency must be in Hz
        ZFactor : float
            used to modify the real and imaginary part of the impedance of the
            entire file loaded
        """

        if not isinstance(filestr, str):
            raise SystemError('The first argument of import_inputtable must ' +
                              'be a string containing the name of the file')

        importedFile = np.loadtxt(self.impedanceFolder + '/' +
                                  filestr, comments=['!', '%', '#'])
        importedFile[:, 0] *= unitFreq
        importedFile[:, 1] *= ZFactor
        importedFile[:, 2] *= ZFactor

        self.table_impedance[filestr] = dict()

        self.table_impedance[filestr]['fr'] = importedFile[:, 0]
        self.table_impedance[filestr]['ReZ'] = importedFile[:, 1]
        self.table_impedance[filestr]['ImZ'] = importedFile[:, 2]
        self.table_impedance[filestr]['type'] = 'inputtable'

    def importInputTableFromList(self, thelist, thestr, unitFreq=1.,
                                 ZFactor=1.):
        """
        This method import an inputtable from a list ([f,ReZ,ImZ]) and store it in
        table_impedance

        Parameters
        ----------
        thelist : list
            The list containing the impedance table in the form [freq, ReZ, ImZ]
        thestr : str
            this str is used a key in the table_impedance dictionary to easily find
            your impedance if necessary.
        unitFreq : float
            the frequency must be in Hz
        ZFactor : float
            used to modify the real and imaginary part of the impedance of the
            entire file loaded
        """

#        if not isinstance(filestr,str):
#            raise SystemError('The first argument of import_inputtable must be a string containing the name of the file')

#        importedFile =  np.loadtxt(self.impedanceFolder + '/' + filestr, comments=['!','%', '#'])
        frequency = thelist[0]*unitFreq
        RealZ = thelist[1]*ZFactor
        ImagZ = thelist[2]*ZFactor

        self.table_impedance[thestr] = dict()

        self.table_impedance[thestr]['fr'] = frequency
        self.table_impedance[thestr]['ReZ'] = RealZ
        self.table_impedance[thestr]['ImZ'] = ImagZ
        self.table_impedance[thestr]['type'] = 'inputtable'

    def importTWCFromFile(self, filestr, unitFreq=1e9, unitRsh=1e3,
                          unitAlpha=1e-6, RshFactor=1.):

        """
        This method import the data for a travelling wave cavity from a file
        and store it in table_impedance

        Parameters
        ----------
        filestr : str
            path of the corresponding file, relative to impedanceFolder, this
            str is used a key in the table_impedance dictionary to easily find
            your impedance if necessary.
        unitFreq : float
            the frequency must be in Hz
        unitRsh : float
            the shunt impedance must be in Ohm
        unitAlpha : float
            alpha must be in sec
        RshFactor : float
            used to modify the shunt impedance of the entire file loaded
        """

        if not isinstance(filestr, str):
            raise SystemError('The first argument of import_inputtable must ' +
                              'be a string containing the name of the file')

        importedFile = np.atleast_2d(np.loadtxt(self.impedanceFolder + '/' +
                                                filestr, comments=['!', '%',
                                                                   '#']))
        importedFile[:, 0] *= unitFreq
        importedFile[:, 1] *= unitRsh * RshFactor
        importedFile[:, 2] *= unitAlpha

        self.table_impedance[filestr] = dict()

        self.table_impedance[filestr]['fr'] = importedFile[:, 0]
        self.table_impedance[filestr]['Rsh'] = importedFile[:, 1]
        self.table_impedance[filestr]['alpha'] = importedFile[:, 2]
        self.table_impedance[filestr]['type'] = 'twc'

    def importImZ_over_f(self, thestr, ImZ_over_f):

        self.table_impedance[thestr] = dict()
        self.table_impedance[thestr]['impedance'] = ImZ_over_f
        self.table_impedance[thestr]['type'] = 'ImZ/f'


# Method to modify the loaded data
# ---------------------------------

    def move_fr_resonatorOrTWC(self, key, freq, delta_freq):

        """
        Shift the frequency freq by delta_freq in the data corresponding to key

        if freq is a list [f1,f2], all the component between f1 and f2 are shifted
        by delta_freq

        Parameters
        ----------
        key : str
            corresponding impedance data
        freq : float or list
            frequency or boundary between which the frequencies will be shifted.
            The frequency must be in Hz.
        delta_freq : float
            frequency shift.
        """

        frequencyToMove_array = self.table_impedance[key]['fr']

        if isinstance(freq, float):
            frequencyToMove_array[frequencyToMove_array == freq] += delta_freq

        elif isinstance(freq, list):
            cond1 = frequencyToMove_array >= freq[0]
            cond2 = frequencyToMove_array <= freq[1]

            frequencyToMove_array[cond1*cond2] += delta_freq

    def set_fr_resonatorOrTWC(self, key, freq, bounds=None):

        """
        Set the frequency to the given freq value in the data corresponding to
        key

        if freq is a list [f1,f2], all the component between f1 and f2 are
        set to freq

        Parameters
        ----------
        key : str
            corresponding impedance data
        freq : float
            Frequency in Hz
        bounds : float or list
            frequency or boundary between which the frequencies will be shifted.
            The frequency must be in Hz. If None all modes are set to freq
        """

        frequencyToMove_array = self.table_impedance[key]['fr']

        if isinstance(bounds, float):
            frequencyToMove_array[frequencyToMove_array == bounds] = freq

        elif isinstance(bounds, list):
            cond1 = frequencyToMove_array >= bounds[0]
            cond2 = frequencyToMove_array <= bounds[1]

            frequencyToMove_array[cond1*cond2] = freq

        elif bounds is None:
            frequencyToMove_array[:] = freq

    def damp_R_resonatorOrTWC(self, key, freq, R_factor):

        """
        Multiply the shunt impedance in the data corresponding to key by R_factor

        if freq is a list [f1,f2], all the component between f1 and f2 are modify

        Parameters
        ----------
        key : str
            corresponding impedance data
        freq : float or list
            frequency or boundary of freq between which Rsh will be multiplied.
            The frequency must be in Hz.
        R_factor : float
            factor by which the shunt impedance is multiplied.
        """

        frequencyToDamp_array = self.table_impedance[key]['fr']
        Rshunt_to_damp_array = self.table_impedance[key]['Rsh']

        if isinstance(freq, float):
            Rshunt_to_damp_array[frequencyToDamp_array == freq] *= R_factor

        elif isinstance(freq, list):
            cond1 = frequencyToDamp_array >= freq[0]
            cond2 = frequencyToDamp_array <= freq[1]

            Rshunt_to_damp_array[cond1*cond2] *= R_factor

    def damp_Q_resonatorOrTWC(self, key, freq, Q_factor):

        """
        Multiply the Q factor in the data corresponding to key by Q_factor

        if freq is a list [f1,f2], all the component between f1 and f2 are modify

        Parameters
        ----------
        key : str
            corresponding impedance data
        freq : float or list
            frequency or boundary of freq between which Q will be multiplied.
            The frequency must be in Hz.
        Q_factor : float
            factor by which the Q factor is multiplied.
        """

        frequencyToDamp_array = self.table_impedance[key]['fr']
        Q_to_damp_array = self.table_impedance[key]['Q']

        if isinstance(freq, float):
            Q_to_damp_array[frequencyToDamp_array == freq] *= Q_factor

        elif isinstance(freq, list):
            cond1 = frequencyToDamp_array >= freq[0]
            cond2 = frequencyToDamp_array <= freq[1]

            Q_to_damp_array[cond1*cond2] *= Q_factor

    def print_resonator_table(self, key):
        s = ''
        if self.table_impedance[key]['type'] == 'resonator':
            fr = self.table_impedance[key]['fr']/1e9
            Rsh = self.table_impedance[key]['Rsh']/1e3
            Q = self.table_impedance[key]['Q']

            # first line
            s += 'Resonator name: '+key+'\n\n'
            s += 'fr (GHz)\t Rsh (kOhm)\t Q\t\t R/Q(Ohm)\t Q/(pi*fr) (ns)\t #b coupled\n\n'

            for it in range(len(fr)):
                s += '%.3f \t\t %.3f \t %.3f \t %.3f \t %.3f \t%.3f\n' % (fr[it], Rsh[it], Q[it], Rsh[it]/Q[it]*1000, Q[it]/np.pi/fr[it], Q[it]/np.pi/fr[it]/25.)
        else:
            s += 'name: '+key+' ... not a resonator'
            print('the function print_resonator_table can be use for resonator type impedance only')
        return s


class impedance2blond(handleImpedance):

    """Class used to convert the impedance scenario to object usable by BLonD


    Attributes
    ----------
    table_impedance : dict
        Dictionary generated by the scenario class containing all the impedance
    wakeList : list
        list containing the BLonD impedance sources to be solved in time
    impedanceList : list
            list containing the BLonD impedance sources to be solved in freq
    impedanceListToPlot : list
            list containing all the BLonD impedance sources (to be plotted)

    Examples
    --------
    >>> from blabla import blibli
    >>>
    """

    def __init__(self, table_impedance):
        self.table_impedance = table_impedance

        self.wakeAndImpListProcess()

    def generateResonator(self, fr, Rsh, Q):

        """
        Generates a BLonD Resonator based on the input data

        Parameters
        ----------
        fr : float,array
            array (or float) of resonant frequency.
        Rsh : float,array
            array (or float) of shunt impedance.
        Q : float,array
            array (or float) of quality factor Q.

        Returns
        -------
        Resonators : Resonators
            BLonD resonators
        """

        return Resonators(Rsh, fr, Q)

    def get_impedance(self, key, freq_array):
        imp = self.table_impedance[key]['impedance']
        imp.imped_calc(freq_array)
        return imp.impedance

    def generateInputTable(self, fr, ReZ, ImZ):

        """
        Generates a BLonD InputTable based on the input data

        Parameters
        ----------
        fr : float,array
            array (or float) of resonant frequency.
        ReZ : float,array
            array (or float) of the real part of the impedance.
        Q : float,array
            array (or float) of the imaginary part of the impedance.

        Returns
        -------
        InputTable : InputTable
            BLonD InputTable
        """

        return InputTable(fr, ReZ, ImZ)

    def generateTWC(self, fr, Rsh, alpha):

        """
        Generates a BLonD TravelingWaveCavity based on the input data

        Parameters
        ----------
        fr : float,array
            array (or float) of resonant frequency.
        Rsh : float,array
            array (or float) of shunt impedance.
        alpha : float,array
            array (or float) of time factor alpha.

        Returns
        -------
        TravelingWaveCavity : TravelingWaveCavity
            BLonD TravelingWaveCavity
        """

        return TravelingWaveCavity(Rsh, fr, alpha)

    def wakeAndImpListProcess(self):

        """
        Reprocess the wake and impedance list from the table_impedance
        """

        self.wakeList = list()
        self.impedanceList = list()
        self.ImZ_over_f_List = list()
        self.impedanceListToPlot = list()

        for key in self.table_impedance:
            imp = self.table_impedance[key]
            if imp['type'] == 'resonator':
                imp['impedance'] = self.generateResonator(imp['fr'], imp['Rsh'], imp['Q'])
                self.wakeList.append(imp['impedance'])
                self.impedanceListToPlot.append(imp['impedance'])
            elif imp['type'] == 'inputtable':
                imp['impedance'] = self.generateInputTable(imp['fr'], imp['ReZ'], imp['ImZ'])
                self.impedanceList.append(imp['impedance'])
                self.impedanceListToPlot.append(imp['impedance'])
            elif imp['type'] == 'twc':
                imp['impedance'] = self.generateTWC(imp['fr'], imp['Rsh'], imp['alpha'])
                self.wakeList.append(imp['impedance'])
                self.impedanceListToPlot.append(imp['impedance'])
            elif imp['type'] == 'ImZ/f':
                self.ImZ_over_f_List.append(imp['impedance'])
            else:
                SystemExit('type of impedance not recognized when importing in BLonD')

    def plot_impedance_source(self, key, saveFile=None, figname=None, suptitle=None,
                       freqArray=None):

        if self.table_impedance[key]['type'] == 'resonator':
            self.plot_resonator(key, saveFile, figname, suptitle, freqArray)

        elif self.table_impedance[key]['type'] == 'inputtable':
            self.plot_inputtable(key, saveFile, figname, suptitle, freqArray)

    def plot_resonator(self, key, saveFile=None, figname=None, suptitle=None,
                       freqArray=None):

        if self.table_impedance[key]['type'] == 'resonator':
            fr = self.table_impedance[key]['fr']
            Rsh = self.table_impedance[key]['Rsh']
            Q = self.table_impedance[key]['Q']

            res2plot = self.generateResonator(fr, Rsh, Q)

            if freqArray is None:
                freqArray = np.linspace(0, 3*np.max(fr), num=10000)

            res2plot.imped_calc(freqArray)
            impedance = res2plot.impedance

            if figname is None:
                figname = key

            plt.figure(figname)
            plt.plot(freqArray, np.abs(impedance))
            plt.xlabel('f (Hz)')
            plt.ylabel('AbsZ ($\\Omega$)')
            if suptitle is None:
                pass
            else:
                plt.suptitle(key)
            if saveFile is not None:
                plt.savefig('./'+saveFile+'.png')

        else:
            print('the function plot_resonator can be use for resonator type impedance only')

    def plot_inputtable(self, key, saveFile=None, figname=None, suptitle=None,
                        freqArray=None):

        if self.table_impedance[key]['type'] == 'inputtable':
            fr = self.table_impedance[key]['fr']
            ReZ = self.table_impedance[key]['ReZ']
            ImZ = self.table_impedance[key]['ImZ']

            if freqArray is None:
                freqArray = fr
            else:
                ReZ = np.interp(freqArray, fr, ReZ)
                ImZ = np.interp(freqArray, fr, ImZ)

            impedance = ReZ + 1j * ImZ

            if figname is None:
                figname = key

            plt.figure(figname)
            plt.plot(freqArray, np.abs(impedance))
            plt.xlabel('f (Hz)')
            plt.ylabel('AbsZ ($\\Omega$)')
            if suptitle is None:
                pass
            else:
                plt.suptitle(key)
            if saveFile is not None:
                plt.savefig('./'+saveFile+'.png')

        else:
            print('the function plot_inputtable can be use for inputtable type impedance only')
