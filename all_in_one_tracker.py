'''
An object including all the necessary to run BLonD
in one cell
'''

import os
import numpy as np
import h5py as hp
from IPython.display import Image, display

from blond.input_parameters.ring import Ring
from blond.input_parameters.rf_parameters import RFStation
from blond.beam.beam import Beam
from blond.beam.distributions import parabolic
from blond.beam.profile import Profile, CutOptions, FitOptions
from blond.trackers.tracker import RingAndRFTracker, FullRingAndRF
from blond.trackers.utilities import separatrix
from blond.plots.plot import Plot
from blond.plots.plot_impedance import plot_induced_voltage_vs_bin_centers, plot_impedance_vs_frequency, plot_wake_vs_time
from blond.monitors.monitors import BunchMonitor
from blond.impedances.impedance_sources import Resonators
from blond.impedances.impedance import InducedVoltageTime, InducedVoltageFreq, TotalInducedVoltage
from blond.llrf.beam_feedback import BeamFeedback

Rs_cavities = 6 * 480e3
f_r_cavities = 200.222e6
Q_cavities = 130

Rs_kickers = 5 * 23e3
f_r_kickers = 1e9
Q_kickers = 1

Rs_vacuum = 100 * 70e3
f_r_vacuum = 1.403e9
Q_vacuum = 1000


def blond_complete_tracker(n_turns, circumference, alpha_0, momentum, particle,
                           harmonic, voltage, phase,
                           n_macroparticles, intensity,
                           initial_bunch_length,
                           initial_bunch_position=None, initial_bunch_energy=None, initial_energy_spread=None,
                           cut_left=None, cut_right=None, n_slices=60,
                           enable_impedance=False,
                           time_or_freq_domain='time', frequency_resolution=None,
                           Rs_list=[Rs_cavities, Rs_kickers, Rs_vacuum],
                           f_r_list=[f_r_cavities, f_r_kickers, f_r_vacuum],
                           Q_list=[Q_cavities, Q_kickers, Q_vacuum],
                           PL_gain=None):

    ring = Ring(circumference, alpha_0, momentum, particle, n_turns=n_turns)

    rf_station = RFStation(ring, harmonic, voltage, phase)

    beam = Beam(ring, n_macroparticles, intensity)

    parabolic(ring, rf_station, beam, initial_bunch_length,
              bunch_position=initial_bunch_position,
              bunch_energy=initial_bunch_energy,
              energy_spread=initial_energy_spread)

    if cut_left is None:
        cut_left = 0
    if cut_right is None:
        cut_right = rf_station.t_rf[0, 0]
    cut_opt = CutOptions(
        cut_left=cut_left, cut_right=cut_right, n_slices=n_slices)
    fit_opt = FitOptions()

    profile = Profile(beam, CutOptions=cut_opt, FitOptions=fit_opt)
    profile.track()

    h5file = './data/saved_data'
    bunchmonitor = BunchMonitor(
        ring, rf_station, beam,
        h5file,
        Profile=profile)
    bunchmonitor.track()

    x_sep = np.linspace(float(cut_left), float(cut_right), 1000)
    y_sep = separatrix(ring, rf_station, x_sep)

    dt_plt = int(ring.n_turns / 10) * 2
    plots = Plot(ring, rf_station, beam,
                 dt_plt, ring.n_turns,
                 cut_left, cut_right,
                 -1.1 * np.nanmax(y_sep), 1.1 * np.nanmax(y_sep),
                 Profile=profile,
                 # show_plots=True,
                 h5file=h5file,
                 separatrix_plot=True)

    plots.alpha = 1

    phase_loop = None
    if PL_gain is not None:
        configuration = {'machine': 'SPS_F',
                         'PL_gain': PL_gain}
        phase_loop = BeamFeedback(ring, rf_station, profile, configuration)

    rf_tracker = RingAndRFTracker(rf_station, beam, BeamFeedback=phase_loop)
    full_tracker = FullRingAndRF([rf_tracker])

    if enable_impedance:
        res_list = []
        for idx_res in range(len(Rs_list)):
            res_list.append(Resonators(
                Rs_list[idx_res], f_r_list[idx_res], Q_list[idx_res]))

        if time_or_freq_domain == 'time':

            ind_voltage_time = InducedVoltageTime(
                beam, profile, res_list)

            total_ind_voltage = TotalInducedVoltage(
                beam, profile, [ind_voltage_time])

            plot_wake_vs_time(ind_voltage_time)

        elif time_or_freq_domain == 'freq':

            ind_voltage_freq = InducedVoltageFreq(
                beam, profile, res_list,
                frequency_resolution=frequency_resolution)

            total_ind_voltage = TotalInducedVoltage(
                beam, profile, [ind_voltage_freq])

            plot_impedance_vs_frequency(ind_voltage_freq)

        total_ind_voltage.induced_voltage_sum()

        plot_induced_voltage_vs_bin_centers(total_ind_voltage)

    for turn in range(ring.n_turns):

        if turn % int(ring.n_turns / 10) == 0:
            print('Tracking %d/%d turns' % (turn, ring.n_turns))

        full_tracker.track()
        profile.track()
        if enable_impedance:
            total_ind_voltage.track()
            if (turn + 1) % dt_plt == 0:
                plot_induced_voltage_vs_bin_centers(
                    total_ind_voltage, figure_index=turn + 1)
        plots.track()
        bunchmonitor.track()

    print('DONE!!!')

    display(Image(filename='./fig/bunch_length.png'))
    display(Image(filename='./fig/bunch_mean_position.png'))
    display(Image(filename='./fig/long_distr_%d.png' %
            (turn + 1)))
    if enable_impedance:

        if time_or_freq_domain == 'time':
            display(Image(filename='./fig/sum_wake_vs_table_times0.png'))
        elif time_or_freq_domain == 'freq':
            display(Image(filename='./fig/sum_imp_vs_freq_fft0.png'))

        display(Image(filename='./fig/induced_voltage_%d.png' %
                (turn + 1)))

    h5data = hp.File(h5file + '.h5', 'r')
    bl = np.array(h5data["/Beam/sigma_dt"][:], dtype=np.double)
    bl *= 4
    pos = np.array(h5data["/Beam/mean_dt"][:], dtype=np.double)
    h5data.close()

    return ring.cycle_time[0:-1], bl[0:-1], pos[0:-1]
