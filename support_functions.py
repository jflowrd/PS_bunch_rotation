'''
Support functions for Longitudinal hands on exercises on Tracking

The list of function to use is

- plot_phase_space_trajectory
- plot_phase_space_distribution
- synchrotron_tune
- separatrix
- run_animation
- oscillation_spectrum
- synchrotron_tune

'''

import numpy as np
import warnings
import matplotlib.pyplot as plt
from matplotlib import animation, rc
from IPython.display import HTML


def plot_phase_space_trajectory(dt_trajectory,
                                energy_trajectory,
                                time_sep=None,
                                separatrix_array=None,
                                figname=None,
                                draw_line=True,
                                xlim=None, ylim=None):
    """
    This support function can be used to plot the trajectory of a particle
    in the longitudinal phase space. The o marker is the starting coordinate,
    the * marker is the end coordinate.

    Parameters
    ----------
    dt_trajectory : np.array or list
        The time coordinates of a distribution of particles in [s]
    energy_trajectory : np.array or list
        The energy coordinates of a distribution of particles in [eV]
    time_sep : np.array
        The time of the separatrix array in [s], as output by the separatrix function
    separatrix_array : np.array
        The separatrix array in [eV], as output by the separatrix function
    figname : str, optional
        The name of your nice figure, by default None
    draw_line : bool, optional
        Hide the trajectory to plot only the start/end coordinates, by default True
    xlim : tuple, optional
        The limits in time for your nice plot in [s], e.g. (-t_rf/2, t_rf), by default None
    ylim : tuple, optional
        The limits in energy for your nice plot [eV], e.g. (-1e6, 1e6), by default None
    """

    dt_trajectory = np.array(dt_trajectory) * 1e9
    energy_trajectory = np.array(energy_trajectory) / 1e6
    plt.figure(figname, figsize=(8, 8))
    if figname is None:
        plt.clf()
    if draw_line:
        alpha = 1
    else:
        alpha = 0
    if dt_trajectory.ndim == 1:
        p = plt.plot(dt_trajectory, energy_trajectory,
                     alpha=alpha)
        plt.plot(dt_trajectory[0], energy_trajectory[0],
                 'o', color=p[0].get_color())
        plt.plot(dt_trajectory[-1], energy_trajectory[-1],
                 '*', color=p[0].get_color())
    else:
        for idx_part in range(dt_trajectory.shape[1]):
            p = plt.plot(dt_trajectory[:, idx_part],
                         energy_trajectory[:, idx_part],
                         alpha=alpha)
            plt.plot(dt_trajectory[0, idx_part],
                     energy_trajectory[0, idx_part],
                     'o', color=p[0].get_color())
            plt.plot(dt_trajectory[-1, idx_part],
                     energy_trajectory[-1, idx_part],
                     '*', color=p[0].get_color())

    if (time_sep is not None) and (separatrix_array is not None):
        plt.plot(time_sep * 1e9, separatrix_array / 1e6, 'r')

    plt.xlabel('Time $\\Delta t$ [ns]')
    plt.ylabel('Energy $\\Delta E$ [MeV]')
    plt.xlim(xlim)
    plt.ylim(ylim)


def plot_phase_space_distribution(
        time_coordinates,
        energy_coordinates,
        time_sep=None,
        separatrix_array=None,
        figname=None,
        xbins=50, ybins=50,
        xlim=None, ylim=None):
    """
    Plot a distribution of particles
    in the longitudinal phase space, and see the longitudinal profiles
    in time and energy.

    Parameters
    ----------
    time_coordinates : np.array or list
        The time coordinates of a distribution of particles in [s]
    energy_coordinates : np.array or list
        The energy coordinates of a distribution of particles in [eV]
    time_sep : np.array
        The time of the separatrix array in [s], as output by the separatrix function
    separatrix_array : np.array
        The separatrix array in [eV], as output by the separatrix function
    figname : str, optional
        The name of your nice figure, by default None
    xbins : int, optional
        The number of bins to generate a nice longitudinal profile in time, by default 50
    ybins : int, optional
        The number of bins to generate a nice longitudinal profile in energy, by default 50
    xlim : tuple, optional
        The limits in time for your nice plot in [s], e.g. (-t_rf/2, t_rf), by default None
    ylim : tuple, optional
        The limits in energy for your nice plot [eV], e.g. (-1, 1), by default None
    """

    plt.figure(figname, figsize=(8, 8))
    plt.clf()
    # Definitions for placing the axes
    left, width = 0.115, 0.63
    bottom, height = 0.115, 0.63
    bottom_h = left_h = left + width + 0.03

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    # rect_txtBox= [left_h, bottom_h, 0.2, 0.2]

    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    axScatter = plt.axes(rect_scatter)

    # txtBox = plt.axes(rect_txtBox)

    global distri_plot
    distri_plot, = axScatter.plot(
        time_coordinates * 1e9, energy_coordinates / 1e6, 'o', alpha=0.5)

    if (time_sep is not None) and (separatrix_array is not None):
        axScatter.plot(time_sep * 1e9, separatrix_array / 1e6, 'r')

    start_xlim = axScatter.get_xlim()
    start_ylim = axScatter.get_ylim()

    if xlim is not None:
        hist_dt = np.histogram(time_coordinates * 1e9, xbins, range=xlim)
    elif time_sep is not None:
        hist_dt = np.histogram(
            time_coordinates * 1e9, xbins,
            range=(np.nanmin(time_sep) * 1e9, np.nanmax(time_sep) * 1e9))
    else:
        hist_dt = np.histogram(time_coordinates * 1e9, xbins, range=xlim)

    global line_dt
    line_dt, = axHistx.plot(hist_dt[1][0:-1] + (
        hist_dt[1][1] - hist_dt[1][0]) / 2, hist_dt[0] / np.max(hist_dt[0]))
    axHistx.axes.get_xaxis().set_ticklabels([])
    axHistx.axes.get_yaxis().set_ticklabels([])
    axHistx.set_xlim(start_xlim)
    # axHistx.set_ylim(start_ylim)
    axHistx.set_ylabel('Bunch profile $\\lambda_{\\Delta t}$')

    if ylim is not None:
        hist_energy = np.histogram(energy_coordinates / 1e6, ybins, range=ylim)
    elif separatrix_array is not None:
        hist_energy = np.histogram(
            energy_coordinates / 1e6, ybins,
            range=(np.nanmin(separatrix_array) / 1e6, np.nanmax(separatrix_array) / 1e6))
    else:
        hist_energy = np.histogram(energy_coordinates / 1e6, ybins, range=ylim)

    global line_energy
    line_energy, = axHisty.plot(hist_energy[0] / np.max(
        hist_energy[0]), hist_energy[1][0:-1] + (hist_energy[1][1] - hist_energy[1][0]) / 2)
    axHisty.axes.get_xaxis().set_ticklabels([])
    axHisty.axes.get_yaxis().set_ticklabels([])
    axHisty.set_ylim(start_ylim)
    axHisty.set_xlabel('Energy spread $\\lambda_{\\Delta E}$')

    axScatter.set_xlabel('Time $\\Delta t$ [ns]')
    axScatter.set_ylabel('Energy $\\Delta E$ [MeV]')
    plt.xlim(xlim)
    plt.ylim(ylim)


def separatrix(time_array, f_rev, eta, beta, energy, charge, voltage, harmonic, acceleration=0):
    """Return the separatrix as an array for plotting purposes (together with the corresponding
    time values).

    Parameters
    ----------
    time_array : np.array
        The input time array in [s]
    f_rev : float
        The revolution frequency in [Hz]
    eta : float
        THe time slippage factor
    beta : float
        The relativistic beta
    energy : float
        The beam total energy in [eV]
    charge : float
        The particle charge in [e]
    voltage : float
        The rf voltage in [V]
    harmonic : float
        The rf harmonic number
    acceleration : float
        The beam energy gain per turn in [eV]

    Returns
    -------
    time_sep, separatrix_array : np.array
        The corresponding time and sepatrix values
    """

    warnings.filterwarnings("once")

    if eta > 0:
        phi_s = np.pi - np.arcsin(acceleration / charge / voltage)
        if acceleration > 0:
            phi_ufp = (np.pi - phi_s)
        else:
            phi_ufp = 2 * np.pi + (np.pi - phi_s)
    else:
        phi_s = np.arcsin(acceleration / charge / voltage)
        if acceleration > 0:
            phi_ufp = np.pi - phi_s
        else:
            phi_ufp = -np.pi - phi_s

    def pot_well(phi):
        return -(np.cos(phi) + phi * np.sin(phi_s)) / np.cos(phi_s)

    f_s0 = np.sqrt(
        -(2 * np.pi * f_rev)**2 * harmonic * eta * charge * voltage * np.cos(phi_s) /
        (2 * np.pi * beta**2 * energy)) / (2 * np.pi)

    sync_tune = f_s0 / f_rev

    sep_fac = sync_tune * beta**2 / (harmonic * np.abs(eta)) * energy
    sqrt_arg = 2 * (pot_well(phi_ufp) - pot_well(2 *
                    np.pi * harmonic * f_rev * time_array))
    sqrt_arg[sqrt_arg < 0] = 0
    separatrix_array = np.sqrt(sqrt_arg) * sep_fac
    separatrix_array[separatrix_array == 0] = np.nan

    separatrix_array = np.append(separatrix_array, -separatrix_array[::-1])
    time_sep = np.append(time_array, time_array[::-1])

    return time_sep[np.isfinite(separatrix_array)], separatrix_array[np.isfinite(separatrix_array)]


def generate_bunch(bunch_position, bunch_length,
                   bunch_energy, energy_spread,
                   n_macroparticles):
    """Generate a nice bunch of particles distributed as a parabola in phase space.

    Parameters
    ----------
    bunch_position : float
        The position in time [s] of the center of mass of the bunch
    bunch_length : float
        The length in time [s] of the bunch
    bunch_energy : float
        The position in energy [eV] of the center of mass of the bunch
        (relative to the synchronous energy)
    energy_spread : float
        The spread in energy [eV] of the bunch
    n_macroparticles : int
        The number of macroparticles to generate.

    Returns
    -------
    particle_dt : np.array
        The distribution of macroparticles in time (s)
    particle_energy : np.array
        The distribution of particles in energy (eV)
    """

    # Generating time and energy arrays
    time_array = np.linspace(bunch_position - bunch_length / 2,
                             bunch_position + bunch_length / 2,
                             100)
    energy_array = np.linspace(bunch_energy - energy_spread / 2,
                               bunch_energy + energy_spread / 2,
                               100)

    # Getting Hamiltonian on a grid
    dt_grid, deltaE_grid = np.meshgrid(
        time_array, energy_array)

    # Bin sizes
    bin_dt = time_array[1] - time_array[0]
    bin_energy = energy_array[1] - energy_array[0]

    # Density grid
    isodensity_lines = ((dt_grid - bunch_position) / bunch_length * 2)**2. + \
        ((deltaE_grid - bunch_energy) / energy_spread * 2)**2.
    density_grid = 1 - isodensity_lines**2.
    density_grid[density_grid < 0] = 0
    density_grid /= np.sum(density_grid)

    # Generating particles randomly inside the grid cells according to the
    # provided density_grid
    indexes = np.random.choice(np.arange(0, np.size(density_grid)),
                               n_macroparticles, p=density_grid.flatten())

    # Randomize particles inside each grid cell (uniform distribution)
    particle_dt = (np.ascontiguousarray(
        dt_grid.flatten()[indexes] + (np.random.rand(n_macroparticles) - 0.5) * bin_dt))
    particle_energy = (np.ascontiguousarray(
        deltaE_grid.flatten()[indexes] + (np.random.rand(n_macroparticles) - 0.5) * bin_energy))

    return particle_dt, particle_energy


def oscillation_spectrum(dt_track, cycle_time=None, fft_zero_padding=0):
    """Compute the spectrum of a particle time oscillation by applying an fft.

    Parameters
    ----------
    dt_track : np.array or list
        The time oscillations of a particle in s, over a few synchrotron periods.
    cycle_time : np.array or list
        The cycle time corresponding to dt_track. If passed, the spectrum is given
        with a frequency array expressed in [Hz]
    fft_zero_padding : int, optional
        The number of points for zero padding, to get a nice spectrum, by default 0

    Returns
    -------
    freq_array : np.array
        The frequency array of the spectrum
    fft_osc : np.array
        The amplitude of the time oscillation spectrum
    """
    n_turns = len(dt_track)

    if cycle_time is None:
        d = 1.0
    else:
        d = (cycle_time[1] - cycle_time[0])

    freq_array = np.fft.rfftfreq(n_turns + fft_zero_padding, d=d)
    fft_osc = np.abs(
        np.fft.rfft(
            dt_track - np.mean(dt_track),
            n_turns + fft_zero_padding) * 2 / (n_turns))

    return freq_array, fft_osc


def synchrotron_tune(dt_track, cycle_time=None, fft_zero_padding=0):
    """Compute the synchrotron tune from a particle time oscillations.

    Parameters
    ----------
    dt_track : np.array or list
        The time oscillations of a particle in s, over a few synchrotron periods.
    cycle_time : np.array or list
        The cycle time corresponding to dt_track. If passed, the spectrum is given
        with a frequency array expressed in [Hz]
    fft_zero_padding : int, optional
        The number of points for zero padding, to get a nice spectrum, by default 0

    Returns
    -------
    oscillation_amplitude : float
        The amplitude of the particle time oscillation in s
    sync_tune : float
        The synchrotron tune of the particle, or synchrotron frequency in [Hz]
        if cycle_time is passed.
    """
    freq_array, spectrum_array = oscillation_spectrum(
        dt_track, cycle_time=cycle_time, fft_zero_padding=fft_zero_padding)

    oscillation_amplitude = np.max(spectrum_array)
    sync_tune = float(
        np.mean(freq_array[spectrum_array == oscillation_amplitude]))

    return oscillation_amplitude, sync_tune


class _TrackAnimation(object):

    def __init__(
            self, time_coordinates, energy_coordinates,
            drift_function, rf_kick_function,
            drift_args, rf_kick_args,
            figname, iterations, framerate,
            xbins=50, ybins=50, xlim=None, ylim=None,
            time_sep=None, separatrix_array=None):

        self.time_coordinates = time_coordinates
        self.energy_coordinates = energy_coordinates
        self.drift_function = drift_function
        self.rf_kick_function = rf_kick_function
        self.drift_args = drift_args
        self.rf_kick_args = rf_kick_args
        self.figname = figname
        self.iterations = iterations
        self.framerate = framerate
        self.xbins = xbins
        self.ybins = ybins
        self.xlim = xlim
        self.ylim = ylim
        self.time_sep = time_sep
        self.separatrix_array = separatrix_array

    def run_animation(self):

        self._init()
        anim = animation.FuncAnimation(
            self.anim_fig, self._animate, init_func=self._init,
            frames=self.iterations, interval=1000 / self.framerate, blit=True)
        return HTML(anim.to_jshtml())

    def save_animation(self, filename):

        self._init()
        anim = animation.FuncAnimation(
            self.anim_fig, self._animate, init_func=self._init,
            frames=self.iterations, interval=1000 / self.framerate, blit=True)
        anim.save(filename, writer='ffmpeg')
        return filename

    def _init(self):

        plot_phase_space_distribution(
            self.time_coordinates,
            self.energy_coordinates,
            figname=self.figname,
            xbins=self.xbins, ybins=self.ybins,
            xlim=self.xlim, ylim=self.ylim,
            time_sep=self.time_sep,
            separatrix_array=self.separatrix_array)

        self.anim_fig = plt.gcf()

        return (line_dt, line_energy, distri_plot)

    def _animate(self, i):

        self.time_coordinates = self.drift_function(
            self.time_coordinates,
            self.energy_coordinates, *self.drift_args)
        self.energy_coordinates = self.rf_kick_function(
            self.energy_coordinates,
            self.time_coordinates, *self.rf_kick_args)

        if self.xlim is not None:
            hist_dt = np.histogram(
                self.time_coordinates * 1e9, self.xbins, range=self.xlim)
        elif self.time_sep is not None:
            hist_dt = np.histogram(
                self.time_coordinates * 1e9, self.xbins,
                range=(np.nanmin(self.time_sep) * 1e9, np.nanmax(self.time_sep) * 1e9))
        else:
            hist_dt = np.histogram(
                self.time_coordinates * 1e9, self.xbins, range=self.xlim)

        line_dt.set_data(
            hist_dt[1][0:-1] +
            (hist_dt[1][1] - hist_dt[1][0]) / 2,
            hist_dt[0] / np.max(hist_dt[0]))

        if self.ylim is not None:
            hist_energy = np.histogram(
                self.energy_coordinates / 1e6, self.ybins, range=self.ylim)
        elif self.separatrix_array is not None:
            hist_energy = np.histogram(
                self.energy_coordinates / 1e6, self.ybins,
                range=(np.nanmin(self.separatrix_array) / 1e6, np.nanmax(self.separatrix_array) / 1e6))
        else:
            hist_energy = np.histogram(
                self.energy_coordinates / 1e6, self.ybins, range=self.ylim)

        line_energy.set_data(
            hist_energy[0] / np.max(hist_energy[0]),
            hist_energy[1][0:-1] + (hist_energy[1][1] - hist_energy[1][0]) / 2)

        distri_plot.set_data(self.time_coordinates * 1e9,
                             self.energy_coordinates / 1e6)

        return (line_dt, line_energy, distri_plot)


def run_animation(time_coordinates,
                  energy_coordinates,
                  drift_function, rf_kick_function,
                  drift_args, rf_kick_args,
                  figname, iterations, framerate,
                  time_sep=None, separatrix_array=None,
                  xbins=50, ybins=50, xlim=None, ylim=None):
    """A routine to animate the motion of particles in the
    longitudinal phase space based on your own tracking equations.

    Parameters
    ----------
    dt_trajectory : np.array or list
        The time coordinates of a distribution of particles in [s]
    energy_trajectory : np.array or list
        The energy coordinates of a distribution of particles in [eV]
    drift_function : function
        Your drift equation of motion with the syntax
        time_coordinates = drift_function(time_coordinates, energy_coordinates, *drift_args)
        where drift_args is a list of arguments for the drift equation of motion (e.g. slippage)
    rf_kick_function : function
        Your kick equation of motion with the syntax
        energy_coordinates = kick_function(energy_coordinates, time_coordinates, *kick_args)
        where kick_args is a list of arguments for the drift equation of motion (e.g. rf voltage)
    figname : str
        The name of your nice animation
    iterations : int
        The number of iterations for the tracking (i.e. number of turns)
    framerate : float
        The framerate of the animation (e.g. 30 fps)
    time_sep : np.array
        The time of the separatrix array in [s], as output by the separatrix function
    separatrix_array : np.array
        The separatrix array in [eV], as output by the separatrix function
    xbins : int, optional
        The number of bins to generate a nice longitudinal profile in time, by default 50
    ybins : int, optional
        The number of bins to generate a nice longitudinal profile in energy, by default 50
    xlim : tuple, optional
        The limits in time for your nice plot in [s], e.g. (-t_rf/2, t_rf), by default None
    ylim : tuple, optional
        The limits in energy for your nice plot [eV], e.g. (-1e6, 1e6), by default None
    """
    trackanim = _TrackAnimation(
        time_coordinates,
        energy_coordinates,
        drift_function, rf_kick_function,
        drift_args, rf_kick_args,
        figname, iterations, framerate,
        time_sep=time_sep, separatrix_array=separatrix_array,
        xbins=xbins, ybins=ybins,
        xlim=xlim, ylim=ylim)

    return trackanim.run_animation()
