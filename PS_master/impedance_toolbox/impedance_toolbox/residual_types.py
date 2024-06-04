import numpy as np

def log_real_imag(source, fit):

    real_diff = (np.log(source[:,1]) - np.log(fit.impedance.real))**2
    imag_diff = (np.log(source[:,2]) - np.log(fit.impedance.imag))**2
    residue_real = np.sqrt(np.nansum(real_diff[np.isfinite(real_diff)]))
    residue_imag = np.sqrt(np.nansum(imag_diff[np.isfinite(imag_diff)]))
    residue = residue_real + residue_imag

    return residue

def log_total(source, fit):

    abs_diff = (np.log(np.abs(source[:,1] + source[:,2]*1j)) - np.log(np.abs(fit.impedance)))**2
    residue = np.sqrt(np.nansum(abs_diff[np.isfinite(abs_diff)]))

    return residue


def lin_real_imag(source, fit):

    residue_real = np.sqrt(np.sum((source[:,1] - fit.impedance.real)**2))
    residue_imag = np.sqrt(np.sum((source[:,2] - fit.impedance.imag)**2))
    residue = residue_real + residue_imag

    return residue


def lin_total(source, fit):

    residue = np.sqrt(np.nansum((np.abs(source[:,1] + source[:,2]*1j) - np.abs(fit.impedance))**2))

    return residue


def log_real_only(source, fit):

    real_diff = (np.log(source[:,1]) - np.log(fit.impedance.real))**2
    residue = np.sqrt(np.nansum(real_diff[np.isfinite(real_diff)]))

    return residue


def lin_real_only(source, fit):

    residue = np.sqrt(np.nansum((np.abs(source[:,1]) - np.abs(fit.impedance.real))**2))

    return residue

