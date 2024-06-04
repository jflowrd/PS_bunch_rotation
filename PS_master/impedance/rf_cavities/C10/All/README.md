# Impedance of 10 MHz cavities averaged

The impedance averaged over all 10 cavities (C10-11 excluded).

The returned impedance corresponds to 1x C10 cavity, with its
characteristics averaged over all the cavities. The averaging
is done by summing all the impedance spectra for all harmonic
configurations. The resulting impedance is then fitted with
resonators. Note that the impedance saved in files is divided
by 10 so that the loaded impedance accounts for a single cavity
(**files to be updated accordingly**).

## Subfolders and usage

- Resonators

The measured data was fitted using one or several resonators (9 presently) to
get the closest representation.

A fit with a single resonator is also avaible for simpler considerations. 

The fit is available for the same set of harmonics as measurements. The fit
should be loaded an multiplied by the number of cavities in a given harmonic.

- Parametric

The fit results are included in a parametrized model to extrapolate the impedance
for all harmonic settings of the cavity.

For harmonics that were not measured, the parameters of the impedance are interpolated
linearly.

## Plots

<details>
  <summary>Fit with resonators for 10 cavities summed</summary>
  <img src="Resonators/fitted_h8.png">
  <img src="Resonators/fitted_h16.png">
  <img src="Resonators/fitted_h21.png">
</details>

<details>
  <summary>Parametric model of averaged cavities</summary>
  <img src="Parametric/all_cavities.png">
</details>

<details>
  <summary>Parametric model for single cavities</summary>
  <img src="Parametric/Rs_single_cavities_summed.png">
  <img src="Parametric/Q_single_cavities_summed.png">
</details>

## History

To be completed

## Further work

- Impedance output to be corrected to provide only single cavity
impedance

## References

To be completed

## Credits

To be completed

