# Impedance of 10 MHz cavities with closed gap relay

A mechanical gap relay is installed on each gap of the cavity to short
circuit the gap and reduce the beam coupling impedance when not in use.
Although greatly reduced, the impedance is not completely negligible.

The impedance was obtained from beam-based measurements on a single gap
for different tuning configuration. The impedance of both gaps
are assumed identical (when both gaps are closed) and the output impedance
corresponds to one full cavity (2 gaps).

The impedance is fitted with two resonators with extra contribution ImZ/n
and a parametric model is built to extrapolate the impedance for all tuning
parameters.

## Subfolders and usage

- Measurements

The beam-based measurements of the impedance for different harmonics.

These are used as input for the fit and it is not advised to use them directly
in particle simulations.

- Resonators

The beam-based measurements are fitted with two resonators and extra ImZ/n.

All contributions should be used to represent well the impedance (resonators
are allowed to have negative shunt impedance to have a better result for the fit).
The fits are done for all the measured cases.

- Parametric

The fit results are included in a parametrized model to extrapolate the impedance
for all harmonic settings of the cavity.

For harmonics that were not measured, the parameters of the impedance are interpolated
linearly.

## Plots

<details>
  <summary>Fit with resonators and ImZ/f h=8</summary>
  <img src="Resonators/fitted_h8.png">
  <img src="Resonators/fitted_h8_real.png">
  <img src="Resonators/fitted_h8_imag.png">
</details>

<details>
  <summary>Fit with resonators and ImZ/f h=16</summary>
  <img src="Resonators/fitted_h16.png">
  <img src="Resonators/fitted_h16_real.png">
  <img src="Resonators/fitted_h16_imag.png">
</details>

<details>
  <summary>Fit with resonators and ImZ/f h=21</summary>
  <img src="Resonators/fitted_h21.png">
  <img src="Resonators/fitted_h21_real.png">
  <img src="Resonators/fitted_h21_imag.png">
</details>

<details>
  <summary>Parametric model with cavities summed</summary>
  <img src="Parametric/all_cavities_shunt.png">
  <img src="Parametric/all_cavities_Q.png">
  <img src="Parametric/all_cavities_fr.png">
  <img src="Parametric/all_cavities_ImZ_over_f.png">
</details>

## History

To be completed

- Impedance with gap relays measured and modelled by G. Favia
- Impedance identified as the driving source for quadrupolar coupled-bunch
instability (ref. needed)

## Further work

To be completed

## References

To be completed

### Measurements
- [G. Favia, Study of the beam-cavity interaction in the CERN PS 10 MHz 
cavities and investigation of hardware solutions to reduce beam loading, 
PhD Thesis](https://cds.cern.ch/record/2286835)

### Other
- [G. Favia et al., Study of the beam-cavity interaction in the CERN PS 10 MHz 
rf system](https://cds.cern.ch/record/2207324/files/mopor012.pdf)

## Credits

To be completed

- Measured data: G. Favia
- Resonator fit: A. Lasheen

