# Impedance of individual 10 MHz cavities

The impedance of individual cavities. The impedance was obtained from
beam-based measurements and already include the contributions of the two
gaps.

The impedance was fitted with resonators for all cavities and all harmonics.
The DC components are removed from the measured data and set to 0 before the
fit.

The C10-11 impedance includes the contribution of the upgraded RF amplifier
that will be installed on all cavities during Long Shutdown 2 (2019-2020).

## Subfolders and usage

- Measurements

The impedance source can be loaded directly from the measured data (impedance
table).

This model may not be the most suitable for particle simulations since the
impedance data is cut at 25 Mhz and will include unphysical wake if the binning
of the bunch profile is too fine.

- Resonators

The measured data was fitted using one or several resonators (9 presently) to
get the closest representation.

A fit with a single resonator is also avaible for simpler considerations. 

The fit is available for the same set of harmonics as measurements. The fit should
be loaded for each cavity individually.

## Plots

<details>
  <summary>Single resonator fit</summary>
  <img src="Rs_single_cavities_overview.png">
  <img src="Q_single_cavities_overview.png">
</details>

<details>
  <summary>C10-11</summary>
  <img src="C10-11/Resonators/fitted_h8.png">
  <img src="C10-11/Resonators/fitted_h16.png">
  <img src="C10-11/Resonators/fitted_h21.png">
</details>

<details>
  <summary>C10-36</summary>
  <img src="C10-36/Resonators/fitted_h8.png">
  <img src="C10-36/Resonators/fitted_h16.png">
  <img src="C10-36/Resonators/fitted_h21.png">
</details>

<details>
  <summary>C10-46</summary>
  <img src="C10-46/Resonators/fitted_h8.png">
  <img src="C10-46/Resonators/fitted_h16.png">
  <img src="C10-46/Resonators/fitted_h21.png">
</details>

<details>
  <summary>C10-51</summary>
  <img src="C10-51/Resonators/fitted_h8.png">
  <img src="C10-51/Resonators/fitted_h16.png">
  <img src="C10-51/Resonators/fitted_h21.png">
</details>

<details>
  <summary>C10-56</summary>
  <img src="C10-56/Resonators/fitted_h8.png">
  <img src="C10-56/Resonators/fitted_h16.png">
  <img src="C10-56/Resonators/fitted_h21.png">
</details>

<details>
  <summary>C10-66</summary>
  <img src="C10-66/Resonators/fitted_h8.png">
  <img src="C10-66/Resonators/fitted_h16.png">
  <img src="C10-66/Resonators/fitted_h21.png">
</details>

<details>
  <summary>C10-76</summary>
  <img src="C10-76/Resonators/fitted_h8.png">
  <img src="C10-76/Resonators/fitted_h16.png">
  <img src="C10-76/Resonators/fitted_h21.png">
</details>

<details>
  <summary>C10-81</summary>
  <img src="C10-81/Resonators/fitted_h8.png">
  <img src="C10-81/Resonators/fitted_h16.png">
  <img src="C10-81/Resonators/fitted_h21.png">
</details>

<details>
  <summary>C10-86</summary>
  <img src="C10-86/Resonators/fitted_h8.png">
  <img src="C10-86/Resonators/fitted_h16.png">
  <img src="C10-86/Resonators/fitted_h21.png">
</details>

<details>
  <summary>C10-91</summary>
  <img src="C10-91/Resonators/fitted_h8.png">
  <img src="C10-91/Resonators/fitted_h16.png">
  <img src="C10-91/Resonators/fitted_h21.png">
</details>

<details>
  <summary>C10-96</summary>
  <img src="C10-96/Resonators/fitted_h8.png">
  <img src="C10-96/Resonators/fitted_h16.png">
  <img src="C10-96/Resonators/fitted_h21.png">
</details>


## History

To be completed

## Further work

- Simple parametric model for individual cavities.
- Fitting with resonators from the open loop impedance with a model
for the direct feedback loop

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

