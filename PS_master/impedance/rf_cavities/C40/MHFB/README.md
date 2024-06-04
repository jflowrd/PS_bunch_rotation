# Impedance of 40 MHz cavities with multi harmonic feedback

The impedance of the cavities with the multi harmonic feedback
(MHFB) enabled.

The model is based on the fit of the measured impedance (see [here](../Individual/README.md)),
which is filtered with independant filters on each revolution harmonic.

The impedance with MHFB is obtained using is applied using

```math
Z_\mathrm{fb}\left(f\right)=\frac{Z_c\left(f\right)}{1+H\left(f\right)Z_c\left(f\right)}
```

where $`Z_c`$ is the cavity impedance and $`H=H_{\mathrm{mhfb}}`$
the transfer function of the MHFB filter which is a
sum of filters at revolution harmonics defined as

```math
H_{\mathrm{mhfb}}\left(f\right)=g\,\sum_{m}{\frac{1}{1 + j Q \left(\frac{\omega}{\omega_{r,m}}-\frac{\omega_{r,m}}{\omega}\right)}e^{-j \omega \tau_m}}
```

where $`g`$ is a gain parameter and $`\tau_m\approx t_\mathrm{rev}`$ is the
delay of the filter adjusted individually for each notch to maximize the impedance reduction
at a given revolution harmonic. The quality factor $`Q`$ 
is adjusted to achieve the requested bandwidth of the notches. 

When the cavity is pulsing, the filter at the central frequency is disabled.
This model uses an extra resonator to flatten the impedance response at the central
frequency (perturbed by the neighboring notches) with
$`H=H_{\mathrm{mhfb}}+H_n`$ where

```math
H_n\left(f\right)=g_n\,\frac{1}{1 + j Q \left(\frac{\omega}{\omega_r}-\frac{\omega_r}{\omega}\right)}e^{-j \omega \tau_n}
```

$`g_n`$ is a gain parameter adjusted to compensate
for the comb filter at the central frequency of the cavity $`\omega_r=h\omega_\mathrm{rev}`$ 
and $`\tau_n`$ is the delay of the notch filter. The quality factor $`Q`$ 
is adjusted to have the flatest transfer function around the central frequency of the cavity. 

The resonator filters can be replaced by IIR filters. This model was built
to reproduce the measured transfer function but a better model would 
consist in using the same filter definition as implemented in the FPGA.

In total, the impedance at 11 revolution harmonics are reduced (the central
harmonic, 5 below, 5 above).

## Usage

- Resonators

Based on C40 Resonator model. The model consists of a modified impedance table.
The impedance at h=84 can be reduced or not (depending if the cavity is active or not).
The model for the central harmonic compensation is only valid when 84*f_rev is close to the
cavity central frequency

## Plots

<details>
  <summary>Feedback with -20dB reduction, 5kHz bandwidth, with the central harmonic damped</summary>
  <img src="Resonators/C40-77_impedance_dB_with_MHFB_-20dB_477kHz_central_False.png">
  <img src="Resonators/C40-77_impedance_with_MHFB_-20dB_477kHz_central_False.png">
</details>

<details>
  <summary>Feedback with -20dB reduction, 5kHz bandwidth, with the central harmonic undamped when pulsing</summary>
  <img src="Resonators/C40-77_impedance_dB_with_MHFB_-20dB_477kHz_central_True.png">
  <img src="Resonators/C40-77_impedance_with_MHFB_-20dB_477kHz_central_True.png">
</details>


## History

To be completed

## Further work

To be completed

- Use more accurate definition of the MHFB filters as implemented in the machine
- Use the DFB model instead of the Resonator as base before filtering

## References

### MHFB
- [F. Bertin et al., Impedance reduction of the High-frequency Cavities in the CERN PS by Multi-harmonic Feedback](https://cds.cern.ch/record/2699578/)

## Credits

To be completed
