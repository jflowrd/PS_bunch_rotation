# Impedance of 10 MHz cavities with 1-turn delay feedback

The impedance of the cavities with the 1-turn delay feedback
(OTFB) enabled.

The model is based on the parametric model of the averaged
C10 cavity (see [here](../All/README.md)), which is filtered with a comb-filter.

The impedance with OTFB is obtained using is applied using

```math
Z_\mathrm{fb}\left(f\right)=\frac{Z_c\left(f\right)}{1+H\left(f\right)Z_c\left(f\right)}
```

where $`Z_c`$ is the cavity impedance and $`H=H_c`$
the transfer function of the OTFB filter which is a
comb-filter defined as

```math
H_c\left(f\right)=g_c\,\frac{a}{1-\left(1-a\right)e^{-j \omega \tau_c}}e^{-j \omega \tau_\mathrm{otfb}}
```

where $`a=2^{-4}`$, $`g_c`$ is a gain parameter and $`\tau_c=t_\mathrm{rev}`$ is the
periodicity of the comb filter, to get beam loading signal at mutiples of $`f_\mathrm{rev}`$, 
and $`\tau_\mathrm{otfb}=t_\mathrm{rev}`$ delay corresponding to one revolution period.

When the cavity is pulsing, the signal at the central frequency
of the cavity is rejected by a notch filter and $`H=H_c+H_n`$.
The notch filter is presently simplified as a resonator with 

```math
H_n\left(f\right)=g_n\,\frac{1}{1 + j Q \left(\frac{\omega}{\omega_r}-\frac{\omega_r}{\omega}\right)}e^{-j \omega \tau_n}
```

where $`g_n`$ is a gain parameter adjusted to compensate
for the comb filter at the central frequency of the cavity $`\omega_r=h\omega_\mathrm{rev}`$ 
and $`\tau_n=t_\mathrm{rev}`$ is the delay of the notch filter. The quality factor $`Q`$ 
is adjusted to have the flatest transfer function around the central frequency of the cavity. 

## Subfolders and usage

- Parametric

The model is based on the parametric model of the averaged
C10 impedance, on which extra parameters were included to
model the OTFB.

The user specifies the cavity harmonic, the revolution period, 
the impedance reduction target and the bandwidth at the first revolution harmonics
around the central frequency. The relevant feedback parameters
(gain, and central notch Q factor) are interpolated from a lookup table.

## Plots

<details>
  <summary>Feedback at h=7 and -15dB impedance reduction on first revolution harmonic sidebands</summary>
  <img src="Parametric/impedance_dB_with_1TFB_h7_-15dB.png">
  <img src="Parametric/impedance_with_1TFB_h7_-15dB.png">
</details>

<details>
  <summary>Feedback at h=21 and -15dB impedance reduction on first revolution harmonic sidebands</summary>
  <img src="Parametric/impedance_dB_with_1TFB_h21_-15dB.png">
  <img src="Parametric/impedance_with_1TFB_h21_-15dB.png">
</details>

## History

To be completed

## Further work

- Central notch is modelled by resonator, to be replaced
by actual notch implementation
- Use the actual implementation like in the FPGA
- Compare transfer function with measurements

## References

To be completed

- [A. Blas, R. Garoby, Design and operational results of a "one-turn-delay feedback"
for beam loading compensation of the CERN PS ferrite cavities](https://cds.cern.ch/record/220107/files/CM-P00059458.pdf)
- [D. Perrelet, New 1-turn feedback for the PS 10 MHz cavities](https://indico.cern.ch/event/254244/contributions/1581263/attachments/444439/616448/Presentation_New_1TFB_LIU_MEETING_DP.pdf)
- [D. Perrelet, New PS one-turn delay feedbacks and further developments](https://indico.cern.ch/event/299470/contributions/686507/attachments/564148/777098/Presentation_1TFB_AVC_TFB_CBFB_LIU-DAY.pdf)

## Credits

To be completed

