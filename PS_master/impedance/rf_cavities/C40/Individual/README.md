# Impedance of 40 MHz cavities fundamental mode

The cavities are tuned for h=84 for proton beam at flat top. The characteristics of
the cavity are (see reference):

| Parameter \ Frequency [MHz] 			| 40		|
| -------------------------------------	|:---------:|
| Maximum VRF [kV]       				| 300		|
| Shunt resistance R [kOhm]				| 600 		|
| Quality factor Q 						| 18200 	|
| Loaded Shunt resistance R [kOhm]		| 330 		|
| Loaded Quality factor Q 				| 10000 	|
| R/Q [Ohm] 							| 33		|

When not in use, the cavity gap can be closed mecanically. Presently both
cavities are usually open during operation.

The impedance model in this database includes the effect of the fast feedback 
loop:

| Parameter \ Frequency [MHz] 			| 40		|
| -------------------------------------	|:---------:|
| Reduced quality factor Q				| 70 		|

Measurements with beam were performed. The measured shunt impedance is higher
than expected and needs to be further investigated (see references,
WARNING: measurements were performed after the release of the AVC loop, therefore
some rf leakage coming from the low-level rf is present in the measured impedance).

## Subfolders and usage

- Measurements

The measurements were performed with multibunch beam and only the amplitude
of the impedance at the revolution harmonics are sampled. The usage of a fitted
model is prefered

- Resonators

The measured impedance is fitted with a single resonator. The shunt impedance
is well fitted but the impedance is broader than the measured one. A more refined
model including the impedance reduction by a fast Direct Feedback can be used.

- DFB

The impedance is fitted by a resonator with a Direct Feedback. The 
impedance is 

```math
Z_\mathrm{fb}\left(f\right)=\frac{Z_c\left(f\right)}{1+H\left(f\right)Z_c\left(f\right)}
```

where $`Z_c`$ is the cavity impedance and $`H`$
the transfer function of the feedback

```math
H\left(f\right)=g\,e^{-j \omega \tau}
```

wher $`g`$ is the feedback gain and $`\tau=n T_{\mathrm{rf}}`$ is a delay
which is an integer multiple of the rf period $`T_{\mathrm{rf}}`$.

The loop gain is 

```math
G = 20\log{\left(R_s g\right)}
```

## Plots

<details>
  <summary>C40-77 Resonator fit (measured with 12 bunches to 72 bunches in steps of 12 bunches)</summary>
  <img src="C77/Resonators/fitted_b12.png">
  <img src="C77/Resonators/fitted_b24.png">
  <img src="C77/Resonators/fitted_b36.png">
  <img src="C77/Resonators/fitted_b48.png">
  <img src="C77/Resonators/fitted_b60.png">
  <img src="C77/Resonators/fitted_b72.png">
</details>

<details>
  <summary>C40-77 Resonator and DFB fit (measured with 12 bunches to 72 bunches in steps of 12 bunches)</summary>
  <img src="C77/DFB/fitted_b12.png">
  <img src="C77/DFB/fitted_b24.png">
  <img src="C77/DFB/fitted_b36.png">
  <img src="C77/DFB/fitted_b48.png">
  <img src="C77/DFB/fitted_b60.png">
  <img src="C77/DFB/fitted_b72.png">
</details>

<details>
  <summary>C40-77 Resonator and DFB fit compared to design values</summary>
  <img src="C77/DFB/final_rshunt.png">
  <img src="C77/DFB/final_Q.png">
  <img src="C77/DFB/final_RoQ.png">
  <img src="C77/DFB/final_gain.png">
  <img src="C77/DFB/final_delay.png">
</details>

<details>
  <summary>C40-78 Resonator fit (measured with 12 bunches to 72 bunches in steps of 12 bunches)</summary>
  <img src="C78/Resonators/fitted_b12.png">
  <img src="C78/Resonators/fitted_b24.png">
  <img src="C78/Resonators/fitted_b36.png">
  <img src="C78/Resonators/fitted_b48.png">
  <img src="C78/Resonators/fitted_b60.png">
  <img src="C78/Resonators/fitted_b72.png">
</details>

<details>
  <summary>C40-78 Resonator and DFB fit (measured with 12 bunches to 72 bunches in steps of 12 bunches)</summary>
  <img src="C78/DFB/fitted_b12.png">
  <img src="C78/DFB/fitted_b24.png">
  <img src="C78/DFB/fitted_b36.png">
  <img src="C78/DFB/fitted_b48.png">
  <img src="C78/DFB/fitted_b60.png">
  <img src="C78/DFB/fitted_b72.png">
</details>

<details>
  <summary>C40-78 Resonator and DFB fit compared to design values</summary>
  <img src="C78/DFB/final_rshunt.png">
  <img src="C78/DFB/final_Q.png">
  <img src="C78/DFB/final_RoQ.png">
  <img src="C78/DFB/final_gain.png">
  <img src="C78/DFB/final_delay.png">
</details>

## Further work

- Understand the discrepancy between the expectation in terms of Rshunt with
respect to the measured one.
- Other tasks to be determined (continue beam based measurement of the impedance ? 
modelling of the HOMs with couplers ? impedance for low vs. high voltage in the gap ? 
simulation based modelling ?)

## History

To be completed

## References

### Design
- [The PS complex as proton pre-injector for the LHC : design and implementation report](https://cds.cern.ch/record/449242)

### Beam based measurements
- [A. Lasheen, G. Favia, Measurements of cavities impedance with beam (40MHz and 80MHz)](https://indico.cern.ch/event/685455/#4-measurements-of-cavities-imp)

### Other
- [Fifty years of the CERN Proton Synchrotron : Volume 1](https://cds.cern.ch/record/1359959)

## Credits

To be completed
