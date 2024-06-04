# Impedance of 80 MHz cavities fundamental mode

In normal operation, two of the three cavities are tuned for h=168 for proton beam at flat top (C80-08, C80-88). 
The third cavity (C80-89) is tuned for h=169 for ion operation (values present for Pb54+ and Xe39+). The characteristics of
the cavity are (see reference):

| Parameter \ Frequency [MHz] 			| 80		|
| -------------------------------------	|:---------:|
| Maximum VRF [kV]       				| 300		|
| Shunt resistance R [kOhm]				| 1260 		|
| Quality factor Q 						| 22600 	|
| R/Q [Ohm] 							| 56		|

When not in use, the cavity gap can be closed mecanically. Presently both proton
cavities are open during operation while the ion one is only open during the ion run.

The impedance model in this database includes the effect of the fast feedback 
loop:

| Parameter \ Frequency [MHz] 			| 80		|
| -------------------------------------	|:---------:|
| Reduced quality factor Q				| 100 		|

Measurements with beam were performed. The measured shunt impedance is higher
than expected and needs to be further investigated (see references, measurements 
were performed before the release of the AVC loop and after the warning timing).

## Subfolders and usage

- Measurements

Using the resonator fit of the measured impedance with beam. Measurements were 
done for different number of bunches, take the relevant input file depending
on your scenario.

Measurements were performed before the release of the AVC loop, and after the
warning timing

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
  <summary>C80-08 Resonator fit (measured with 12 bunches to 72 bunches in steps of 12 bunches)</summary>
  <img src="C08/Resonators/fitted_b12.png">
  <img src="C08/Resonators/fitted_b24.png">
  <img src="C08/Resonators/fitted_b36.png">
  <img src="C08/Resonators/fitted_b48.png">
  <img src="C08/Resonators/fitted_b60.png">
  <img src="C08/Resonators/fitted_b72.png">
</details>

<details>
  <summary>C80-08 Resonator and DFB fit (measured with 12 bunches to 72 bunches in steps of 12 bunches)</summary>
  <img src="C08/DFB/fitted_b12.png">
  <img src="C08/DFB/fitted_b24.png">
  <img src="C08/DFB/fitted_b36.png">
  <img src="C08/DFB/fitted_b48.png">
  <img src="C08/DFB/fitted_b60.png">
  <img src="C08/DFB/fitted_b72.png">
</details>

<details>
  <summary>C80-08 Resonator and DFB fit compared to design values</summary>
  <img src="C08/DFB/final_rshunt.png">
  <img src="C08/DFB/final_Q.png">
  <img src="C08/DFB/final_RoQ.png">
  <img src="C08/DFB/final_gain.png">
  <img src="C08/DFB/final_delay.png">
</details>

<details>
  <summary>C80-88 Resonator fit (measured with 12 bunches to 72 bunches in steps of 12 bunches)</summary>
  <img src="C88/Resonators/fitted_b12.png">
  <img src="C88/Resonators/fitted_b24.png">
  <img src="C88/Resonators/fitted_b36.png">
  <img src="C88/Resonators/fitted_b48.png">
  <img src="C88/Resonators/fitted_b60.png">
  <img src="C88/Resonators/fitted_b72.png">
</details>

<details>
  <summary>C80-88 Resonator and DFB fit (measured with 12 bunches to 72 bunches in steps of 12 bunches)</summary>
  <img src="C88/DFB/fitted_b12.png">
  <img src="C88/DFB/fitted_b24.png">
  <img src="C88/DFB/fitted_b36.png">
  <img src="C88/DFB/fitted_b48.png">
  <img src="C88/DFB/fitted_b60.png">
  <img src="C88/DFB/fitted_b72.png">
</details>

<details>
  <summary>C80-88 Resonator and DFB fit compared to design values</summary>
  <img src="C88/DFB/final_rshunt.png">
  <img src="C88/DFB/final_Q.png">
  <img src="C88/DFB/final_RoQ.png">
  <img src="C88/DFB/final_gain.png">
  <img src="C88/DFB/final_delay.png">
</details>


<details>
  <summary>C80-89 Resonator fit (measured with 12 bunches to 72 bunches in steps of 12 bunches)</summary>
  <img src="C89/Resonators/fitted_b12.png">
  <img src="C89/Resonators/fitted_b24.png">
  <img src="C89/Resonators/fitted_b36.png">
  <img src="C89/Resonators/fitted_b48.png">
  <img src="C89/Resonators/fitted_b60.png">
  <img src="C89/Resonators/fitted_b72.png">
</details>

<details>
  <summary>C80-89 Resonator and DFB fit (measured with 12 bunches to 72 bunches in steps of 12 bunches)</summary>
  <img src="C89/DFB/fitted_b12.png">
  <img src="C89/DFB/fitted_b24.png">
  <img src="C89/DFB/fitted_b36.png">
  <img src="C89/DFB/fitted_b48.png">
  <img src="C89/DFB/fitted_b60.png">
  <img src="C89/DFB/fitted_b72.png">
</details>

<details>
  <summary>C80-89 Resonator and DFB fit compared to design values</summary>
  <img src="C89/DFB/final_rshunt.png">
  <img src="C89/DFB/final_Q.png">
  <img src="C89/DFB/final_RoQ.png">
  <img src="C89/DFB/final_gain.png">
  <img src="C89/DFB/final_delay.png">
</details>

## Further work

- Understand the discrepancy between the expectation in terms of Rshunt with
respect to the measured one.
- Include the impedance with mechanical short
- Other tasks to be determined (continue beam based measurement of the impedance ? 
modellilng of the HOMs with couplers ? impedance for low vs. high voltage in the gap ? 
simulation based modelling ?)

## References

### Design
- [The PS complex as proton pre-injector for the LHC : design and implementation report](https://cds.cern.ch/record/449242)
- [E. Jensen et al., The PS 80 MHz cavities](https://accelconf.web.cern.ch/accelconf/e98/PAPERS/TUP02H.PDF)

### Beam based measurements
- [A. Lasheen, G. Favia, Measurements of cavities impedance with beam (40MHz and 80MHz)](https://indico.cern.ch/event/685455/#4-measurements-of-cavities-imp)

### Other
- [Fifty years of the CERN Proton Synchrotron : Volume 1](https://cds.cern.ch/record/1359959)

## Credits

To be completed
