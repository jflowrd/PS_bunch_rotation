# Impedance of Wall Current Monitors (WCM)

Two WCM are located in straight section 03.

The model is from the impedance website \[1\] and is shown in \[2\]. It was re-run for longitudinal wakefield simulation.

\[2\] Heavily documents its use and shows its ferrite data.

The model does not include the neighboring vacuum flanges.

**The present CST model provide a large shunt impedance and needs to be checked.** 

<details>
  <summary>WCM03 in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss03.jpg">
</details>

## Subfolders and usage

- Circuit Model

The WCM is modeled using a simplified parallel RLC circuit model with parameters from \[4\].

    - R = 6  # Ohm
    - C = 0.65e-11  # Farad
    - L = 10e-6  # Henri

The gap capacitance value was adjusted to have a higher cut-off frequency of ~4GHz (-3dB).

- CST Raw Data

The wakefield is obtained with CST simulations.

The output result is the impedance in frequency domain and
can directly be used in particle simulations. The impedance
is cut at high frequency, which can lead to unphysical wake.
For very short bunches the resonator fit can be preferred.

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/cst_model_images/WCM_CST_Model_View1.png">
</details>

<details>
	<summary>Warning</summary>
	<p>When running the model in CST Eigenmode you must deselect 'Consider losses in post processing only' button.
	Then you need to go into "Specials..." menu and select the "Materials" => "Evaluation Frequency: " and enter the frequency at which to evaluate the ferrites at (i.e. what mu value to use).</p> 
	<p>For this device it is absolutely necessary to know this frequency (~10 MHz), otherwise your resonance frequency and r-shunt will be wrong.
	The problem is to know what this frequency will be, since it is essentially at which frequency the main mode (at the gap) is resonating at, which you would try to get from eigenmode simulations.</p> 
	<p>I (Branko Popovic) used wakefield results to know that the main mode is at 10 MHz.
	I don't know a better way to do this without running an wakefield simulation where you can identify the resonance mode you are looking for.</p> 
</details>

- Resonators

The result of the CST simulation is fitted with resonators to get
the more global features of the impedance.

The broadband component of the impedance is fitted with 3 resonators
while resonances at low frequency are fitted with 4 resonators.

## Plots

<details>
  <summary>Circuit model</summary>
  <img src="circuit_model/abs_impedance.png">
</details>

<details>
  <summary>Resonator fit of CST model</summary>
  <img src="Resonators/fitted_impedance.png">
  <img src="Resonators/fitted_impedance_zoom.png">
  <img src="Resonators/fitted_impedance_realimag.png">
</details>

## History

To be completed

## Further work

To be completed

- Shunt impedance is large in CST model, to be checked
- Eigenmode data?

## References

To be completed

- \[1\] http://impedance.web.cern.ch/impedance/PS.htm
- \[2\] http://jeroen.web.cern.ch/jeroen/WCMIV/
- \[3\] ECR: 1290662 https://edms.cern.ch/ui/#!master/navigator/document?D:1887109572:1887109572:subDocs
- \[4\] [J. Belleman, A New Wall Current Monitor for the CERN Proton Synchrotron](https://cds.cern.ch/record/2313362/files/mopg41.pdf)
- [BI/TB Meetings 21/05/2015, Status and plans for PS WCM (Satellite and Ghost)](https://indico.cern.ch/event/390422/contributions/926910/)

## Credits

To be completed
