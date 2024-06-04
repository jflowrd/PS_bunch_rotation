# Impedance of Internal Dumps

There are 2 internal dumps in the PS (SD 47 & SD48).

The internal dump was modeled using CST and is stored on
[EDMS](https://edms.cern.ch/ui/#!master/navigator/project?P:100327819:100512435:subDocs).

Both internal dumps are replaced during LS2 by a new design.

<details>
  <summary>Internal dumps in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss47.jpg">
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss48.jpg">
</details>

## Usage

- CST Raw Data

The wakefield and eigenmodes are obtained using CST.

The PreLS2 model was built from ST0490320_01.
The eigenmode and wakefield give comparable results. However, a resonance
at higher frequency than the maximum frequency considered in simulations
seem to contribute with broandband inductive impedance.

<details>
  <summary>PreLS2 CST model view</summary>
  <img src="PreLS2/cst_raw_data/cst_model_images/internal_beam_dump_Model_View1.png">
  <img src="PreLS2/cst_raw_data/cst_model_images/internal_beam_dump_Model_View2.png">
  <img src="PreLS2/cst_raw_data/cst_model_images/internal_beam_dump_Model_View3.png">
  <img src="PreLS2/cst_raw_data/cst_model_images/internal_beam_dump_Model_View4.png">
  <img src="PreLS2/cst_raw_data/cst_model_images/internal_beam_dump_Model_View5.png">
</details>

PostLS2 model info to be completed.

<details>
  <summary>PostLS2 CST model view</summary>
  <img src="PostLS2/cst_raw_data/cst_model_images/NewDump_CST_Model_View1.png">
  <img src="PostLS2/cst_raw_data/cst_model_images/NewDump_CST_Model_View2.png">
  <img src="PostLS2/cst_raw_data/cst_model_images/NewDump_CST_Model_View3.png">
</details>

- Resonators

The result of the CST simulation is fitted with resonators to get
the more global features of the impedance.

For the PreLS2 model, the impedance is fitted with 8 resonators and an extra contribution of 
constance ImZ/f.
The eigenmode being in excellent agreement with the wakefield simulation
and including higher frequency components, the eigenmode is recommended over
the resonator fit.

For the PostLS2 model, many resonances are included but makes the convergence
of the fit more difficult. *The eigenmode calculations from CST could
allow to avoid this issue.*

## Plots

<details>
  <summary>Comparison eigen/wake PreLS2</summary>
  <img src="PreLS2/cst_raw_data/plot_comparison/cst_comparison.png">
</details>

<details>
  <summary>Resonator fit PreLS2</summary>
  <img src="PreLS2/Resonators/fitted_impedance_wake.png">
</details>

<details>
  <summary>Resonator fit PostLS2</summary>
  <img src="PreLS2/Resonators/fitted_impedance_wake.png">
</details>

## Further work

To be completed

## References

To be completed

- [B. Popovic, Impedance Working Group meeting #14](https://indico.cern.ch/event/671318/contributions/2746286/attachments/1538182/2410894/IWG14_Preliminary_PS_Dump_Simulations_091017_FINAL.pdf)

## Credits

- CST Model: B. Popovic
