# Impedance of KFA45

The KFA45 is the injection kicker for the proton beam.

The first estimation was done using a Tsutsui model available on the following [PS impedance webpage](http://impedance.web.cern.ch/impedance/)
and on the [IRIS PS repository](https://gitlab.cern.ch/IRIS/PS_IW_model/-/tree/master/Impedances/Longitudinal).

More recent evaluation consisted in modelling the kicker
using CST. The model was built from old drawings available on
DFS (\\cern.ch\dfs\Workspaces\o\Old Drawings\Complexe_PS\PS\)
and CATIA files (reference needed).

The present CST model is stored on
[EDMS](https://edms.cern.ch/ui/#!master/navigator/document?P:100327819:100578961:subDocs).

Transition pieces are included in the CST model and their original purpose
was to damp low frequency resonances (see [KFA45](https://indico.cern.ch/event/830722/contributions/3479602/attachments/1870503/3077852/IWG_PS_KFA45_17_MeasSims_FINAL.pdf)).
The kicker is connected to HV ports. When ports are open, low frequency resonances are obtained
in wire measurements and are reproduced in CST wire simulations.
The coupling with HV ports is not included in the CST wakefield model.

<details>
  <summary>KFA45 in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss45.jpg">
</details>

## Subfolders and usage

- Tsutsui

The model was computed using ImpedanceWake2D with the following parameters
									
    - a=0.075 half width [m]
    - b=0.0265 half heigth (center-ferrite border) [m]
    - d=0.040 half heigth (center-PEC border including ferrite) [m]
    - L=0.884 Kicker's length [m]

- CST Raw Data

The wakefield is obtained with CST simulations. The model was
built based on 2D drawings.

The model include the vacuum flanges at either end of the magnet
for realstic boundary conditions.

The output result is the impedance in frequency domain and
can directly be used in particle simulations. The impedance
is noisy, to obtain a smoother response the fit with resonators
can be preferred.

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/cst_model_images/KFA45_CST_Model_View1.png">
  <img src="cst_raw_data/cst_model_images/KFA45_CST_Model_View2.png">
  <img src="cst_raw_data/cst_model_images/KFA45_CST_Model_View3.png">
  <img src="cst_raw_data/cst_model_images/KFA45_CST_Model_View4.png">
  <img src="cst_raw_data/cst_model_images/KFA45_CST_Model_View5.png">
  <img src="cst_raw_data/cst_model_images/KFA45_CST_Model_View6.png">
  <img src="cst_raw_data/cst_model_images/KFA45_CST_Model_View7.png">
</details>

- Resonators

The result of the CST simulation is fitted with resonators to get
the more global features of the impedance.

The broadband component of the impedance is fitted with 2 resonators.

## Plots

<details>
  <summary>Comparing Tsutsui model with CST</summary>
  <img src="Tsutsui/comparison_tsu_cst.png">
</details>

<details>
  <summary>Resonator fit</summary>
  <img src="Resonators/fitted_broadband.png">
  <img src="Resonators/fitted_broadband_realimag.png">
</details>

## History

To be completed

## Further work

To be completed

- Including coupling to HV ports in CST wakefield simulations
- Assess whether low frequency resonances are actually present and seen
by the beam in operational mode
- Confirm 8C11 material data

## References

To be completed

### Simulations

- [S. Persichelli, The beam coupling impedance model of CERN Proton Synchrotron](https://cds.cern.ch/record/2027523)

- [B. Popovic, 35th Impedance Working Group meeting](https://indico.cern.ch/event/844161/contributions/3562633/attachments/1911343/3158330/IWG_KFA45_17_with_cable_measurements_short_rev2.pdf)
- [B. Popovic, 33th Impedance Working Group meeting](https://indico.cern.ch/event/830722/contributions/3479602/attachments/1870503/3077852/IWG_PS_KFA45_17_MeasSims_FINAL.pdf)
- [B. Popovic, LIU-PS Beam Dynamics WG meeting #27](https://indico.cern.ch/event/801198/contributions/3329723/attachments/1803818/2943005/PS_Imped_Model_Update_280219_FINAL.pdf)
- [B. Popovic, 24th Impedance Working Group Meeting](https://indico.cern.ch/event/764129/contributions/3171804/attachments/1732002/2799672/KFA45_Ferrite_Study_240818_Rev3.pdf)
- [D. Ventura, 24th Impedance Working Group Meeting](https://indico.cern.ch/event/764129/contributions/3171804/attachments/1732002/2799708/KFA45_newFerrite.pdf)
- [B. Popovic, Longitudinal limitations with LIU-PS RF upgrades and mitigation strategy](https://indico.cern.ch/event/750790/contributions/3108016/attachments/1719965/2776247/RG_PS_Impedance_Model_Meeting.pdf)
- [B. Popovic, SPS Injection Losses Review](https://indico.cern.ch/event/672967/contributions/2753563/attachments/1567554/2471183/PS_Impedance_Status_Injection_Losses_Meeting_301117.pdf)
- [B. Popovic, LIU-PS Beam Dynamics WG meeting \#1](https://indico.cern.ch/event/662292/contributions/2704441/attachments/1516648/2367108/BP_PS_Impedance_Model_Update_310817_Final.pdf)

## Credits

To be completed

- Tsutsui: BE-ABP-HSC Serena Persichelli, Elias Metral, Nicolo Biancacci or Benoit Salvant
- Model: B. Popovic
- 8C11 Material: S. Persichelli
