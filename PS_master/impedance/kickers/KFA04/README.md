# Impedance of KFA04

The KFA04 is part of the Multi-Turn Extraction system (MTE).

The first estimation was done using a Tsutsui model available on the following [PS impedance webpage](http://impedance.web.cern.ch/impedance/)
and on the [IRIS PS repository](https://gitlab.cern.ch/IRIS/PS_IW_model/-/tree/master/Impedances/Longitudinal).

More recent evaluation consisted in modelling the kicker
using CST. The model was built from old drawings available on
DFS (\\cern.ch\dfs\Workspaces\o\Old Drawings\Complexe_PS\PS\MPS-2A68).

The present CST model is stored on
[EDMS](https://edms.cern.ch/ui/#!master/navigator/document?P:100327819:100451893:subDocs).

<details>
  <summary>KFA04 in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/big/ss04i-a.jpg">
</details>

## Subfolders and usage

- Tsutsui

The model was computed using ImpedanceWake2D with the following parameters

    - a=0.056 half width [m]
    - b=0.037 half heigth (center-ferrite border) [m]
    - d=0.040 half heigth (center-PEC border including ferrite) [m]
    - L=0.615 Kicker's length [m]

- CST Raw Data

The wakefield is obtained with CST simulations. The model was
built based on 2D drawings.

The output result is the impedance in frequency domain and
can directly be used in particle simulations. The impedance
is noisy, to obtain a smoother response the fit with resonators
can be preferred.

<details>
  <summary>CST model view</summary>
  <img src="cst_raw_data/cst_model_images/KFA04_CST_Model_View1.png">
  <img src="cst_raw_data/cst_model_images/KFA04_CST_Model_View2.png">
  <img src="cst_raw_data/cst_model_images/KFA04_CST_Model_View3.png">
  <img src="cst_raw_data/cst_model_images/KFA04_CST_Model_View4.png">
</details>

- Resonators

The result of the CST simulation is fitted with resonators to get
the more global features of the impedance.

The broadband component of the impedance is fitted with 2 resonators.
Low and intermediate frequency resonances are also fitted with
2 extra resonators.

## Plots

<details>
  <summary>Comparing Tsutsui model with CST</summary>
  <img src="Tsutsui/comparison_tsu_cst.png">
</details>

<details>
  <summary>Resonator fit</summary>
  <img src="Resonators/fitted_broadband.png">
  <img src="Resonators/fitted_broadband_realimag.png">
  <img src="Resonators/fitted_broadband_extra_resonances.png">
</details>

## History

To be completed

## Further work

To be completed

- Confirm 8C11 material data

## References

To be completed

### Simulations

- [S. Persichelli, The beam coupling impedance model of CERN Proton Synchrotron](https://cds.cern.ch/record/2027523)
- [B. Popovic, LIU-PS Beam Dynamics WG meeting \#27](https://indico.cern.ch/event/801198/contributions/3329723/attachments/1803818/2943005/PS_Imped_Model_Update_280219_FINAL.pdf)

## Credits

To be completed

- Tsutsui: BE-ABP-HSC Serena Persichelli, Elias Metral, Nicolo Biancacci or Benoit Salvant
- Model: B. Popovic
- 8C11 Material: S. Persichelli
