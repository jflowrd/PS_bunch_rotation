# Impedance of BFA21

The BFA21 is part of the Continuous Transfer (CT) extraction
(five turn extration system, old scheme).

The BFA21 is composed of two elements, the staircase (S) and the pedestal (P).

The CT extraction scheme will not be used after the LS2 (2019-2020), and the
equipment will be removed. Both BFA21P and BFA21S are removed during LS2.
More details in the following [ECR](https://edms.cern.ch/ui/file/1981131/0.1/PS-LJ-EC-0007-00-10.pdf).

The first estimation was done using a Tsutsui model available on the following [PS impedance webpage](http://impedance.web.cern.ch/impedance/)
and on the [IRIS PS repository](https://gitlab.cern.ch/IRIS/PS_IW_model/-/tree/master/Impedances/Longitudinal).

Due to the removal, no new CST model was designed. *Need to determine whether
BFA09 and BFA21 have the same design*.

<details>
  <summary>BFA21 in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss21.jpg">
</details>

## Subfolders and usage

- Tsutsui

The model was computed using ImpedanceWake2D with the following parameters

BFA21S	

    - a=0.079 half width [m]
    - b=0.02625 half heigth (center-ferrite border) [m]
    - d=0.035 half heigth (center-PEC border including ferrite) [m]
    - L=0.43 Kicker's length [m]
    
BFA21P

    - a=0.079 half width [m]
    - b=0.02625 half heigth (center-ferrite border) [m]
    - d=0.020 half heigth (center-PEC border including ferrite) [m]
    - L=0.54 Kicker's length [m]

- Resonators

The Tsutsui model is fitted with resonators for more flexible calculations.

## Plots

<details>
  <summary>Tsutsui model</summary>
  <img src="Tsutsui/tsutsui_model.png">
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

- Check whether BFA09 and BFA21 share the same design.

## References

To be completed

### Simulations

- [S. Persichelli, The beam coupling impedance model of CERN Proton Synchrotron](https://cds.cern.ch/record/2027523)

## Credits

To be completed

- Tsutsui: BE-ABP-HSC Serena Persichelli, Elias Metral, Nicolo Biancacci or Benoit Salvant
- 8C11 Material: S. Persichelli
