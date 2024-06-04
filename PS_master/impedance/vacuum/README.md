# Impedance of Vacuum equipment

The impedance linked to vacuum equipment (aperture changes, cavity gaps,
ground loops...)

The aperture changes and cavity gaps are mostly responsible for high frequency
resonances, leading to constant inductive impedance ImZ/n at low frequency
and can drive microwave instabilities.

The ground loops with rf bypasses are responsible for low frequency resonances
mostly contributing to energy loss at the single bunch level.

<details>
  <summary>Individual Impedance</summary>
  <img src="overview/overview_indiv_vacuum.png">
  <img src="overview/overview_indiv_vacuum_zoom.png">
</details>

<details>
  <summary>Summed resistive impedance</summary>
  <img src="overview/overview_real_summed.png">
  <img src="overview/overview_real_summed_zoom.png">
</details>

<details>
  <summary>Summed reactive impedance, high freq. inductive contribution</summary>
  <img src="overview/overview_imag_summed_highfreqonly_zoom.png">
</details>

<details>
  <summary>Summed reactive impedance, all contributions</summary>
  <img src="overview/overview_imag_summed.png">
  <img src="overview/overview_imag_summed_zoom_x.png">
  <img src="overview/overview_imag_summed_zoom_xy.png">
</details>

## Equipment around main Magnet Units (MU)

There are 100 MU sections in the PS.
Each section has an upstream assembly and a downstream assembly.

The upstream assembly are bellow sections and are listed as follows:

- Upstream Assemblies (99):
    - Standard (67) > **Modeled, without the flange**
    - Enlarged (25)
    - 'Exotic' Variations (7)
        - Wide, special, etc

Internal & External are referring to which side the armature sits (internal=>inner side)

- Pumping Manifold Types (100):
    - Internal (29) & External (37) > **Modeled**
    - Enlarged Internal (10) & External (8)
    - 'Exotic' Variations (16)
        - Wide, Narrow etc.

Almost each downstream section (99) has pumping manifold slot:

- Pickup Assembly (62)
    - Standard (40) > **Modeled**
	- Large (12)
- Group Pumping Coupling (8)
- Actuator of SEM Grid Detector (3)
- Cover/Empty (36) > **Modeled**

Not all variations are modeled. The present approximation is to take the most represented element
as a sufficient representation of the elements that are not yet modeled.

The MUs are connected to the straight sections with isolated vacuum flanges. The isolation of the
flanges is responsible of ground loops leading to low frequency resonance. To reduce the resonance,
RF bypasses are installed on the flange.

New bumper magnets and vacuum chambers will be installed in straight sections 40, 41, 43, 44
as described in [ECR](https://edms.cern.ch/ui/#!master/navigator/document?D:100485334:100485334:subDocs).
The steps in the vacuum chambers is responsible of impedance

## Sector valves

There are 10 identical sector valves located after straight sections 00,10,20,30,40,50,60,70,80,90.

## Flanges and step transitions

There are 200+ flanges (**survey to be refreshed**). The flanges lead to small
cavity gaps between vacuum chambers. Their contribution is partially included
in the modeling of the equipment around the MUs (**survey to be refreshed**).

Around each main magnet units, the flanges are isolated to reduce the beam induced current in ground loops. The insulation
lead to resonances at low frequency, which are mitigated using RF-bypasses. The impedance
contribution at low frequency is summed independently from the impedance of
the cavity gap at high frequency.

The transitions from one vacuum chamber geometry to another is responsible
of impedance which is presently approximated as a constant ImZ/n contribution.

## Standalone elements

The flanges next to the equipment in the straight sections are partially included. Some of these
are directly modeled with the equipment, while some are not.

The Y-chambers close to the septa are not included.

A better list of these elements should be established.

## Vacuum chamber

The impedance of the vacuum chamber itself is included in the resistive wall impedance model. 

## References

To be completed

- [B. Popovic, 29th Impedance Working Group](https://indico.cern.ch/event/797962/contributions/3316173/attachments/1794685/2925021/IWG_AOB_Inj_Bumpers_120219_IWG.pdf)
- [B. Popovic, Longitudinal limitations with LIU-PS RF upgrades and mitigation strategy](https://indico.cern.ch/event/750790/contributions/3108016/attachments/1719965/2776247/RG_PS_Impedance_Model_Meeting.pdf)
- [B. Popovic, SPS Injection Losses Review](https://indico.cern.ch/event/672967/contributions/2753563/attachments/1567554/2471183/PS_Impedance_Status_Injection_Losses_Meeting_301117.pdf)
- [B. Popovic, LIU-PS Beam Dynamics WG meeting #1](https://indico.cern.ch/event/662292/contributions/2704441/attachments/1516648/2367108/BP_PS_Impedance_Model_Update_310817_Final.pdf)
- [D. Boussard, Observation of microwave longitudinal instabilities in the CPS](https://cds.cern.ch/record/872559/files/cer-002556319.pdf)

## Credits

To be completed
