# PS longitudinal impedance model

The cutoff frequency of the main elliptic vacuum chamber (146mm x 70mm) for the fundamental longitudinal mode is 2.54 GHz.

## Impedance sources included

- RF cavities
    - 10 MHz (with fast feedback, with/without one-turn delay feedback)
    - 20 MHz (with fast feedback, with/without multi-harmonic feedback)
    - 40 MHz (with fast feedback, with/without multi-harmonic feedback)
    - 80 MHz (with fast feedback, with/without multi-harmonic feedback)
    - 200 MHz
    - Finemet cavity
- Kickers
    - KFA
		- KFA04
        - KFA13
        - KFA21
        - KFA45
        - KFA71
        - KFA79
    - BFA
		- BFA09
		- BFA21
- Septa
    - SMH42
- Vacuum equipments
    - Sector valves
    - MU sections, downstream assembly (pumping manifolds, empty and with CODD inserted) 
    - MU sections, upstream assembly (bellows)
- Beam Instrumentation
    - Vertical BGI (Simulations Done, Measurements Forthcoming)
    - Stripline BPM (SD72)
    - Wall Current Monitor (WCM) (SD03)
    - Wire Scanners (Model Availabilty?)
- Miscellaneous
    - Internal beam dumps (before and after LS2)
    - Tranverse Feedback Kicker (TFB Kicker) (SD97)
- Space charge
- Resistive wall

## Missing impedance sources to be included

- Septa
    - TPS15 (Added during LS1, was modeled in CST)
    - SMH16 (To be Approximated)
    - SEH23 (To be Approximated)
    - SMH26 (To be Approximated)
    - SMH57 (To be Approximated)
    - SMH61 (To be Approximated)
    - SEH31 (To be removed during LS2, was modeled in CST)
- Vacuum equipments
    - Variations of upstream/downstream assembly
    - Flanges/Step transitions (need better survey)
    - Possible missing element to be identified (e.g. y-chambers)
- Beam instrumentation
    - Horizontal BGI (Model Availabilty?)
    - SEM grids (Model Availabilty?)

## History

### During LS2

Added:

- Vertical BGI

Upgraded:

- New amplifiers on C10
- Multi-harmonic feedback on C20/C40/C80

Replaced:

- New beam dump in place of the old ones

Removed:

- CT equipement to be removed during LS2
    - Kickers: BFA21S, BFA21P and BFA09S
    - Septa: SES31, SEH31
    - Orbit bumpers: magnet type D205 in SS27 and D210 in SS35

## Useful links

### General info

- [Proton Synchrotron Straight Section Pictures](http://psring.web.cern.ch/psring/pscomplex.shtml)
- [Jeroen Belleman's web pages](https://jeroen.web.cern.ch/jeroen/)

### Layout

- [EDMS CAD Navigator](https://edms.cern.ch/nav/SMTI137:450976)
- [PS Ring in Layout DB](https://layout.cern.ch/elements?id=2180074&parentId=0&circuitId=0&version=LS2&navigator=MACHINE&tab=BASIC_ELEMENT)
- [MAD Sequence from Layout DB](https://layout.cern.ch/reports/mad?machineId=2180074&versionId=34464591)
- [PS Ring on GIS Portal](https://gis.cern.ch/gisportal/Machine.htm?pp=[42723712])
- [PS Ring on Google maps](https://www.google.com/maps/@46.2326176,6.0478399,2a,75y,241h,83.59t/data=!3m6!1e1!3m4!1s01xLBETOwHzt8Frx6D995A!2e0!7i13312!8i6656)
- [Straight sections drawings](https://edms5.cern.ch/cdd/plsql/c4w_folders.display_contents?cookie=2404844&p_folder_id=97676)
- [Magnet units drawings](https://edms5.cern.ch/cdd/plsql/c4w_folders.display_contents?cookie=2404846&p_folder_id=97574)

### Impedance

- [CST models and measurements on EDMS](https://edms.cern.ch/project/CERN-0000197012)
- [CERN impedance webpage - Transverse](https://impedance.web.cern.ch/impedance/)
- [PS optics repository - Gitlab](https://gitlab.cern.ch/acc-models/acc-models-ps)
- [PS transverse impedance model - Gitlab](https://gitlab.cern.ch/IRIS/PS_IW_model)
- [PS transverse impedance model 2 - Gitlab](https://gitlab.cern.ch/IRIS-2/PS-IW-model)

### Optics

- [PS optics repository - Web](http://acc-models.web.cern.ch/ps/)



