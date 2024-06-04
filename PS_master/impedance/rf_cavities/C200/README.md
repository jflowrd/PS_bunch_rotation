# Impedance of 200 MHz cavities

Those cavities are used for longitudinal emittance blow-up and bunching before 
SPS extraction (pill-box with ceramic gap).

There are six cavities in total in the same location in the ring. The cavities
are tuned to h=420 for proton beam at 14 GeV (for fixed target beam in the SPS).
In operation the Q of the cavity is reduced by 3 PIN lines.

The design values of the cavities (see reference hand written note,
with a central frequency of fr=199.6 MHz and R/Q=28.5 Ohm):

| Case 							| Q		| Rs [kOhm]	|
| ----------------------------- |:-----:|:---------:|
| Naked cavity  				| 1900	| 54.15 	|
| 2 magnetic loop (operational)	| 1000	| 28.5 		|
| 1 line 50 Ohm 				| 328 	| 9.3 		|
| 2 lines 50 Ohm 				| 200	| 5.7 		|
| 3 lines 50 Ohm (damped)		| 130 	| 3.7 		| 

The CST simulations including the input coupler and the 3 PIN lines gave the 
following results:

| Parameter \ Frequency [MHz] 			| 199.9		|
| -------------------------------------	|:---------:|
| Shunt resistance R [kOhm]				| 1.092		|
| Quality factor Q 						| 40 		|
| R/Q [Ohm] 							| 27.3		|

According to measurements, the loaded quality factor is Q=134 when damped by
the PIN diodes and Q=971 in operational mode. These are the values presently
kept in the impedance model.

<details>
  <summary>C200 in the ring</summary>
  <img src="http://cern.ch/psring/psring/pictures/fullsize/ss06.jpg">
</details>

## Usage

- Resonator model based on CST simulations (main harmonic)

Using the values as obtained by S. Persichelli (see reference), tuned
to the operationnal frequency. The Q value was selected to the be the measured
ones.

## Further work

Tasks to be determined (beam measurements of the impedance ?)

## References

### Simulations and measurements
- [S. Persichelli, The PS transverse impedance model](https://indico.cern.ch/event/326091/#2-the-ps-transverse-impedance)
- [S. Persichelli, The beam coupling impedance model of CERN Proton Synchrotron](http://cds.cern.ch/record/2027523)

### Design
- [D. Boussard, The PS 200 MHz RF system : present situation and future prospects](https://cds.cern.ch/record/118981/)

### Other
- Hand written note from Jamsek to Cappi, 08/01/1985
- [Fifty years of the CERN Proton Synchrotron : Volume 1](https://cds.cern.ch/record/1359959)
