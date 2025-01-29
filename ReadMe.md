## Outline
The intended use of this repository is to calculate the aerodynamic parameters of a reentry vehicle with a given geometry, assuming Newtonian flow. To do this, the geometry must first be specified in the input file. 

## Manual
To run a simulation batch, run some version of the following command from the Example subdirectory: 

```bash run.sh -i input.flap -d runs -p pitch -s -5.0 -e 5.0 -n 21```

## TODO: 
- [X] Write bash scripts for each geometry case 
- [ ] Write documentation for github repo 
- [ ] Generate SWERVE plots/write matplotlib scripts
- [ ] Connect MC-NEW output to pytrajlib input
- [ ] Test hypersonic flow convergence w/ empirical data/Newtonian approx
- [ ] Look for analytic/empirical CAV moments
- [ ] Look into whether MC-NEW can handle reaction mass control
- [ ] Compare SWERVE moments from zero roll to 45 deg roll
- [ ] Generate stability plots that show how time constant changes with different geometry parameters (flap thickness, length, radii, etc.)
- [ ] Generate response plots from trajectory variations
- [ ] Generate accuracy numbers w/ Earth-GRAM wind/density errors
- [ ] Attempt to structure non-axisymmetric geometries
- [ ] Configure wedge geometry?
- [ ] Generate roll time-constants for wedge
- [ ] Formulate 6-DOF guidance law?

## BUGS:
- [ ] Running a scan over Mach number leaves the Mach number in the input files unchanged (problem w/ modify_input.py?)


