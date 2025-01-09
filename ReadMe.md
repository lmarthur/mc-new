## Outline
The intended use of this repository is to calculate the aerodynamic parameters of a reentry vehicle with a given geometry, assuming Newtonian flow. To do this, the geometry must first be specified in the input file. 

## TODO: 
- Write bash scripts for each geometry case
- Write documentation for github repo
- Write f2py wrapper for MC-NEW
- Implement test suite via f2py wrapper
- Install code on Engaging cluster
- Generate SWERVE plots/write matplotlib scripts
- Generate conic, biconic, and SWERVE geometries
- Connect MC-NEW output to pytrajlib input
- Test hypersonic flow convergence w/ empirical data/Newtonian approx
- Generate conic moments w/ perturbed CoM
- Look for analytic/empirical CAV moments
- Look into whether MC-NEW can handle reaction mass control
- Attempt to structure non-axisymmetric geometries
- Compare SWERVE moments from zero roll to 45 deg roll
- Generate stability plots that show how time constant changes with different geometry parameters (flap thickness, length, radii, etc.)
- Generate response plots from trajectory variations
    - Instantaneous lateral position anomaly
    - Step function acceleration anomaly
    - Density-coupled acceleration anomaly
    - Simulated ablation error?
- Generate accuracy numbers w/ Earth-GRAM wind/density errors
- Configure wedge geometry?
- Generate roll time-constants for wedge
- Formulate 6-DOF guidance law?
- Try learned guidance law (PySR? NN?)


