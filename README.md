# Coupled cluster doubles on the simple pairing model
Finding the ground state energy of the pairing model with coupled cluster doubles, see [project report](https://github.com/stiandb/QuantumComputing/blob/master/PDF/QCProject.pdf).

### Run CCD code
Run the CCD_pairing.py script as follows
```
python CCD_pairing.py {h} {p} {xi} {g} {alpha}
```
where
* h - number of hole states. That is, the number of particles in your system.
* p - number of particle states. That is, the number of virtual orbitals in your basis.
* xi - The level spacing.
* g - Interaction strength.
* alpha - Iterative mixing constant a=(0,1]. Optional, it is set to 0.5 when avoided.


### Reproduce results from project
To reproduce the results in the project report. Run the plotter.py script as follows
```
python plotter.py {h} {p}
```
where
* h - Number of hole states.
* p - Number of particle states.

