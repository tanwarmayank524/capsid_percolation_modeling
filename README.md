# capsid_percolation_modeling
 Collection of codes with Monte Carlo simulations for modeling capsid percolation

 ## volume fraction distribution
volume_fraction_distribution gives the amount of volume and the fraction of volume around a reference atom occupied by the rest of the atoms for a periodic ab initio molecular dynamics simulation. This distribution quantifies the local environment by treating atoms as van der Waals spheres and calculating how much volume is occupied by surrounding atoms for a reference atom as a function of distance. The code gives the volume occupied and the fraction relative to the total volume. The script can be run with:

```python volume_fraction_distribution.py *XDATCAR.xlsx atom1 grid OUTPUT.png```

Here, atom1 is the reference atom for creating the volume fraction distribution, grid is the gridpoint size for generating the distributions, and the OUTPUT file gives volume and volume fraction distributions. *XDATCAR.xlsx is the input file.


# Authors
Mayank Tanwar
email: tanwa008@umn.edu
GitHub: tanwarmayank524

Collaborators: Jonathan Zajac (zajac028@umn.edu) and Prof. Sapna Sarupria (sarupria@umn.edu)
