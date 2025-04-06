# capsid_percolation_modeling
 Collection of codes with Monte Carlo simulations for modeling virus shell (capsid) fragmentation

![](https://github.com/tanwarmayank524/capsid_percolation_modeling/blob/main/mygif1.gif)

 ## Background
In this work, we use graph theory to model the stability of capsids. The graph geometry gives the shape of the initial capsid, the nodes represent the individual subunits of the capsid, and the edges between the nodes quantify the interactions between the subunits.
The interactions between the subunits are obtained with molecular dynamics simulations, and the corresponding energetics are fed into the Monte Carlo code to simulate the fragmentation.

The Monte Carlo code looks at different aspects of the fragmentation process, including the effect of:
1. The energy of interaction between subunits
2. The size and shape of the initial capsid
3. The different modeling approaches to include the surroundings
4. The effect of surroundings interactions


## code 1
code1.py gives the evolution of the fragmentation process in the graph for a single run. The script can be run with:

```python code1.py grid_x grid_y number_of_runs fragmentation_threshold_fraction```

Here, grid_x(e.g., 6) and grid_y(e.g., 10) are the number of grid points in the x and y directions in the periodic graph. 
number_of_runs(e.g., 200) x the number of nodes gives the number of times you want to run the single simulation. 
fragmentation_threshold_fraction(e.g., 3) gives the termination condition when empty nodes are more than a certain fraction of total nodes.      


## code 2
code2.py gives the average number of remaining graph edges after running a certain number of independent Monte Carlo simulations. The script can be run with:

```python code2.py grid_x grid_y number_of_runs fragmentation_threshold_fraction number_of_simulations```

number_of_simulations(e.g., 100) gives the number of independent Monte Carlo simulations for obtaining the average behavior of the fragmentation process.


This work is in progress, and more details will be added soon!

# Authors
Mayank Tanwar
email: tanwa008@umn.edu
GitHub: tanwarmayank524

Collaborators: Jonathan Zajac (zajac028@umn.edu) and Prof. Sapna Sarupria (sarupria@umn.edu)
