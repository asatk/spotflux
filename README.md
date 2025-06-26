# Finite Differences Flux Transport Model (C, python)

Requires gcc, make, and a (python) conda environment with numpy and matplotlib

### Change Simulation Behavior (TEMPORARY)
edit file src/constants.c
Note: this behavior will be changed. constants.c will hold default parameter
values and a configuration file in the top level of the directory will be used
to tweak simulation behavior.

### Compile Simulation
make clean
make

### Run Simulation
./main

### Generate Plots
python plot.py
