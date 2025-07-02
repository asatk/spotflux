# Finite Differences Flux Transport Model (C, python)

Requires gcc, make, and a (python) conda environment with numpy and matplotlib

### Configure Simulation
Copy the default parameters at `src/constants.c` to a config file, e.g.,
`config.c`:

```cp src/constants.c config.c```

Edit this configutation file with the parameters you want changed. Documentation
of these parameters is located at `src/constants.h`.

### Compile Simulation
```make clean```

```make```

If not using the default simulation parameters, specify your configuration file
using instead:
```make CONSTS=config.c```

### Run Simulation
```./main```

A 2D map of the surface magnetic field is saved every `freq` timesteps to file
`bfld.dat`. This behavior can be changed in `config.c`.

### Generate Plots
```python plot.py```

Python will read in the raw data file `bfld.dat` produced by `main`. It will
also read in the parameter set given in a config file which is by default
`config.c`. The config file identified at the beginning of `plot.py` as `cfname`
should correspond to the data file you are analyzing.

Note: the last line of the script will save a .gif file of all the
frames. You can change the name or comment this line out as it takes some time
to render the .gif image especially if there are a lot of frames.

### Analysis
```script.py```

A basic Python script is provided for users to quickly grasp how the data are
formatted and read in to a Python application. This is a barebones script from
which more meaningful analyses can build off.
