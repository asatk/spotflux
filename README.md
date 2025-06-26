# Finite Differences Flux Transport Model (C, python)

Requires gcc, make, and a (python) conda environment with numpy and matplotlib

### Configure Simulation
Copy the default parameters at `src/constants.c`:

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

### Generate Plots
```python plot.py```
