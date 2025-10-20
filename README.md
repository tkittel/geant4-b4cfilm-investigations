# Simple B4C film simulations with Geant4

Very simple modelling of B4C thin film detectors, in order to gauge and parameterise the probability of a neutron absorption in the film leading to a detection event in the counting gas.

The simulations require a conda environment (see [conda.yml](conda.yml)), and can be launched with the script:

```sh
./bin/rung4sim --vis # to visualise
./bin/rung4sim -n 1e7 --out results.json #produce results
./bin/rung4sim --help # see more options
./bin/analyse results.json #analyse results
```

The [data](data/) directory contains some already generated output files (note that the json data includes details of the configuration with which it was produced).

The [recipes](recipes/) directory contains a function implementing the resulting parameterisation, with files available for both C and Python.

For more details, refer to the code in [src/film.cc](src/film.cc) and [bin/analyse](bin/analyse)

*Project initiated by T. Kittelmann, October 2025.*
