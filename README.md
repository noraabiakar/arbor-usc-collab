# Granule and Basket Cells
This repository contains the Neuron and Arbor implementations and simulations of
1) A single basket cell
2) A single granule cell

## Morphologies
1) The basket cell is a single soma model
2) The granule cell reads the morphology from **output0_updated.swc** in the morphologies directory

## Density Mechanisms
The basket and granule cells both have the following mechanisms inserted on the soma and dendrites:
1) borgka.mod
2) cagk.mod
3) cat.mod
4) ccanl.mod
5) gskch.mod
6) ichan2.mod
7) lca.mod
8) nca.mod

## Parameters
Various parameters can be set and will be read by Arbor and Neuron in param_basket.json and param_granule.json

## Stimulus
A series of spikes is delivered to a synapse on the basket/granule cell.
The spike times are read from param_basket.json and param_granule.json.

## Running the examples

### Neuron:
```
$ cd neuron
$ nrnivmodl modfiles
$ python2 test_granulecell.py ../param_granule.json
OR
$ python2 test_basketcell.py ../param_basket.json
```

Running the examples will generate a plot

### Arbor:
**Installing Arbor**:
* Install Arbor from the following branch:
https://github.com/noraabiakar/??
* For help installing Arbor refer to:
https://arbor.readthedocs.io/en/latest/install.html

**Compiling the example**:
```
$ cd arbor/granule
OR
$ cd arbor/basket

$ mkdir build && cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/arbor/installation -DCMAKE_BUILD_TYPE=release
$ make
```

**Running the example**:
```
$ ./granule ../../../param_granule.json
OR
$ ./basket ../../../param_basket.json
```

**Plotting the results**:

Running the examples will generate **voltages.json**. To plot the results, you can use the provided **tsplot** script:
```
$ python2 tsplot.py arbor/granule/build/voltages.json
OR
$ python2 tsplot.py arbor/basket/build/voltages.json
```
