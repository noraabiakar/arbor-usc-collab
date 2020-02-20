## Running the network

**Installing Arbor**:
* Install Arbor from the following branch:
https://github.com/noraabiakar/usc_collab/tree/testing/new-revpot
* For help installing Arbor refer to:
https://arbor.readthedocs.io/en/latest/install.html

**Compiling the example**:
```
$ cd network
$ mkdir build && cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/arbor/installation/lib/cmake/install -DCMAKE_BUILD_TYPE=release
$ make
```

**Running the example**:
```
$ ./arbor-network ../param_sim.json
```

**Plotting the results**:

Running the examples will generate **voltages.json**. To plot the results, you can use the provided **tsplot** script:
```
$ python2 tsplot.py arbor/network/build/voltages.json
```
