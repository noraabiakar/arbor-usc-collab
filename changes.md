# Changes
This document lists the changes made to the original granule and basket cell models from USC


####The mechanism files 
The mechanisms used in the models are implemented in the following 8 `.mod` files
1) borgka.mod
2) cagk.mod
3) cat.mod
4) ccanl.mod
5) gskch.mod
6) ichan2.mod
7) lca.mod
8) nca.mod

Arbor has it's own `modcc` compiler that doesn't support all the same features `nrnivmodl` supports. 
One of these features is the `TABLE` command which generates a look-up table to speed-up a simulation 
in Neuron. We removed this command from all the mechanism files to guarantee accurate comparison 
of Arbor and Neuron's voltage measurements. This has the positive effect of generating more 
accurate results in Neuron that match Arbor's almost perfectly, but it slows down the Neuron 
simulation around 1.5x.

The original mechanism files can be found in `neuron/modfiles` and the edited modfiles without `TABLE` 
are in `neuron/modfiles/arb`. We suggest using `neuron/modfiles/arb` for arbor-neuron verification. 

####Fast forwarding to steady state
The original `negative_init.hoc` tries to fast forward the simulation by starting the it from a 
negative time point `t = -1e10` and taking a single big time step `dt=1e9`. This is to ensure that
the voltage has reached a steady state before starting the real simulation. This was removed
for a more accurate comparison between Arbor and Neuron and to study how they progress towards 
the steady state. The effect of fast-forwarding the simulation can be mimicked in Arbor if necessary. 

####NMDA/GABA synapses
In the provided setup of the granule and baskets cells, `AMPA`, `NMDA` and `GABA` synapse groups are 
defined. Using the provided setup and mechanisms, the 3 groups are identical, and we can only use 
1 of them at a time. Therefore, we simplified the code by having only one set of synapse groups `AMPA`.
This doesn't affect the results of the simulations, and can be easily reversed in both Arbor and Neuron
when `NMDA` and `GABA` become relevant. 

####Parameterizable variables
Many of the input variables to the simulation were hard coded in the python files. We moved many of the 
variable assignments into `.json` files, namely: `param_granule.json` and `basket_json.json`. This makes 
it easier to perform some parameter sweeps in the simulations as well as makes it easier to guarantee that
Arbor and Neuron have the same values as they both read the same file to initialize the simulation. 

To see a full list of the parameterizable variables check `param_granule.json` or `param_basket.json`. 

####Input Spikes
In the provided setup, random spikes are generated as stimuli to the simulation. Because we can't guarantee
that the same random numbers will be generated in Arbor and Neuron, we hard code the spike times in the
`param_granule.json` and `param_basket.json` files and feed them to both the Arbor and Neuron simulations. 


####Parameters setting on segment points 
In the original code, some of the mechanism parameters are set at specific points along a segment:
```
for sec in cell.granuleCellLayer:
    if len(cell.granuleCellLayer[sec]) > 0:
        for norm_dist in cell.granuleCellLayer[sec]:
            sec(norm_dist).gnatbar_ichan2 = 0.018*cell.a1
	    ...
```

Here `gnatbar_ichan2` is set at the point `norm_dist` on the section `sec`. In the default version 
of the original code, `sec` will always have `nseg=1`, so `gnatbar_ichan2` will be set on the entire 
section. But, if `sec` has `nseg>1`, this will result in `sec` having a different value of 
`gnatbar_ichan2` on different segments along `sec`. 
In Arbor, mechanism parameters can only be set at the granularity of a segment, therefore we change 
the previous code to the following: 

```
for sec in cell.granuleCellLayer:
    if len(cell.granuleCellLayer[sec]) > 0:
        sec.gnatbar_ichan2 = 0.018*cell.a1
	    ...
```
