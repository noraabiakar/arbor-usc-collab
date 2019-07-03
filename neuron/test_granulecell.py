# Test multiple synaptic activation of cell models

from neuron import h
import cell
import numpy as np
import time as cookie
import pylab as plt
import pickle
import json
import sys

with open(sys.argv[1]) as json_file:
    in_param = json.load(json_file)

np.random.seed(149)

h.load_file("stdrun.hoc")

##################
# Creating cells #
##################
ID = 0
location = (0,0)
synvars = {}
synvars['type'] = "E2"
fname_morph = in_param["morph_file"]
modeltype = 'Multi'
celltype = 'granulecell'

neuron_DG = cell.Cell(ID,location,synvars,celltype,fname_morph,in_param, modeltype)

################################
# Create spike times for input #
################################
num_input = 1
tstop = in_param["run_time"] # unts: ms

vecstims = [h.VecStim() for ii in range(num_input)]
evecs = [h.Vector() for ii in range(num_input)]

for ii in range(num_input):
    for spike in in_param["spikes"]:
        if spike > tstop:
            break
        evecs[ii].append(spike)
    vecstims[ii].play(evecs[ii])

#####################
# Connecting inputs #
#####################
w_MEA_av = in_param["weight"]
net_cons = []
for ii in range(num_input):
    choice = in_param["syn_id"]
    nc = h.NetCon(vecstims[ii], neuron_DG.synGroups['AMPA'][in_param["syn_layer"]][choice])
    neuron_DG.synGroups['AMPA'][in_param["syn_layer"]][choice].tau1 = in_param["tau1_syn"]
    neuron_DG.synGroups['AMPA'][in_param["syn_layer"]][choice].tau2 = in_param["tau2_syn"]
    neuron_DG.synGroups['AMPA'][in_param["syn_layer"]][choice].e = in_param["e_syn"]
    nc.weight[0] = w_MEA_av
    nc.delay = 0
    net_cons.append(nc)

################################
# Setting up vectors to record #
################################
v = h.Vector()
v.record(neuron_DG.c.soma[0](0.5)._ref_v)

t = h.Vector()
t.record(h._ref_t)

#########################
# Setting up simulation #
#########################
h.v_init = in_param["vinit"]
h.t = 0
h.dt = in_param["dt_neuron"]
h.celsius = in_param["temp"]
h("tstep = 0")
h("period = 2")
h.tstop = tstop
h.steps_per_ms = 1/in_param["dt_neuron"]
h.load_file('hoc_files/negative_init.hoc')

print h.secondorder


##################
# Run simulation #
##################
print("Starting...!")
ST = cookie.time()
h.run()
ET = cookie.time()-ST
print("Finished in %f seconds" % ET)

########
# Plot #
########
_=plt.plot(t,v)
_=plt.xlabel('Time (ms)')
_=plt.ylabel('Somatic Voltage (mV)')
plt.show()

