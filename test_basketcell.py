# Test multiple synaptic activation of cell models

from neuron import h
import cell
import numpy as np
import time as cookie
import pylab as plt
import pickle
np.random.seed(149)

h.load_file("stdrun.hoc")

##################
# Creating cells #
##################
ID = 0
location = (0,0)
synvars = {}
synvars['type'] = "E2"
fname_morph = 'morphologies/output0_updated.swc'
modeltype = 'Single'
celltype = 'basketcell'

neuron_DG = cell.Cell(ID,location,synvars,celltype,fname_morph,modeltype)

################################
# Create spike times for input #
################################
tstop = 10000 # unts: ms
num_input = 1
frequency = 5 # units: Hz

vecstims = [h.VecStim() for ii in range(num_input)]
evecs = [h.Vector() for ii in range(num_input)]
vec = []

for ii in range(num_input):
    intervals = []
    mu = 1000./frequency # Convert to ms
    elapsed_time = 0
    flag = 1
    while flag:
        roll = np.random.uniform(0,1)
        interval = -mu*np.log(roll)
        
        intervals.append(interval)
        elapsed_time += interval
        
        if elapsed_time > tstop:
            flag = 0
            spikes = np.cumsum(intervals)[:-1]
            for spike in spikes:
                evecs[ii].append(spike)
                vec.append(spike)
            
            vecstims[ii].play(evecs[ii])

for v in vec:
    print(v)

#####################
# Connecting inputs #
#####################
w_MEA_av = 1.438995e-04
net_cons = []
for ii in range(num_input):
    nc = h.NetCon(vecstims[ii], neuron_DG.synGroups['AMPA']['soma'][0])
    neuron_DG.synGroups['AMPA']['soma'][0].tau1 = 0.1
    neuron_DG.synGroups['AMPA']['soma'][0].tau2 = 16.423438
    neuron_DG.synGroups['AMPA']['soma'][0].e = 0
    nc.weight[0] = w_MEA_av
    nc.delay = 0
    net_cons.append(nc)

################################
# Setting up vectors to record #
################################
v = h.Vector()
v.record(neuron_DG.soma(0.5)._ref_v)

t = h.Vector()
t.record(h._ref_t)

#########################
# Setting up simulation #
#########################
h.v_init = -70
h.t = 0
h.dt = 0.025
h.celsius = 35.0
h("tstep = 0")
h("period = 2")
h.tstop = tstop
h("steps_per_ms = 10")
h.load_file('hoc_files/negative_init.hoc')

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
