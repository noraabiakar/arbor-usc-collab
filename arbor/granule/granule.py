import sys    
import arbor
from arbor import mechanism as mech
from arbor import location as loc
import matplotlib.pyplot as plt

with open(sys.argv[1]) as json_file:
    params = json.load(json_file)

vinit = params["vinit"]
temp  = params["temp"]
tsim  = params["run_time"]
dt    = params["dt_arbor"]
spikes= params["spikes"]
plot_to_file=False
seed = 42

# load morphology from swc file
tree = arbor.load_swc(params["morph_file"])
morph = arbor.morphology(tree, spherical_root=True)

# define named regions and locations
labels = { 'all': '(all)',
           'soma':       '(tag 1)',
           'dend':       '(tag 3)',
           'granuleCellLayer': '(intersect (region "dend") (z_dist_from_root_gt(0    )) (z_dist_from_root_le(6.68 )))',
           'innerThird'      : '(intersect (region "dend") (z_dist_from_root_gt(6.68 )) (z_dist_from_root_le(20.04)))',
           'middleThird'     : '(intersect (region "dend") (z_dist_from_root_gt(20.04)) (z_dist_from_root_le(40.08)))',
           'outerThird'      : '(intersect (region "dend") (z_dist_from_root_gt(40.08)) (z_dist_from_root_le(66.9 )))',
           'soma_syn'        : '(location 0 0.5)',
           'granule_syn'     : '(uniform (region "granuleCellLayer") 0 0 ' + str(seed) + ')',
           'inner_syn'       : '(uniform (region "innerThird")       0 0 ' + str(seed) + ')',
           'middle_syn'      : '(uniform (region "middleThird")      0 0 ' + str(seed) + ')',
           'outer_syn'       : '(uniform (region "outerThird")       0 0 ' + str(seed) + ')'}

ldict = arbor.label_dict(labels)

# construct cell
cell = arbor.cable_cell(morph, ldict)

# set cell properties
cell.set_properties(Vm=vinit, tempK=temp+273.15)
cell.set_ion(  'k', rev_pot=params["ek"])
cell.set_ion( 'ca', ext_con=params["cao"])

cell.paint('soma', mech('ichan2', {"gnatbar":0.120*params["gnatbar_ichan2"],
                                   "gkfbar":0.016*params["gkfbar_ichan2"],
                                   "gksbar":0.006*params["gksbar_ichan2"],
                                   "gl":0.00004*params["gl_ichan2"],
                                   "el":params["el_ichan2"]}))
cell.paint('soma', mech('borgka', {"gkabar" : 0.001   *params["gkabar_borgka"]}))
cell.paint('soma', mech('nca'   , {"gncabar": 0.001   *params["gncabar_nca "]}))
cell.paint('soma', mech('lca'   , {"glcabar": 0.005   *params["glcabar_lca"]}))
cell.paint('soma', mech('cat'   , {"gcatbar": 0.000037*params["gcatbar_cat"]}))
cell.paint('soma', mech('gskch' , {"gskbar" : 0.001   *params["gskbar_gskch"]}))
cell.paint('soma', mech('cagk'  , {"gkbar"  : 0.0006  *params["gkbar_cagk"]}))
cell.paint('soma', mech('ccanl' , {"catau"  : 10      *params["catau_ccanl"], "caiinf": 5.0e-6*params["caiinf_ccanl"]}))
cell.paint('soma', cm=1.0*params["cm_mult"]/100.0, rL=410*params["ra_mult"])

cell.paint('granuleCellLayer', mech('ichan2', {"gnatbar":0.018*params["gnatbar_ichan2"],
                                               "gkfbar":0.004,
                                               "gksbar":0.006,
                                               "gl":0.00004*params["gl_ichan2"],
                                               "el":params["el_ichan2"]}))
cell.paint('granuleCellLayer', mech('nca'   , {"gncabar": 0.003   *params["gncabar_nca "]}))
cell.paint('granuleCellLayer', mech('lca'   , {"glcabar": 0.0075}))
cell.paint('granuleCellLayer', mech('cat'   , {"gcatbar": 0.000075}))
cell.paint('granuleCellLayer', mech('gskch' , {"gskbar" : 0.0004}))
cell.paint('granuleCellLayer', mech('cagk'  , {"gkbar"  : 0.0006  *params["gkbar_cagk"]}))
cell.paint('granuleCellLayer', mech('ccanl' , {"catau"  : 10      *params["catau_ccanl"], "caiinf": 5.0e-6*params["caiinf_ccanl"]}))
cell.paint('granuleCellLayer', cm=1.0*params["cm_mult"]/100.0, rL=410*params["ra_mult"])

cell.paint('innerThird', mech('ichan2', {"gnatbar":0.013*params["gnatbar_ichan2"],
                                         "gkfbar":0.004,
                                         "gksbar":0.006,
                                         "gl":0.000063*params["gl_ichan2"],
                                         "el":params["el_ichan2"]}))
cell.paint('innerThird', mech('nca'   , {"gncabar": 0.001   *params["gncabar_nca "]}))
cell.paint('innerThird', mech('lca'   , {"glcabar": 0.0075}))
cell.paint('innerThird', mech('cat'   , {"gcatbar": 0.00025}))
cell.paint('innerThird', mech('gskch' , {"gskbar" : 0.0002}))
cell.paint('innerThird', mech('cagk'  , {"gkbar"  : 0.001   *params["gkbar_cagk"]}))
cell.paint('innerThird', mech('ccanl' , {"catau"  : 10      *params["catau_ccanl"], "caiinf": 5.0e-6*params["caiinf_ccanl"]}))
cell.paint('innerThird', cm=1.6*params["cm_mult"]/100.0, rL=410*params["ra_mult"])

cell.paint('middleThird', mech('ichan2', {"gnatbar":0.008*params["gnatbar_ichan2"],
                                          "gkfbar":0.001,
                                          "gksbar":0.006,
                                          "gl":0.000063*params["gl_ichan2"],
                                          "el":params["el_ichan2"]}))
cell.paint('middleThird', mech('nca'   , {"gncabar": 0.001   *params["gncabar_nca "]}))
cell.paint('middleThird', mech('lca'   , {"glcabar": 0.0005}))
cell.paint('middleThird', mech('cat'   , {"gcatbar": 0.0005}))
cell.paint('middleThird', mech('gskch' , {"gskbar" : 0.0}))
cell.paint('middleThird', mech('cagk'  , {"gkbar"  : 0.0024  *params["gkbar_cagk"]}))
cell.paint('middleThird', mech('ccanl' , {"catau"  : 10      *params["catau_ccanl"], "caiinf": 5.0e-6*params["caiinf_ccanl"]}))
cell.paint('middleThird', cm=1.6*params["cm_mult"]/100.0, rL=410*params["ra_mult"])

cell.paint('outerThird', mech('ichan2', {"gnatbar":0.0*params["gnatbar_ichan2"],
                                         "gkfbar":0.001,
                                         "gksbar":0.008,
                                         "gl":0.000063*params["gl_ichan2"],
                                         "el":params["el_ichan2"]}))
cell.paint('outerThird', mech('nca'   , {"gncabar": 0.001   *params["gncabar_nca "]}))
cell.paint('outerThird', mech('lca'   , {"glcabar": 0.0}))
cell.paint('outerThird', mech('cat'   , {"gcatbar": 0.001}))
cell.paint('outerThird', mech('gskch' , {"gskbar" : 0.0}))
cell.paint('outerThird', mech('cagk'  , {"gkbar"  : 0.0024  *params["gkbar_cagk"]}))
cell.paint('outerThird', mech('ccanl' , {"catau"  : 10      *params["catau_ccanl"], "caiinf": 5.0e-6*params["caiinf_ccanl"]}))
cell.paint('outerThird', cm=1.6*params["cm_mult"]/100.0, rL=410*params["ra_mult"])

# add synapses
syn = mech('exp2syn', {"tau1": params["tau1_syn"], "tau2": params["tau2_syn"], "e": params["e_syn"]})
if params["syn_layer"] == "soma":
    cell.place('soma_syn', syn)
if params["syn_layer"] == "granuleCellLayer":
    cell.place('granule_syn', syn)
if params["syn_layer"] == "innerLayer":
    cell.place('inner_syn', syn)
if params["syn_layer"] == "middleLayer":
    cell.place('middle_syn', syn)
if params["syn_layer"] == "outerLayer":
    cell.place('outer_syn', syn)

cell.compartments_on_samples()

# make the model
model = arbor.single_cell_model(cell)

model.add_ion('nca', valence=2, int_con=1.0, ext_con=1.0, rev_pot=0);
model.add_ion('lca', valence=2, int_con=1.0, ext_con=1.0, rev_pot=params["elca"]);
model.add_ion('tca', valence=2, int_con=1.0, ext_con=1.0, rev_pot=params["etca"]);
model.add_ion('nat', valence=1, int_con=1.0, ext_con=1.0, rev_pot=params["enat"]);
model.add_ion('kf' , valence=1, int_con=1.0, ext_con=1.0, rev_pot=params["ekf"]);
model.add_ion('ks' , valence=1, int_con=1.0, ext_con=1.0, rev_pot=params["eks"]);
model.add_ion('sk' , valence=1, int_con=1.0, ext_con=1.0, rev_pot=params["esk"]);

model.probe('voltage', where=loc(0,0), frequency=50000)

model.run(tsim)


# Plot the recorded voltages over time.
fig, ax = plt.subplots()
for t in model.traces:
    ax.plot(t.time, t.value)

legend_labels = ['{}: {}'.format(s.variable, s.location) for s in model.traces]
ax.legend(legend_labels)
ax.set(xlabel='time (ms)', ylabel='voltage (mV)', title='Purkinje cell demo')
plt.xlim(0, tsim)
#plt.ylim(-80,50)
ax.grid()

if plot_to_file:
    fig.savefig("voltages.png", dpi=300)
else:
    plt.show()