# Granule cell class
from neuron import h
import numpy as np
import pickle
h.load_file("hoc_files/importCell.hoc")

# The following layers are dictionaries whose keys are NEURON sections.
# Each entry of the dictionaries contains a list. The list contains
# the normalized distances along the section that lie within a layer.
def makeSecDict():
    SecList = {}
    SecList['soma'] = {}
    SecList['granuleCellLayer'] = {}
    SecList['innerThird'] = {}
    SecList['middleThird'] = {}
    SecList['outerThird'] = {}
    return SecList

def makeLayerDict(cell):
    LayerDict = {}
    LayerDict['Apical'] = {}
    LayerDict['Apical']['soma'] = cell.soma
    LayerDict['Apical']['granuleCellLayer'] = cell.granuleCellLayer
    LayerDict['Apical']['innerThird'] = cell.innerThird
    LayerDict['Apical']['middleThird'] = cell.middleThird
    LayerDict['Apical']['outerThird'] = cell.outerThird

    return LayerDict

# Since all the morphology is defined by HOC code, we need a pointer to the HOC object.
def loadMorph(morphFileName):
    param = {}
    param['c'] = h.mkcell(morphFileName)
    return param

# Initilize the list of synapses for the cell
def makeSynGroups(cell):
    SynGroups = {}
    SynGroups['AMPA'] = {}
    SynGroups['AMPA']['soma'] = []
    SynGroups['AMPA']['granuleCellLayer'] = []
    SynGroups['AMPA']['innerThird'] = []
    SynGroups['AMPA']['middleThird'] = []
    SynGroups['AMPA']['outerThird'] = []

    # SynGroups['NMDA'] = {}
    # SynGroups['NMDA']['soma'] = []
    # SynGroups['NMDA']['granuleCellLayer'] = []
    # SynGroups['NMDA']['innerThird'] = []
    # SynGroups['NMDA']['middleThird'] = []
    # SynGroups['NMDA']['outerThird'] = []
    #
    # SynGroups['GABA'] = {}
    # SynGroups['GABA']['soma'] = []
    # SynGroups['GABA']['granuleCellLayer'] = []
    # SynGroups['GABA']['innerThird'] = []
    # SynGroups['GABA']['middleThird'] = []
    # SynGroups['GABA']['outerThird'] = []

    return SynGroups

# Defines the major axis of the morphology
def getNewAxis():
    new_axis = {}
    new_axis['new_axis'] = [np.cos(np.pi/2),0,np.sin(np.pi/2)]
    return new_axis

# Function to return soma
def getSoma(cell):
    if cell.modeltype == 'Multi':
        return cell.c.soma[0]
    else:
        return cell.soma

# Function to return the "center" of the morphology
# The reference point is set to the somatic location
def getCenter(soma):
    soma.push()
    center = np.array((h.x3d(0),h.y3d(0),h.z3d(0)))
    h.pop_section()
    return center

# Function to return to a dendrite lists organized by type (apical or basal)
def getDendTypeList(cell):
    dendTypeList = {}
    dendTypeList['Apical'] = getApicDend(cell)

    return dendTypeList

# Function to return apical dendrites
def getApicDend(cell):
    return cell.c.dend

# Function to return bounds of the layers
def getBounds(maxExtent):
    bounds = {}
    bounds['Apical'] = {}
    bounds['Apical']['soma'] = (0,0)
    bounds['Apical']['granuleCellLayer'] = (0,0.1*maxExtent['Apical'])
    bounds['Apical']['innerThird'] = (0.1*maxExtent['Apical'],0.3*maxExtent['Apical'])
    bounds['Apical']['middleThird'] = (0.3*maxExtent['Apical'],0.6*maxExtent['Apical'])
    bounds['Apical']['outerThird'] = (0.6*maxExtent['Apical'],maxExtent['Apical'])

    return bounds

# Function to make the lists containing the locations of the segments
# The list is organized as [ [x], [y], [z] ]
def makeSegLocDict(cell):
    SegLocDict = {}
    SegLocDict['Apical'] = {}
    SegLocDict['Apical']['soma'] = [ [], [], [] ]
    SegLocDict['Apical']['granuleCellLayer'] = [ [], [], [] ]
    SegLocDict['Apical']['innerThird'] = [ [], [], [] ]
    SegLocDict['Apical']['middleThird'] = [ [], [], [] ]
    SegLocDict['Apical']['outerThird'] = [ [], [], [] ]
    return SegLocDict

# Function to specify the biophysics of the cell
def getBiophysics(cell, in_param):
    #cell.c.soma[0].L *= 2

    # Now, insert the proper biophysics for each section.
    for sec in cell.c.all:
        # sec.insert('hh') #**
        # sec.gnabar_hh= in_param["hh_gnabar"] #**
        # sec.gkbar_hh = in_param["hh_gkbar"] #**
        # sec.gl_hh = in_param["hh_gl"] #**
        # sec.ena = in_param["hh_ena"] #**
        # sec.ek = in_param["hh_ek"] #**
        # sec.Ra = in_param["ra"] #**
        # sec.cm = in_param["cm"] #**

        # sec.insert('ccanl')
        # sec.catau_ccanl=10*cell.k1
        # sec.caiinf_ccanl=5.0e-6*cell.l1
        sec.Ra = 410 * in_param["ra_mult"]

    for sec in cell.c.somatic:
        sec.insert('ichan2')
        sec.gnatbar_ichan2 = 0.120 * in_param["gnatbar_ichan2"]
        sec.gkfbar_ichan2  = 0.016 * in_param["gkfbar_ichan2"]
        sec.gksbar_ichan2  = 0.006 * in_param["gksbar_ichan2"]
        sec.gl_ichan2    = 0.00004 * in_param["gl_ichan2"]
        sec.el_ichan2     =          in_param["el_ichan2"]

        sec.insert('borgka')
        sec.gkabar_borgka  = 0.001 * in_param["gkabar_borgka"]

        sec.insert('nca')
        sec.gncabar_nca    = 0.001 * in_param["gncabar_nca"]

        sec.insert('lca')
        sec.glcabar_lca    = 0.005 * in_param["glcabar_lca"]

        sec.insert('cat')
        sec.gcatbar_cat = 0.000037 * in_param["gcatbar_cat"]

        sec.insert('gskch')
        sec.gskbar_gskch   = 0.001 * in_param["gskbar_gskch"]

        sec.insert('cagk')
        sec.gkbar_cagk    = 0.0006 * in_param["gkbar_cagk"]

        sec.cm               = 1.0 * in_param["cm_mult"]

        sec.cao       = in_param["cao"]
        sec.ek        = in_param["ek"]
        sec.enat      = in_param["enat"]
        sec.ekf       = in_param["ekf"]
        sec.eks       = in_param["eks"]
        sec.elca      = in_param["elca"]
        sec.etca      = in_param["etca"]
        sec.esk       = in_param["esk"]

    # for sec in cell.c.dend:
    #     sec.insert('ichan2')
    #     sec.insert('nca')
    #     sec.insert('lca')
    #     sec.insert('cat')
    #     sec.insert('gskch')
    #     sec.insert('cagk')

    for sec in cell.granuleCellLayer:
        if len(cell.granuleCellLayer[sec]) > 0:
            sec.insert('ichan2')
            sec.insert('nca')
            sec.insert('lca')
            sec.insert('cat')
            sec.insert('gskch')
            sec.insert('cagk')

            #for norm_dist in cell.granuleCellLayer[sec]:
            sec.gnatbar_ichan2 = 0.018   * in_param["gnatbar_ichan2"]
            sec.gkfbar_ichan2  = 0.004
            sec.gksbar_ichan2  = 0.006
            sec.gl_ichan2      = 0.00004 * in_param["gl_ichan2"]
            sec.el_ichan2      =           in_param["el_ichan2"]
            sec.gncabar_nca    = 0.003   * in_param["gncabar_nca"]
            sec.glcabar_lca    = 0.0075
            sec.gcatbar_cat    = 0.000075
            sec.gskbar_gskch   = 0.0004
            sec.gkbar_cagk     = 0.0006  * in_param["gkbar_cagk"]
            sec.cm             = 1.0     * in_param["cm_mult"]

            sec.cao       = in_param["cao"]
            sec.enat      = in_param["enat"]
            sec.ekf       = in_param["ekf"]
            sec.eks       = in_param["eks"]
            sec.elca      = in_param["elca"]
            sec.etca      = in_param["etca"]
            sec.esk       = in_param["esk"]


    for sec in cell.innerThird:
        if len(cell.innerThird[sec]) > 0:
            sec.insert('ichan2')
            sec.insert('nca')
            sec.insert('lca')
            sec.insert('cat')
            sec.insert('gskch')
            sec.insert('cagk')

            # for norm_dist in cell.innerThird[sec]:
            sec.gnatbar_ichan2 = 0.013    * in_param["gnatbar_ichan2"]
            sec.gkfbar_ichan2  = 0.004
            sec.gksbar_ichan2  = 0.006
            sec.gl_ichan2      = 0.000063 * in_param["gl_ichan2"]
            sec.el_ichan2      =            in_param["el_ichan2"]
            sec.gncabar_nca    = 0.001    * in_param["gncabar_nca"]
            sec.glcabar_lca    = 0.0075
            sec.gcatbar_cat    = 0.00025
            sec.gskbar_gskch   = 0.0002
            sec.gkbar_cagk     = 0.001    * in_param["gkbar_cagk"]
            sec.cm             = 1.6      * in_param["cm_mult"]

            sec.cao       = in_param["cao"]
            sec.enat      = in_param["enat"]
            sec.ekf       = in_param["ekf"]
            sec.eks       = in_param["eks"]
            sec.elca      = in_param["elca"]
            sec.etca      = in_param["etca"]
            sec.esk       = in_param["esk"]


    for sec in cell.middleThird:
        if len(cell.middleThird[sec]) > 0:
            sec.insert('ichan2')
            sec.insert('nca')
            sec.insert('lca')
            sec.insert('cat')
            sec.insert('gskch')
            sec.insert('cagk')

            # for norm_dist in cell.middleThird[sec]:
            sec.gnatbar_ichan2 = 0.008    * in_param["gnatbar_ichan2"]
            sec.gkfbar_ichan2  = 0.001
            sec.gksbar_ichan2  = 0.006
            sec.gl_ichan2      = 0.000063 * in_param["gl_ichan2"]
            sec.el_ichan2      =            in_param["el_ichan2"]
            sec.gncabar_nca    = 0.001    * in_param["gncabar_nca"]
            sec.glcabar_lca    = 0.0005
            sec.gcatbar_cat    = 0.0005
            sec.gskbar_gskch   = 0.0
            sec.gkbar_cagk     = 0.0024   * in_param["gkbar_cagk"]
            sec.cm             = 1.6      * in_param["cm_mult"]

            sec.cao       = in_param["cao"]
            sec.enat      = in_param["enat"]
            sec.ekf       = in_param["ekf"]
            sec.eks       = in_param["eks"]
            sec.elca      = in_param["elca"]
            sec.etca      = in_param["etca"]
            sec.esk       = in_param["esk"]


    for sec in cell.outerThird:
        if len(cell.outerThird[sec]) > 0:
            sec.insert('ichan2')
            sec.insert('nca')
            sec.insert('lca')
            sec.insert('cat')
            sec.insert('gskch')
            sec.insert('cagk')

            # for norm_dist in cell.outerThird[sec]:
            sec.gnatbar_ichan2 = 0.0      * in_param["gnatbar_ichan2"]
            sec.gkfbar_ichan2  = 0.001
            sec.gksbar_ichan2  = 0.008
            sec.gl_ichan2      = 0.000063 * in_param["gl_ichan2"]
            sec.el_ichan2      =            in_param["el_ichan2"]
            sec.gncabar_nca    = 0.001    * in_param["gncabar_nca"]
            sec.glcabar_lca    = 0.0
            sec.gcatbar_cat    = 0.001
            sec.gskbar_gskch   = 0.0
            sec.gkbar_cagk     = 0.0024   * in_param["gkbar_cagk"]
            sec.cm             = 1.6      * in_param["cm_mult"]

            sec.cao       = in_param["cao"]
            sec.enat      = in_param["enat"]
            sec.ekf       = in_param["ekf"]
            sec.eks       = in_param["eks"]
            sec.elca      = in_param["elca"]
            sec.etca      = in_param["etca"]
            sec.esk       = in_param["esk"]

    # for sec in cell.c.all:
    #     sec.cao       = in_param["cao"]
    #     sec.ek        = in_param["ek"]
    #     sec.enat      = in_param["enat"]
    #     sec.ekf       = in_param["ekf"]
    #     sec.eks       = in_param["eks"]
    #     sec.elca      = in_param["elca"]
    #     sec.etca      = in_param["etca"]
    #     sec.esk       = in_param["esk"]

# Function to specify the biophysics of the reduced cell model
def getReducedBiophysics(cell):
    with open('granulecell_marascoProp.pickle') as f:
        props = pickle.load(f)

    cell.soma.nseg = 1
    cell.soma.L = 23.319360733032227
    cell.soma.diam = 11.659700393676758
    cell.soma.Ra = 410
    cell.soma.cm = 9.8

    secDict = {}
    secDict['granuleCellLayer'] = cell.granuleCellLayer
    secDict['innerThird'] = cell.innerThird
    secDict['middleThird'] = cell.middleThird
    secDict['outerThird'] = cell.outerThird

    for layer in props:
        secDict[layer].nseg = props[layer]['nseg']
        secDict[layer].L = props[layer]['L']
        secDict[layer].diam = props[layer]['d']
        secDict[layer].Ra = props[layer]['Ra']

    cell.CmMult = 9.8
    cell.a1 = 7.0				#gnatbar_ichan2
    cell.b1 = 2.25				#gkfbar_ichan2
    cell.c1 = 1.0				#gksbar_ichan2
    cell.d1 = 9.0				#gkabar_borgka
    cell.e1 = 1/1.36			#gncabar_nca
    cell.f1 = 0.5				#glcabar_lca
    cell.g1 = 2.0				#gcatbar_cat
    cell.h1 = 1.0				#gskbar_gskch
    cell.i1 = 1/5.0			#gkbar_cagk
    cell.j1 = 7.2538			#gl_ichan2
    cell.k1 = 1.0				#catau_ccanl
    cell.l1 = 1.0				#caiinf_ccanl

    cell.soma.insert('ccanl')
    cell.soma.catau_ccanl=10*cell.k1
    cell.soma.caiinf_ccanl=5.0e-6*cell.l1
    cell.soma.insert('ichan2')
    cell.soma.gnatbar_ichan2 = 0.12*cell.a1
    cell.soma.gkfbar_ichan2=0.016*cell.b1
    cell.soma.gksbar_ichan2=0.006*cell.c1
    cell.soma.insert('borgka')
    cell.soma.gkabar_borgka=0.001*cell.d1
    cell.soma.insert('nca')
    cell.soma.gncabar_nca=0.001*cell.e1
    cell.soma.insert('lca')
    cell.soma.glcabar_lca=0.005*cell.f1
    cell.soma.insert('cat')
    cell.soma.gcatbar_cat=0.000037*cell.g1
    cell.soma.insert('gskch')
    cell.soma.gskbar_gskch=0.001*cell.h1*5
    cell.soma.insert('cagk')
    cell.soma.gkbar_cagk=0.0006*cell.i1*5
    cell.soma.gl_ichan2=0.00004*cell.j1

    for layer in secDict:
        secDict[layer].insert('ccanl')
        secDict[layer].catau_ccanl=10*cell.k1
        secDict[layer].caiinf_ccanl=5.0e-6*cell.l1
        secDict[layer].insert('ichan2')
        secDict[layer].insert('nca')
        secDict[layer].insert('lca')
        secDict[layer].insert('cat')
        secDict[layer].insert('gskch')
        secDict[layer].insert('cagk')

    cell.granuleCellLayer.gnatbar_ichan2 = 0.018*cell.a1*props['granuleCellLayer']['fact']
    cell.granuleCellLayer.gkfbar_ichan2=0.004*props['granuleCellLayer']['fact']
    cell.granuleCellLayer.gksbar_ichan2=0.006*props['granuleCellLayer']['fact']
    cell.granuleCellLayer.gncabar_nca=0.003*cell.e1*props['granuleCellLayer']['fact']
    cell.granuleCellLayer.glcabar_lca=0.0075*props['granuleCellLayer']['fact']
    cell.granuleCellLayer.gcatbar_cat=0.000075*props['granuleCellLayer']['fact']
    cell.granuleCellLayer.gskbar_gskch=0.0004*props['granuleCellLayer']['fact']
    cell.granuleCellLayer.gkbar_cagk=0.0006*cell.i1*props['granuleCellLayer']['fact']
    cell.granuleCellLayer.gl_ichan2=0.00004*cell.j1*props['granuleCellLayer']['fact']
    cell.granuleCellLayer.cm=1.0*cell.CmMult*props['granuleCellLayer']['fact']

    cell.innerThird.gnatbar_ichan2 = 0.013*cell.a1*props['innerThird']['fact']
    cell.innerThird.gkfbar_ichan2=0.004*props['innerThird']['fact']
    cell.innerThird.gksbar_ichan2=0.006*props['innerThird']['fact']
    cell.innerThird.gncabar_nca=0.001*cell.e1*props['innerThird']['fact']
    cell.innerThird.glcabar_lca=0.0075*props['innerThird']['fact']
    cell.innerThird.gcatbar_cat=0.00025*props['innerThird']['fact']
    cell.innerThird.gskbar_gskch=0.0002*props['innerThird']['fact']
    cell.innerThird.gkbar_cagk=0.001*cell.i1*props['innerThird']['fact']
    cell.innerThird.gl_ichan2=0.000063*cell.j1*props['innerThird']['fact']
    cell.innerThird.cm=1.6*cell.CmMult*props['innerThird']['fact']

    cell.middleThird.gnatbar_ichan2 = 0.008*cell.a1*props['middleThird']['fact']
    cell.middleThird.gkfbar_ichan2=0.001*props['middleThird']['fact']
    cell.middleThird.gksbar_ichan2=0.006*props['middleThird']['fact']
    cell.middleThird.gncabar_nca=0.001*cell.e1*props['middleThird']['fact']
    cell.middleThird.glcabar_lca=0.0005*props['middleThird']['fact']
    cell.middleThird.gcatbar_cat=0.0005*props['middleThird']['fact']
    cell.middleThird.gskbar_gskch=0.0*props['middleThird']['fact']
    cell.middleThird.gkbar_cagk=0.0024*cell.i1*props['middleThird']['fact']
    cell.middleThird.gl_ichan2=0.000063*cell.j1*props['middleThird']['fact']
    cell.middleThird.cm=1.6*cell.CmMult*props['middleThird']['fact']

    cell.outerThird.gnatbar_ichan2 = 0.0*cell.a1*props['outerThird']['fact']
    cell.outerThird.gkfbar_ichan2=0.001*props['outerThird']['fact']
    cell.outerThird.gksbar_ichan2=0.008*props['outerThird']['fact']
    cell.outerThird.gncabar_nca=0.001*cell.e1*props['outerThird']['fact']
    cell.outerThird.glcabar_lca=0.0*props['outerThird']['fact']
    cell.outerThird.gcatbar_cat=0.001*props['outerThird']['fact']
    cell.outerThird.gskbar_gskch=0.0*props['outerThird']['fact']
    cell.outerThird.gkbar_cagk=0.0024*cell.i1*props['outerThird']['fact']
    cell.outerThird.gl_ichan2=0.000063*cell.j1*props['outerThird']['fact']
    cell.outerThird.cm=1.6*cell.CmMult*props['outerThird']['fact']

    cell.soma.enat = 45
    cell.soma.ekf = -90
    cell.soma.eks = -90
    cell.soma.ek = -90
    cell.soma.elca = 130
    cell.soma.etca = 130
    cell.soma.esk = -90
    cell.soma.el_ichan2 = -73
    cell.soma.cao = 2
    for layer in secDict:
        secDict[layer].enat = 45
        secDict[layer].ekf = -90
        secDict[layer].eks = -90
        secDict[layer].ek = -90
        secDict[layer].elca = 130
        secDict[layer].etca = 130
        secDict[layer].esk = -90
        secDict[layer].el_ichan2 = -73
        secDict[layer].cao = 2

# Function to create a synapse at the chosen segment in a section
def createSyn(synvars,sec_choice,seg_choice):
    if synvars['type'] == "E2-NMDA2":
        syn = h.Exp2Syn(sec_choice(seg_choice))
        nmda = h.Exp2NMDA_Wang(sec_choice(seg_choice))
        nmda_flag = 1
    if synvars['type'] == "E2":
        syn = h.Exp2Syn(sec_choice(seg_choice))
    if synvars['type'] == "E2_Prob":
        syn = h.E2_Prob(sec_choice(seg_choice))
        syn.P = synvars['P']
    if synvars['type'] == "E2_STP_Prob":
        syn = h.E2_STP_Prob(sec_choice(seg_choice))
    if synvars['type'] == "STDPE2":
        syn = h.STDPE2(sec_choice(seg_choice))
    if synvars['type'] == "STDPE2_Clo":
        syn = h.STDPE2_Clo(sec_choice(seg_choice))
    if synvars['type'] == "STDPE2_STP"	:
        syn = h.STDPE2_STP(sec_choice(seg_choice))
    if synvars['type'] == "STDPE2_Prob":
        syn = h.STDPE2_Prob(sec_choice(seg_choice))
        syn.P = synvars['P']
    #initializes different variables depending on synapse
    if (synvars['type'] == "STDPE2_STP")|(synvars['type'] == "E2_STP_Prob"):
        syn.F1 = synvars['F1']
    if  (synvars['type'] == "STDPE2_Clo" )|( synvars['type'] == "STDPE2_STP")|( synvars['type'] == "STDPE2")| (synvars['type'] == "STDPE2_Prob"):
        syn.wmax = synvars['wmax']
        syn.wmin = synvars['wmin']
        syn.thresh = synvars['thresh']
    if  (synvars['type'] == "E2_Prob" )|( synvars['type'] == "E2_STP_Prob")|(synvars['type'] == "STDPE2_STP") | (synvars['type'] == "STDPE2_Prob"):
        h.use_mcell_ran4(1)
        syn.seed = self.ranGen.randint(1,4.295e9)
    syn.tau1 = 0.5
    syn.tau2 = 0.6
    syn.e = 0

    return syn

# Function to add synapses to the reduced cell model
def addReducedSynapses(cell):
    dendLayers = {}
    dendLayers['granuleCellLayer'] = cell.granuleCellLayer
    dendLayers['innerThird'] = cell.innerThird
    dendLayers['middleThird'] = cell.middleThird
    dendLayers['outerThird'] = cell.outerThird
    for syntype in cell.synGroups:
        # soma
        syn = createSyn(cell.synvars,cell.soma,0.5)
        cell.synGroups[syntype]['soma'].append(syn)

        for layer in dendLayers:
            seg_pos = np.linspace(0.5/dendLayers[layer].nseg,1-0.5/dendLayers[layer].nseg,dendLayers[layer].nseg)
            for ii in range(dendLayers[layer].nseg):
                syn = createSyn(cell.synvars,dendLayers[layer],seg_pos[ii])
                cell.synGroups[syntype][layer].append(syn)

# End of file
