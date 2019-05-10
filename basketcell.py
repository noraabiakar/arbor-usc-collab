# Begin Python script
from neuron import h

'''
Single-compartment implementation of a dentate basket cell.  Values for biophysics taken from 
Santhakumar et al 2004.
'''

# The following layers are dictionaries whose keys are NEURON sections.
# Each entry of the dictionaries contains a list. The list contains
# the normalized distances along the section that lie within a layer.
def makeSecDict():
	SecList = {}
	SecList['soma'] = {}
	
	return SecList

def makeLayerDict(cell):
	LayerDict = {}
	LayerDict['Apical'] = {}
	LayerDict['Apical']['soma'] = cell.soma
	
	return LayerDict

# Since all the morphology is defined by HOC code, we need a pointer to the HOC object.
def loadMorph(morphFileName):
	param = {}
	param['c'] = mksoma()
	return param

class mksoma():
	def __init__(self):
		self.soma = h.Section()
		h.pt3dclear(sec=self.soma)
		h.pt3dadd(0,0,-10,15,sec=self.soma)
		h.pt3dadd(0,0,10,15,sec=self.soma)

# Initilize the list of synapses for the cell
def makeSynGroups(cell):
	SynGroups = {}
	SynGroups['AMPA'] = {}
	SynGroups['AMPA']['soma'] = []
	
	SynGroups['NMDA'] = {}
	SynGroups['NMDA']['soma'] = []
	
	SynGroups['GABA'] = {}
	SynGroups['GABA']['soma'] = []
	
	return SynGroups

# Defines the major axis of the morphology
def getNewAxis():
	new_axis = {}
	new_axis['new_axis'] = [0,0,0]
	return new_axis

# Function to return the nseg resolution
def getNsegRes():
	return 50

# Function to return soma
def getSoma(cell):
	return cell.c.soma

# Function to return the "center" of the morphology
# The reference point is set to the somatic location
def getCenter(soma):
	center = (0,0,0)
	return center

# Function to return to a dendrite lists organized by type (apical or basal)
def getDendTypeList(cell):
	dendTypeList = {}
	dendTypeList['Apical'] = [cell.c.soma]
	
	return dendTypeList

# Function to return apical dendrites
def getApicDend(cell):
	return cell.c.dend

# Function to return bounds of the layers
def getBounds(maxExtent):
	bounds = {}
	bounds['Apical'] = {}
	bounds['Apical']['soma'] = (0,0)
	
	return bounds

# Function to make the lists containing the locations of the segments
# The list is organized as [ [x], [y], [z] ]
def makeSegLocDict(cell):
	SegLocDict = {}
	SegLocDict['Apical'] = {}
	SegLocDict['Apical']['soma'] = [ [], [], [] ]
	
	return SegLocDict

# Function to specify the biophysics of the cell
def getBiophysics(cell):
	cell.c.soma.nseg = 1
	cell.c.soma.L = 20
	cell.c.soma.diam = 15
	
	cell.c.soma.insert('ccanl')
	cell.c.soma.insert('borgka')
	cell.c.soma.insert('nca')
	cell.c.soma.insert('lca')
	cell.c.soma.insert('gskch')
	cell.c.soma.insert('cagk')
	
	cell.c.soma.catau_ccanl = 10
	cell.c.soma.caiinf_ccanl = 5.0e-6
	cell.c.soma.gkabar_borgka = 0.00015
	cell.c.soma.gncabar_nca = 0.0008
	cell.c.soma.glcabar_lca = 0.005
	cell.c.soma.gskbar_gskch = 0.000002
	cell.c.soma.gkbar_cagk = 0.0002

	cell.c.soma.insert('ichan2')
	cell.c.soma.cm = 1.4
	cell.c.soma.gnatbar_ichan2 = 0.12
	cell.c.soma.gkfbar_ichan2 = 0.013
	cell.c.soma.gl_ichan2 = 0.00018

	cell.c.soma.Ra = 100
	cell.c.soma.enat = 55
	cell.c.soma.ekf = -90
	cell.c.soma.ek = -90
	cell.c.soma.elca = 130
	cell.c.soma.esk = -90
	cell.c.soma.el_ichan2 = -65
	cell.c.soma.cao = 2

# Function to specify the biophysics of the reduced cell model
def getReducedBiophysics(cell):
	cell.soma.nseg = 1
	cell.soma.L = 20
	cell.soma.diam = 15
	
	cell.soma.insert('ccanl')
	cell.soma.insert('borgka')
	cell.soma.insert('nca')
	cell.soma.insert('lca')
	cell.soma.insert('gskch')
	cell.soma.insert('cagk')
	
	cell.soma.catau_ccanl = 10
	cell.soma.caiinf_ccanl = 5.0e-6
	cell.soma.gkabar_borgka = 0.00015
	cell.soma.gncabar_nca = 0.0008
	cell.soma.glcabar_lca = 0.005
	cell.soma.gskbar_gskch = 0.000002
	cell.soma.gkbar_cagk = 0.0002

	cell.soma.insert('ichan2')
	cell.soma.cm = 1.4
	cell.soma.gnatbar_ichan2 = 0.12
	cell.soma.gkfbar_ichan2 = 0.013
	cell.soma.gl_ichan2 = 0.00018

	cell.soma.Ra = 100
	cell.soma.enat = 55
	cell.soma.ekf = -90
	cell.soma.ek = -90
	cell.soma.elca = 130
	cell.soma.esk = -90
	cell.soma.el_ichan2 = -65
	cell.soma.cao = 2

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
	for syntype in cell.synGroups:	
		# soma
		syn = createSyn(cell.synvars,cell.soma,0.5)
		cell.synGroups[syntype]['soma'].append(syn)

# End file
