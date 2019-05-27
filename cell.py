from neuron import h
pc = h.ParallelContext()
import random
import numpy as np

class Cell:
    def __init__(self,cellID,locn,synvars,celltype,morphFileName,modeltype='Multi'):
        # Sets the synapse properties
        self.synvars = synvars

        # Sets ID of cell
        self.ID = cellID

        # Sets location of cell
        self.location = locn

        # Create random number generator for cell
        self.ranGen = random.Random()
        self.ranGen.seed(cellID)

        # Initialize the list of NetCons for the cell
        self.nc_list = []

        # Set cell type by importing correct module
        # Based on cell type module, the called functions may perform slightly
        # different actions, but it is more related to the names of the layers
        # Make sure each module has the same function names
        self.celltype = __import__(celltype)
        self.type = celltype
        self.modeltype = modeltype

        # Defining nseg using a distance-rule to try to standardize the segment
        # lengths across sections
        self.nseg_res = self.celltype.getNsegRes()

        # Create synapse lists
        self.synGroups = self.celltype.makeSynGroups(self)

        # If the modeltype flag is Multi, instantiate a multicompartmental model
        # Else, make a single compartment model...

        if modeltype == 'Multi':
            # Load in morphology
            self.morphFileName = morphFileName
            self.setAttr(self.celltype.loadMorph(morphFileName))

            # Create section dictionaries
            self.setAttr(self.celltype.makeSecDict())

            # Need to create a list that has all of the section dictionaries
            # for the layer organization section
            #self.layerDict = self.makeAttrDict(self.celltype.makeSecDict())
            self.layerDict = self.celltype.makeLayerDict(self)

            ##############################
            # Organize sections by layer #
            ##############################
            # Setting new axis
            self.setAttr(self.celltype.getNewAxis())

            # Get somatic location
            soma = self.celltype.getSoma(self)
            self.center = self.celltype.getCenter(soma)

            # Create a data structure to hold locations of segments
            self.segLocDict = self.celltype.makeSegLocDict(self)

            # Run for each dendrite type (apical or basal)
            self.dendTypeList = self.celltype.getDendTypeList(self)

            # Create the dictionary with keys being distance from final swc coordinate
            # to soma and the values being pointer to section
            self.dist2sec = self.make_sec_dict(self.dendTypeList)

            # Add soma to dist2sec
            self.dist2sec['Apical'][0.0] = soma

            self.distKeys = {}
            for dendType in self.layerDict:
                self.distKeys[dendType] = sorted(self.dist2sec[dendType].keys())

            # Find maximum extent
            self.maxExtent = self.getMaxExtent(self.dendTypeList)

            # Set nsegs for each section
            self.norm_seg(self.nseg_res)
            #self.d_lambda()

            # Set boundaries for layers
            self.bounds = self.celltype.getBounds(self.maxExtent)

            # Separating segments of sections into different layers
            self.orgSeg(self.dendTypeList,self.layerDict,self.bounds,self.segLocDict)

            # Creating the data structure for the soma, which doesn't need to be organized
            # into layers, but this is necessary for placing synapses
            self.soma[soma] = [0.5]
            for layer in self.layerDict['Apical']:
                if layer == 'soma':
                    self.layerDict['Apical'][layer][soma] = [0.5]
                else:
                    self.layerDict['Apical'][layer][soma] = []

            for ii in range(3):
                self.segLocDict['Apical']['soma'][ii].append(self.center[ii])

            # Create synapses
            all = []
            all.append(self.celltype.getSoma(self))
            for d in self.dendTypeList['Apical']:
                all.append(d)
            
            for dendType in self.layerDict:
                for synType in self.synGroups:
                    self.addSegmentSynapse(self.layerDict[dendType],all,self.synGroups[synType],synType)

            # Add biophysics
            self.celltype.getBiophysics(self)

        elif modeltype == 'Single':
            # Not the best implementation, I feel, because the soma compartment won't be
            # stored as self.c.soma.
            self.soma = h.Section()
            if 'granule' in celltype:
                self.granuleCellLayer = h.Section()
                self.innerThird = h.Section()
                self.middleThird = h.Section()
                self.outerThird = h.Section()

                self.granuleCellLayer.connect(self.soma,1,0)
                self.innerThird.connect(self.granuleCellLayer,1,0)
                self.middleThird.connect(self.innerThird,1,0)
                self.outerThird.connect(self.middleThird,1,0)

            if 'ca3' in celltype:
                self.axon = h.Section()
                self.axon.connect(self.soma,1,0)

                self.oriensProximal = h.Section()
                self.oriensProximal.connect(self.soma,1,0)

                self.oriensDistal = h.Section()
                self.oriensDistal.connect(self.oriensProximal,1,0)

                self.lucidum = h.Section()
                self.lucidum.connect(self.soma,1,0)

                self.radiatum = h.Section()
                self.radiatum.connect(self.lucidum,1,0)

                self.lacunosumMEC = h.Section()
                self.lacunosumMEC.connect(self.radiatum,1,0)

                self.lacunosumLEC = h.Section()
                self.lacunosumLEC.connect(self.lacunosumMEC,1,0)

            # Each cell class needs to have this function to define their own conductance
            # levels.
            self.celltype.getReducedBiophysics(self)

            # Each cell class needs to have a single compartment version to add all the
            # different synapses
            self.celltype.addReducedSynapses(self)

        else:
            print 'Please choose \'Single\' or \'Multi\' for the modeltype option'

    '''#########################################################
    # 1. Group of functions that set/get attributes to a class #
    #########################################################'''

    ###########################################
    # Function to set attributes to the class #
    ###########################################
    def setAttr(self,dictionary):
        for k,v in dictionary.iteritems():
            setattr(self,str(k),v)

    ###########################################
    # Function to get attribute from the class #
    ###########################################
    def getAttr(self,attribute):
        return getattr(self,attribute),synType

    #############################################
    # Function to append attributes into a list #
    #############################################
    def makeAttrDict(self,attributes):
        attrDict = {}
        for attribute in attributes:
            attrDict[attribute] = self.getAttr(attribute)

        return attrDict

    '''############################################
    # 2. Code for organizing sections into layers #
    ############################################'''

    ################################################################
    # Function that creates a dictionary for the sections based on #
    # distance of section's final swc coordinate from soma         #
    ################################################################
    def make_sec_dict(self,dendTypeList):
        d = {}
        for dendType in dendTypeList:
            d[dendType] = {}
            dend = dendTypeList[dendType]
            for ii in range(len(dend)):
                dend[ii].push()
                idx = h.n3d()-1
                dist = np.sum((np.array([h.x3d(idx),h.y3d(idx),h.z3d(idx)])-self.center)**2)**0.5
                d[dendType][dist] = dend[ii]
                h.pop_section()

        return d


    ###########################################################
    # Function that assigns a number of segments to a section #
    # based on the specified segment resolution               #
    ###########################################################
    def norm_seg(self,nseg_res):
        for dendType in self.dendTypeList:
            for sec in self.dendTypeList[dendType]:
                if sec.L >= nseg_res:
                    sec.nseg = int(np.ceil(sec.L/nseg_res))
                else:
                    sec.nseg = 1

    ########################################################################
    # Function that assigns number of segments based on                    #
    # d_lambda rule                                                        #
    # http://www.neuron.yale.edu/neuron/static/docs/d_lambda/d_lambda.html #
    ########################################################################
    def d_lambda(self):
        freq = 100
        d_lambda = 0.1
        for dendType in self.dendTypeList:
            for sec in self. dendTypeList[dendType]:
                sec.push()
                if h.n3d() < 2:
                    lambda_f = 1e5*np.sqrt(sec.diam/(4*np.pi*freq*sec.Ra*sec.cm))
                else:
                    x1 = h.arc3d(0)
                    d1 = h.diam3d(0)
                    lam = 0
                    for ii in range(1,int(h.n3d())):
                        x2 = h.arc3d(ii)
                        d2 = h.diam3d(ii)
                        lam += (x2-x1)/np.sqrt(d1+d2)
                        x1 = x2
                        d1 = d2
                    lam *= np.sqrt(2)*1e-5*np.sqrt(4*np.pi*freq*sec.Ra*sec.cm)
                    lambda_f = sec.L/lam

                nseg = (int((sec.L/(d_lambda*lambda_f)+0.9)/2)*2 + 1)/2
                if nseg%2 == 0:
                     nseg += 1

                sec.nseg = nseg
                h.pop_section()

    ###################################################################
    # Function that calculates the maximum distance that a morphology #
    # extends from the soma along its major axis                      #
    ###################################################################
    def getMaxExtent(self,dendTypeList):
        maxExtent = {}
        for dendType in dendTypeList:
            dend = dendTypeList[dendType]
            values = []
            for i in range(len(dend)):
                dend[i].push()
                for j in range(int(h.n3d())):
                    swc = np.array([h.x3d(j),h.y3d(j),h.z3d(j)])
                    distance = abs(np.dot(self.new_axis,swc-self.center))
                    values.append(distance)
                h.pop_section()
            maxExtent[dendType] = max(values)

        return maxExtent

    #####################################################################
    # Function that determines the segments of the sections that belong #
    # to a particular layer and the coordinates of the segments         #
    #####################################################################
    def orgSeg(self,dendTypeList,layerDict,bounds,segLocDict):
        for dendType in self.dendTypeList:
            dendList = self.dendTypeList[dendType]
            # Calculate which segments on which sections belong to which layer
            self.organizeSec(dendList,self.layerDict[dendType],self.bounds[dendType])

            # Dictionary of dictionary of tuples that will be used to give
            # locations to synapses
            for i in range(len(dendList)):
                sec = dendList[i]

                # Figure out which swc points the segment lies between
                for layer in self.segLocDict[dendType]:
                    self.find_segloc(sec,self.layerDict[dendType][layer],self.segLocDict[dendType][layer])

    #####################################################################
    # Function that determines the segments of the sections that belong #
    # to a particular layer                                             #
    #####################################################################
    def organizeSec(self,dend,layerDict,bounds):
        for i in range(len(dend)):
            sec = dend[i]
            for layer in layerDict:
                layerDict[layer][sec] = []

            # Interpolating the points in between consecutive pairs of swc
            # coordinates in order to increase the spatial resolution of
            # the algorithm.
            x,y,z = self.interpSWC(sec)

            # Calculate extent and normalized distance of points along section
            extent, norm_dist = self.secExtent(x,y,z,sec)

            # Figure out which points lie in which layer of the DG and the corresponding
            # normalized distances
            for ii in range(len(extent)):
                for layer in bounds:
                    if bounds[layer][0] <= extent[ii] < bounds[layer][1]:
                        layerDict[layer][sec].append(norm_dist[ii])

            # Calculate which segment the points are closest to
            for layer in layerDict:
                self.point2seg(sec,layerDict[layer])

    #######################################################################
    # Function that linearly interpolates points between swcs to increase #
    # the resolution of the segment organizing function (organizeSec)     #
    #######################################################################
    def interpSWC(self,sec):
        sec.push()
        x = []
        y = []
        z = []
        for ii in range(int(h.n3d())-1):
            x_current = h.x3d(ii)
            y_current = h.y3d(ii)
            z_current = h.z3d(ii)
            x_next = h.x3d(ii+1)
            y_next = h.y3d(ii+1)
            z_next = h.z3d(ii+1)
            x += list(np.linspace(x_current,x_next,sec.nseg+2))[:-1]
            y += list(np.linspace(y_current,y_next,sec.nseg+2))[:-1]
            z += list(np.linspace(z_current,z_next,sec.nseg+2))[:-1]
        x.append(h.x3d(ii+1))
        y.append(h.y3d(ii+1))
        z.append(h.z3d(ii+1))
        h.pop_section()
        return x,y,z

    #######################################################################
    # Function that determines the distance of swcs of a section from the #
    # soma along the major axis                                           #
    #######################################################################
    def secExtent(self,x,y,z,sec):
        extent = []
        norm_dist = []
        dist = 0
        for ii in range(len(x)):
            if ii > 0:
                x2 = (x[ii]-x[ii-1])**2
                y2 = (y[ii]-y[ii-1])**2
                z2 = (z[ii]-z[ii-1])**2
                dist += ( x2 + y2 + z2 )**0.5
            norm_dist.append(dist/sec.L)
            loc = np.array([x[ii],y[ii],z[ii]])
            extent.append(abs(np.dot(self.new_axis,loc-self.center)))
        return extent,norm_dist

    ######################################################################
    # Function that determines which segment the interpolated points are #
    # closest to                                                         #
    ######################################################################
    def point2seg(self,sec,reg):
        seg_pos = np.linspace(0.5/sec.nseg,1-0.5/sec.nseg,sec.nseg)
        for ii in range(len(reg[sec])):
            idx = np.abs(seg_pos-reg[sec][ii]).argmin()
            reg[sec][ii] = seg_pos[idx]

        reg[sec] = list(set(reg[sec]))
        reg[sec].sort()


    ##########################################################
    # Function that extracts the coordinates of the segments #
    ##########################################################
    def find_segloc(self,sec,reg,loclist):
        sec.push()
        if len(reg[sec]) > 0:
            for seg in reg[sec]:
                for ii in range(int(h.n3d())-1):
                    if (h.arc3d(ii)/sec.L) <= seg < (h.arc3d(ii+1)/sec.L):
                        swc1 = np.array([h.x3d(ii),h.y3d(ii),h.z3d(ii)])
                        swc2 = np.array([h.x3d(ii+1),h.y3d(ii+1),h.z3d(ii+1)])
                        f = (seg-h.arc3d(ii)/sec.L)/((h.arc3d(ii+1)-h.arc3d(ii))/sec.L)
                        break
                point = f*(swc2-swc1)+swc1
                loclist[0].append(point[0])
                loclist[1].append(point[1])
                loclist[2].append(point[2])
        h.pop_section()

    '''#############################
    # 3. Code for connecting cells #
    #############################'''
    def connect_pre(self, syn, wt, dly):
        # This function attaches my axon (or other generator of APs), to the synapse
        # that's passed to me from another cell.

        # The section of the cell that is going to be the AP generator has to be the
        # currectly accessed section, so do that...
        self.celltype.getSoma(self).push()
        nc = h.NetCon(self.celltype.getSoma(self)(1)._ref_v, syn, 0, dly, wt)
        self.nc_list.append(nc)
        h.pop_section()
        return nc

    def do_connect(self,inputs,dummy,params,pre,sharedData,nmda_flag):
        synweight = params[0]
        stype = params[1]
        layer = params[2]
        sdelay = params[3]
        vel = params[4]
        AMPA_NMDA_ratio = params[5]
        tau = params[6]
        # The GC-CA3 projection has a slightly different form. It contains both the presynaptic
        # cell IDs and the mossy fiber distances traveled to reach the postsynaptic neuron
        if len(inputs) > 0:
            if isinstance(inputs[0],np.ndarray):
                IDs = inputs[0]
                distances = inputs[1]
            else:
                IDs = inputs

            for iSyn in range(len(IDs)):
                inCellNum = IDs[iSyn]

                # Different initialization for STDP
                if 'STDP' in self.synvars['type']:
                    weight = synweight
                else:
                    weight = [0,0]
                    weight[0] = self.synvars['wmin']
                    weight[1] = synweight

                # Need to choose a segment in the right layer
                choice = self.ranGen.randint(0,len(self.synGroups[stype][layer])-1)
                self.synGroups[stype][layer][choice].tau2 = tau
                #self.synGroups[stype][layer][choice].tau1 = tau_dict[post][pre][0]

                nc = pc.gid_connect(inCellNum, self.synGroups[stype][layer][choice])
                if isinstance(inputs[0],np.ndarray):
                    dist = distances[iSyn]
                else:
                    if isinstance(sharedData[inCellNum],tuple):
                        dist = self.euclidean(self.location,sharedData[inCellNum])
                    else:
                        dist = self.euclidean(self.location,sharedData[inCellNum][0])

                delay = sdelay+dist/vel
                nc.weight[0] = synweight
                nc.delay = delay
                self.nc_list.append(nc)
                if nmda_flag == 1:
                    nc = pc.gid_connect(inCellNum, self.synGroups['NMDA'][layer][choice])
                    nc.weight[0] = AMPA_NMDA_ratio*weight
                    nc.delay = delay
                    self.nc_list.append(nc)

    '''###########################################
    # 4. Code for populating cells with synapses #
    ###########################################'''
    #######################################################################
    # This function will place a synapse on each segment of the model and #
    # append the synapse to the corresponding list based on which layer   #
    # the segment lies                                                    #
    #######################################################################
    def addSegmentSynapse(self,layerDict,list,synList,synType):

        for layer in sorted(layerDict.keys()):
            for sec_choice in list:
                dict = layerDict[layer]
                for seg_choice in dict[sec_choice]:
                    print layer, sec_choice, seg_choice
                    nmda_flag = 0
                    if self.synvars['type'] == "E2-NMDA2":
                        syn = h.Exp2Syn(sec_choice(seg_choice))
                        nmda = h.Exp2NMDA_Wang(sec_choice(seg_choice))
                        nmda_flag = 1
                    if self.synvars['type'] == "E2":
                        syn = h.Exp2Syn(sec_choice(seg_choice))
                    if self.synvars['type'] == "E2_Prob":
                        syn = h.E2_Prob(sec_choice(seg_choice))
                        syn.P = self.synvars['P']
                    if self.synvars['type'] == "E2_STP_Prob":
                        syn = h.E2_STP_Prob(sec_choice(seg_choice))
                    if self.synvars['type'] == "STDPE2":
                        syn = h.STDPE2(sec_choice(seg_choice))
                    if self.synvars['type'] == "STDPE2_Clo":
                        syn = h.STDPE2_Clo(sec_choice(seg_choice))
                    if self.synvars['type'] == "STDPE2_STP"	:
                        syn = h.STDPE2_STP(sec_choice(seg_choice))
                    if self.synvars['type'] == "STDPE2_Prob":
                        syn = h.STDPE2_Prob(sec_choice(seg_choice))
                        syn.P = self.synvars['P']
                    #initializes different variables depending on synapse
                    if (self.synvars['type'] == "STDPE2_STP")|(self.synvars['type'] == "E2_STP_Prob"):
                        syn.F1 = self.synvars['F1']
                    if  (self.synvars['type'] == "STDPE2_Clo" )|( self.synvars['type'] == "STDPE2_STP")|( self.synvars['type'] == "STDPE2")| (self.synvars['type'] == "STDPE2_Prob"):
                        syn.wmax = self.synvars['wmax']
                        syn.wmin = self.synvars['wmin']
                        syn.thresh = self.synvars['thresh']
                    if  (self.synvars['type'] == "E2_Prob" )|( self.synvars['type'] == "E2_STP_Prob")|(self.synvars['type'] == "STDPE2_STP") | (self.synvars['type'] == "STDPE2_Prob"):
                        h.use_mcell_ran4(1)
                        syn.seed = self.ranGen.randint(1,4.295e9)
                    syn.tau1 = 0.5
                    syn.tau2 = 0.6
                    syn.e = 0
                    if nmda_flag == 1:
                        if synType == 'NMDA':
                            synList[layer].append(nmda)
                        else:
                            synList[layer].append(syn)
                    else:
                        synList[layer].append(syn)

    '''###########################################
    # 5. Code to save out data from each segment #
    ###########################################'''
    def setdata(self):
        start_ends = []
        volts = []
        for dendType in self.dendTypeList:
            for sec in self.dendTypeList[dendType]:
                sec.push()
                seg_pos = np.linspace(0.5/sec.nseg,1-0.5/sec.nseg,sec.nseg)
                boundary_pos = np.linspace(0,1,sec.nseg+1)
                for jj in range(len(seg_pos)):
                    pos = seg_pos[jj]
                    tmp_v = h.Vector()
                    tmp_v.record(sec(pos)._ref_v,0.25)
                    start_pos = boundary_pos[jj]
                    end_pos = boundary_pos[jj+1]
                    for ii in range(int(h.n3d())-1):
                        if (h.arc3d(ii)/sec.L) <= start_pos <= (h.arc3d(ii+1)/sec.L):
                            swc1 = np.array([h.x3d(ii),h.y3d(ii),h.z3d(ii)])
                            swc2 = np.array([h.x3d(ii+1),h.y3d(ii+1),h.z3d(ii+1)])
                            f = (start_pos-h.arc3d(ii)/sec.L)/((h.arc3d(ii+1)-h.arc3d(ii))/sec.L)
                            break
                    start_point = tuple(f*(swc2-swc1)+swc1)
                    for ii in range(int(h.n3d())-1):
                        if (h.arc3d(ii)/sec.L) <= end_pos <= (h.arc3d(ii+1)/sec.L):
                            swc1 = np.array([h.x3d(ii),h.y3d(ii),h.z3d(ii)])
                            swc2 = np.array([h.x3d(ii+1),h.y3d(ii+1),h.z3d(ii+1)])
                            f = (end_pos-h.arc3d(ii)/sec.L)/((h.arc3d(ii+1)-h.arc3d(ii))/sec.L)
                            break
                    end_point = tuple(f*(swc2-swc1)+swc1)
                    volts.append(tmp_v)
                    start_ends.append([start_point,end_point])
                h.pop_section()
        return [volts, start_ends]

    '''############################
    # 6. Code to test EPSPs/IPSPs #
    ############################'''
    def testSingleSyn(self,syn_type,layer,idx,t1,t2,weight):
        Input = h.VecStim()
        evec = h.Vector([200])
        Input.play(evec)
        self.synGroups[syn_type][layer][idx].tau1 = t1
        self.synGroups[syn_type][layer][idx].tau2 = t2
        if syn_type == 'GABA':
            self.synGroups[syn_type][layer][idx].e = -75
        else:
            self.synGroups[syn_type][layer][idx].e = -0
        self.syn_con = h.NetCon(Input,self.synGroups[syn_type][layer][idx])
        self.syn_con.delay = 0
        self.syn_con.weight[0] = weight
        self.syn_con.threshold = 10.0

    '''##################
    # 7. Static methods #
    ##################'''
    @staticmethod
    def euclidean(loc1,loc2):
        return ((loc1[0]-loc2[0])**2 + (loc1[1]-loc2[1])**2)**0.5
