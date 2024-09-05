import pandas as pd
import numpy as np
import hypernetx as hnx
import hypernetx.algorithms.hypergraph_modularity as hmod
from collections import defaultdict
import copy
import math
import csv
from collections import Counter
from scipy.stats import binom 

#from h_louvain_utils import last_step


class hLouvain:

    """
    hLouvain algorithm main class

    Parameters
    ----------

    HG : hypernetx.hypergraph.

    hmod_tau : str | float, optional, default= "infinity"
        "infinity" or float greater or equal to 0
        Hyperparameter for hypergraph modularity 
    
    resolution : float, optional, default = 1 (corresponding to the standard settings for modularity calculation) 
        float greater or equal than 0
        Hyperparameter for hypergraph modularity, controlling the weight of degree tax. 
        The classical definition of modularity is retrieved when the resolution parameter is set to 1.

    Note
    ----
    For 'tau', any float >= 0 can be used. 
    Basically it corresponds to hyperedge weights (c/d)^tau when c > d/2 and 0 otherwise.
    Default is "infinity" (strict modularity).
    Other natural choices are 0 -  majority modularity, 1 - linear modularity, 2 - quadratic modularity.

    """
    
    def __init__(
        self,
        HG,
        hmod_tau = "infinity",
        resolution = 1,
        random_seed = 123
    ) -> None:
        self.HG = HG
        self.hmod_tau = hmod_tau
        self.resolution = resolution
        self.random_seed = random_seed
        self.G = hmod.two_section(HG)
        
        # setting the order of nodes based on their apperance in edges
        self.h_nodes = []
        for e in self.HG.edges:
            E = self.HG.edges[e]
            for node in E:
                if node not in self.h_nodes:
                    self.h_nodes.append(node)
      
        
        self.edge_sizes = [self.HG.size(i) for i in self.HG.edges]
        self.startHGdict, self.startA = self._build_HG_dict_from_HG()
        self.HGdict = {},
        self.dts = defaultdict(float) #dictionary containing degree taxes for given volumes
        self.wdc = self._calculate_wdc()

        self.edge_weights = [self.HG.edges[j].weight for j in self.HG.edges]
        self.uniq = len(Counter(self.edge_weights)) == 1 
        self.changes = 0
        self.communities = len(self.startA)
        self.phase = 1
        self.iteration = 1
        self.change_mode = "communities"
        self.d_weights = self.d_weights_calculation()
        self.total_weight = sum(self.d_weights.values())
        self.no_edges = len(self.HG.edges)              
        self.forsing_alpha_1 = True
        self.modified  = False
        self.improved = False

    def _calculate_wdc(self):
        Ctr = Counter(self.edge_sizes)
        wdc = {}
        for d in Ctr.keys():
            wdc[d]={}
            for c in range(d+1):
                if self.hmod_tau == "infinity":
                    wdc[d][c] = 1 if c == d else 0
                else:
                    wdc[d][c] = (c/d)**self.hmod_tau if c > d / 2 else 0
        return wdc
    
    def h_modularity(self, A, tau = "infinity", resolution = 1):

        """
        Computes modularity of hypergraph HG with respect to partition A.

        Parameters
        ----------
        A : list of sets
            Partition of the vertices in HG
        tau : str | float, optional, default= "infinity"
            "infinity" or float greater or equal to 0
            Hyperparameter for hypergraph modularity 
        resolution : float, optional, default = 1 (corresponding to the standard settings for modularity calculation) 
            float greater or equal than 0
            Hyperparameter for hypergraph modularity, controlling the weight of degree tax. 
            The classical definition of modularity is retrieved when the resolution parameter is set to 1.

        Note
        ----
        For 'tau', any float >= 0 can be used. 
        Basically it corresponds to hyperedge weights (c/d)^tau when c > d/2 and 0 otherwise.
        Default is "infinity" (strict modularity).
        Other natural choices are 0 -  majority modularity, 1 - linear modularity, 2 - quadratic modularity.

        Returns
        -------
        : float
        The modularity function for partition A on HG
        """
        
        ## calculating wdc
        Ctr = Counter(self.edge_sizes)
       
        wdc = {}
        for d in Ctr.keys():
            wdc[d]={}
            for c in np.arange(int(np.floor(d / 2 + 1)), d + 1):
                if tau == "infinity":
                    wdc[d][c] = 1 if c == d else 0
                else:
                    wdc[d][c] = (c/d)**tau if c > d / 2 else 0
        

        ## all same edge weights?
        #uniq = (len(Counter(self.HG.edges.properties["weight"]))==1)
        uniq = self.uniq

        ## Edge Contribution
        HG_id = self.HG.incidence_dict
        d = hmod.part2dict(A)
        L = [[d[i] for i in HG_id[x]] for x in HG_id]

        ## all same weight
        if uniq:
            _ctr = Counter([(Counter(l).most_common(1)[0][1], len(l)) for l in L])
            EC = sum([wdc[k[1]][k[0]] * _ctr[k] for k in _ctr.keys() if k[0] > k[1] / 2])
        else:
            _keys = [(Counter(l).most_common(1)[0][1], len(l)) for l in L]
            _vals = list(self.HG.edges.dataframe["weight"])
            _df = pd.DataFrame(zip(_keys, _vals), columns=["key", "val"])
            _df = _df.groupby(by="key").sum()
            EC = sum(
                [wdc[k[1]][k[0]]* v.iloc[0] for (k, v) in _df.iterrows() if k[0] > k[1] / 2]
            )

        ## Degree Tax
        if uniq:
            VolA = [sum([self.HG.degree(i) for i in k]) for k in A]
            

        else:
            ## this is the bottleneck
            VolA = np.repeat(0, 1 + np.max(list(d.values())))
            m = np.max(self.edge_sizes)
            Ctr = np.repeat(0, 1 + m)
            S = 0
            for e in self.HG.edges:
                w = self.HG.edges[e].weight
                Ctr[self.HG.size(e)] += w
                S += w
                for v in self.HG.edges[e]:
                    VolA[d[v]] += w

        VolV = np.sum(VolA)
        VolA = [x / VolV for x in VolA]
        DT = 0

        if uniq:
            for d in Ctr.keys():
                Cnt = Ctr[d]
                for c in np.arange(int(np.floor(d / 2 + 1)), d + 1):
                    for Vol in VolA:
                        DT += Cnt * wdc[d][c] * binom.pmf(c, d, Vol)
            return EC/ self.no_edges - (resolution * DT / self.no_edges)
        else:
            for d in range(len(Ctr)):
                Cnt = Ctr[d]
                if Cnt > 0:
                    for c in np.arange(int(np.floor(d / 2 + 1)), d + 1):
                        for Vol in VolA:
                            DT += Cnt * wdc[d][c] * binom.pmf(c, d, Vol)
            return EC/S - (resolution * DT/S)




    #precalculation of d_weights (for speed up of degree tax change calculation)

    def d_weights_calculation(self):

        #uniq = (len(Counter(self.HG.edges.properties["weight"]))==1)
        uniq = self.uniq

        if uniq:        
            Ctr = Counter(self.edge_sizes)
            

        else:
            
            m = np.max(self.edge_sizes)
            Ctr2 = np.repeat(0,1+m)
            for e in self.HG.edges:
                w = self.HG.edges[e].weight
                Ctr2[self.HG.size(e)] += w  

            Ctr = Counter([len(self.HG.edges[e]) for e in self.HG.edges])
            for k in range(len(Ctr2)):
                if Ctr2[k] > 0:
                    Ctr[k] = Ctr2[k]
    
        
        return Ctr

    def _setHGDictAndA(self):
        self.HGdict = copy.deepcopy(self.startHGdict)
        A = copy.deepcopy(self.startA)
        return A
    
    
    def _hg_neigh_dict(self):
        """
        Optimizing the access to nodes neighbours by additional dictionary based on hnx neighbours() function

        Parameters
        HG : an HNX hypergraph

        Returns
        -------
        : dict
        a dictionary with {node: list of neighboring nodes}
        """

        result = {}
        for v in self.h_nodes:
            result[v] = self.HG.neighbors(v)
        return result


    def combined_modularity(self,A, alpha=0.5, hmod_tau="infinity", resolution=1):
        
        # Exclude empty clusters if exist
        A = [a for a in A if len(a) > 0]
        
        # Calculation of hypergraph modularity 
        hyper = self.h_modularity(A, hmod_tau, resolution) 
        
        # Calculation of modularity of 2-section graph (based on iGraph)
        d = hmod.part2dict(A)
        partition  = [d[i] for i in self.h_nodes]
        twosect = self.G.modularity(partition,weights='weight',resolution=resolution)  # weighting is enabled
        
        return (1-alpha)*twosect+alpha*hyper




    def _ec_loss_h(self, c, s):
        
        ec_h = 0  # ec for hypergraph
        
        #Edge-by-edge accumulation of edge contribution loss (for hypergraph and 2section)
        for e in self.HGdict["v"][s]["e"].keys():
            d = self.HGdict["e"][e]["size"]
            change = self.HGdict["v"][s]["e"][e] 
            if d > 1:  #we ignore hyperedges with one vertex
                w = self.HGdict["e"][e]["weight"]
                old_c = self.HGdict["c"][c]["e"][e] # counting vertices of edge e in community c
                new_c = old_c - change #updating after taking the supernode     
                ec_h += w * (self.wdc[d][old_c]- self.wdc[d][new_c])
        return ec_h / self.total_weight
    

    def _ec_loss_2s(self, c, s):
        
        ec_2s = 0  #ec for 2section graph
        
        #Edge-by-edge accumulation of edge contribution loss (for hypergraph and 2section)
        for e in self.HGdict["v"][s]["e"].keys():
            d = self.HGdict["e"][e]["size"]
            if d > 1:  #we ignore hyperedges with one vertex
                change = self.HGdict["v"][s]["e"][e] 
                new_c = self.HGdict["c"][c]["e"][e]  - change #updating after taking the supernode     
                if new_c > 0:
                    w = self.HGdict["e"][e]["weight"]
                    ec_2s += w * change * new_c/(d-1)
        return 2*ec_2s / self.total_volume 


    def _ec_gain_h(self, c, s):
        
        '''
        Calculating edge contribution gain when joinig supernode s to destination community c
        '''

        ec_h = 0
        #Edge-by-edge accumulation of edge contribution change 
        for e in self.HGdict["v"][s]["e"].keys(): #taking all edges of supernode "s"
            d = self.HGdict["e"][e]["size"]
            input = self.HGdict["v"][s]["e"][e]
            if d > 1:  #we ignore hyperedges with one vertex
                
                if (input > d/2) or (e in self.HGdict["c"][c]["e"].keys()):
                    w = self.HGdict["e"][e]["weight"] #weight of the edge
                    old_c = self.HGdict["c"][c]["e"][e]  # current number of vertices of edge e in community c
                    new_c = old_c + input #updated number of vertices after moving a supernode s to community c 
                    ec_h += w * (self.wdc[d][new_c]- self.wdc[d][old_c]) # edge contribution gain (hypergraph)

        return ec_h / self.total_weight
    
    def _ec_gain_2s(self, c, s):
        

        ec_2s = 0

        #Edge-by-edge accumulation of edge contribution change 
        for e in self.HGdict["v"][s]["e"].keys(): #taking all edges of supernode "s"
            d = self.HGdict["e"][e]["size"]
            if d > 1:  #we ignore hyperedges with one vertex          
                if e in self.HGdict["c"][c]["e"].keys():   
                    input = self.HGdict["v"][s]["e"][e]
                    w = self.HGdict["e"][e]["weight"] #weight of the edge
                    old_c = self.HGdict["c"][c]["e"][e]  # current number of vertices of edge e in community c
                    ec_2s += w * input * old_c/(d-1) #cummulating edge contribution gain (2section)
        return 2*ec_2s / self.total_volume 


    def _degree_tax_hyper(self, volume):
        
        volume = volume/self.total_volume
        
        DT = 0
        for d in self.d_weights.keys():
            x = 0
            for c in np.arange(int(np.floor(d / 2)) + 1, d + 1):
                x += self.wdc[d][c] * binom.pmf(c, d, volume)
            DT += x * self.d_weights[d]
        return DT / self.total_weight

    def _combined_delta(self,delta_h, delta_2s,alpha=0.5):
        return (1-alpha)*delta_2s+alpha*delta_h

    

    def _build_HG_dict_from_HG(self):
        
        '''
        Procedure for building dictionary-based data structures for hypergraph with supernodes
        Procedure produces the initial structure based on the input hypergraph with all supernodes with just only node
        '''

        HGdict = {}
        HGdict["v"]={}  #subdictionary storing information on vertices
        HGdict["c"]={}  #subdictionary storing information on current communities
        A = []          #current partition as a list of sets
        name2index = {} 
        
        # preparing structures for vertices and initial communities
        
        for i in range(len(self.h_nodes)):
            HGdict["v"][i]={}
            HGdict["v"][i]["name"] = [self.h_nodes[i]]
            HGdict["v"][i]["e"] = defaultdict(int) #dictionary  edge: edge_dependent_vertex_weight   
            HGdict["v"][i]["strength"]=0   #strength = weigthed degree
            HGdict["c"][i]={}
            HGdict["c"][i]["e"] = defaultdict(int) #dictionary  edge: edge_dependent_community_weight
            HGdict["c"][i]["strength"]=0
            A.append({self.h_nodes[i]})
            name2index[self.h_nodes[i]] = i

        # preparing structures for edges 
            
        HGdict["e"] = {}
        edge_dict = copy.deepcopy(self.HG.incidence_dict)
        
        HGdict["total_volume"] = 0    # total volume of all vertices             
        HGdict["total_weight"] = sum(self.HG.edges[j].weight for j in self.HG.edges) #sum of edges weights
        
        for j in range(len(self.HG.edges)):
            HGdict["e"][j] = {}
            HGdict["e"][j]["v"] = defaultdict(int)  #dictionary  vertex: edge_dependent_vertex_weight
            HGdict["e"][j]["weight"] = copy.deepcopy(self.HG.edges[j].weight)
            HGdict["e"][j]["size"] = copy.deepcopy(self.HG.size(j))
            HGdict["e"][j]["volume"] = HGdict["e"][j]["weight"]*HGdict["e"][j]["size"]    
            for t in edge_dict[j]:
                HGdict["e"][j]["v"][name2index[t]] += 1 
                HGdict["v"][name2index[t]]["e"][j] += 1
                HGdict["c"][name2index[t]]["e"][j] += 1
                HGdict["v"][name2index[t]]["strength"] += HGdict["e"][j]["weight"]
                HGdict["c"][name2index[t]]["strength"] += HGdict["e"][j]["weight"]  
                HGdict["total_volume"] += HGdict["e"][j]["weight"]

                
        return HGdict, A


    def _node_collapsing(self,HGd,A):

        '''
        Procedure for node collapsing
        Procedure produces supervertices based on communities from the previous phase
        '''
        
        # reindexing the partition from previous phase (skipping empty indexes)
        newA = [a for a in A if len(a) > 0]
        
        # building new vertices and communities (by collapsing nodes to supernodes)
        
        HGdict = {}
        HGdict["v"]={}
        HGdict["c"]={}
        i = 0
        for j in range(len(A)):
            if len(A[j])>0:
                HGdict["v"][i]={}
                HGdict["c"][i]={}
                HGdict["v"][i]["e"]={}
                HGdict["c"][i]["e"]=defaultdict(int)
                for e in HGd["c"][j]["e"].keys():
                    if HGd["c"][j]["e"][e] > 0:
                        HGdict["v"][i]["e"][e] = HGd["c"][j]["e"][e]
                        HGdict["c"][i]["e"][e] = HGd["c"][j]["e"][e]
                HGdict["v"][i]["strength"] = HGd["c"][j]["strength"]
                HGdict["c"][i]["strength"] = HGd["c"][j]["strength"]
                i+=1        

        #updating edges based on new supernodes indexes and weights
        
        HGdict["e"] = HGd["e"]
        HGdict["total_volume"] = HGd["total_volume"]
        HGdict["total_weight"] = HGd["total_weight"]
        for j in range(len(HGdict["e"])):
            HGdict["e"][j]["v"] = {}
            HGdict["e"][j]["v"] = defaultdict(int)  
        for v in HGdict["v"].keys():
            for e in HGdict["v"][v]["e"].keys():
                HGdict["e"][e]["v"][v] = HGdict["v"][v]["e"][e] 
                
        return HGdict, newA


    def _next_maximization_iteration(self, L, DL, A1, D, 
                alpha = 0.5):
        
        '''
        Function for hLouvain modularity maximization iteration (checking all vertices in the random order)
        '''
        
        current_alpha = alpha
        A = copy.deepcopy(A1)

        local_length = len(A)
        step = math.ceil(local_length/10)
        self.improved = False
        
        for sn in list(np.random.RandomState(self.random_seed).permutation(L)): #permute supernodes (L contains initial partition, so supernodes)
            #check the current cluster number (using the first node in supernode)
            c = D[list(sn)[0]] 
            #check the supernode index
            si = DL[list(sn)[0]] 
            
            # matrix with modularity deltas
            
            M = []       
            superneighbors = list(set(np.concatenate([[D[v] for v in self.neighbors_dict[w]] for w in sn])))
            if len(superneighbors) > 1 or c not in superneighbors: #checking only clusters with some neighbors of current supernodes
                               
                # calculating loss in degree tax for hypergraph when taking vertex si from community c
                vols = self.HGdict["v"][si]["strength"]
                volc = self.HGdict["c"][c]["strength"]

               

                if not self.dts[volc]:
                    self.dts[volc] = self._degree_tax_hyper(volc)
                if volc > vols and not self.dts[volc - vols]:
                    self.dts[volc-vols] = self._degree_tax_hyper(volc - vols)
                
                dt_h_loss = self.dts[volc] - self.dts[volc-vols] 
                
                # calculating loss in edge contribution for hypergraph and 2section when taking vertex si from community c
                
                ec_h_loss = 0
                if current_alpha > 0: 
                    ec_h_loss = self._ec_loss_h( c, si)
                
                ec_2s_loss = 0
                if current_alpha < 1:
                    ec_2s_loss = self._ec_loss_2s( c, si)
                                
                for i in superneighbors:   
                    if c == i:
                        M.append(0) 
                    else:
                        # gain in degree tax for hypergraph when joining vertex si to community i
                
                        voli = self.HGdict["c"][i]["strength"]
                   
                   
                        if not self.dts[voli]:
                            self.dts[voli] = self._degree_tax_hyper(voli)
                        if not self.dts[voli + vols]:
                            self.dts[voli + vols] = self._degree_tax_hyper(voli + vols)
                        dt_h_gain = self.dts[voli + vols] - self.dts[voli]
                        
                        # change in degree tax for 2section graph
                        delta_dt_2s =2*vols*(vols+voli-volc)/(self.total_volume**2)
                        
                        # gains in edge contribution when joining si to community i
                        ec_h_gain = 0
                        if current_alpha > 0:
                            ec_h_gain = self._ec_gain_h(i, si)
                        ec_2s_gain = 0
                        if current_alpha < 1:
                            ec_2s_gain = self._ec_gain_2s(i, si)
                        
                        #calulating deltas
                        delta_h = ec_h_gain - ec_h_loss - (self.resolution * (dt_h_gain - dt_h_loss))  # self.resolution
                        delta_2s = ec_2s_gain - ec_2s_loss - (self.resolution * delta_dt_2s)  #self.resolution
                        
                        M.append(self._combined_delta(delta_h, delta_2s, current_alpha))
                # make a move maximizing the gain 
                if max(M) > 0:
                    i = superneighbors[np.argmax(M)]
                    for v in list(sn):
                        A[c] = A[c] - {v}
                        A[i] = A[i].union({v})
                        D[v] = i

                    self.changes+=1
                    self.new_changes+=1
                    if len(A[c]) == 0:
                        self.communities-=1

                    self.improved = True
                    self.modified = True


                    if self.change_mode == "modifications":
                        if self.new_changes >= self.after_changes:
                            self.level+=1
                            self.new_changes = 0
                            if self.level > len(self.alphas) - 1:
                                self.alphas.append(self.alphas[-1])
                            current_alpha = self.alphas[self.level]

                    if self.change_mode == "communities":
                        if self.communities <= self.community_threshold:
                            self.level+=1
                            self.community_threshold = self.community_threshold/self.community_factor
                            if self.level > len(self.alphas) - 1:
                                self.alphas.append(self.alphas[-1])
                            current_alpha = self.alphas[self.level]

                                        
                    for e in self.HGdict["v"][si]["e"].keys():
                        self.HGdict["c"][c]["e"][e] -= self.HGdict["v"][si]["e"][e]
                        if self.HGdict["c"][c]["e"][e] == 0:
                            del self.HGdict["c"][c]["e"][e]
                        self.HGdict["c"][i]["e"][e] += self.HGdict["v"][si]["e"][e]
                        self.HGdict["c"][c]["strength"] -= self.HGdict["v"][si]["e"][e]*self.HGdict["e"][e]["weight"]
                        self.HGdict["c"][i]["strength"] += self.HGdict["v"][si]["e"][e]*self.HGdict["e"][e]["weight"]   


        return A, D


    def _next_maximization_phase(self, L):

        '''
        Function for single hLouvain pass (maximization phase with the node collapsing at the end)
        '''
        DL = hmod.part2dict(L)
        A1 = L[:]  
        D = hmod.part2dict(A1) #current partition as a dictionary (D will change during the phase)   
        self.modified = False

        if self.back == True:
            L = self.LPrevious
            self.HGdict = self.HGdictPrevious
            DL = self.DLPrevious
            D = self.DPrevious
            A1 = self.APrevious
            self.back = False
        
        
        self.total_volume = self.HGdict["total_volume"]

        if self.change_mode == "iter":
            if self.level > len(self.alphas) - 1:
                self.alphas.append(self.alphas[-1])
        
        while True:

            A2, D =  self._next_maximization_iteration(L, DL, A1, D, self.alphas[self.level])

            if self.improved == False:

                if self.iteration != 1 and self.forsing_alpha_1 == True:
                    self.HGdictPrevious = copy.deepcopy(self.HGdict)
                    self.DPrevious = copy.deepcopy(D)
                    self.APrevious = copy.deepcopy(A2)
                    self.DLPrevious = copy.deepcopy(DL)
                    self.LPrevious = copy.deepcopy(L)


                self.HGdict, newA = self._node_collapsing(self.HGdict,A2)
                

                break

            self.iteration += 1
            if self.change_mode == "iter":
                self.level+=1
                if self.level > len(self.alphas) - 1:
                    self.alphas.append(self.alphas[-1])
                
            A1 = A2
        
        return newA




    def h_louvain_community(self,alphas = [1], change_frequency = 0.5, change_mode = "communities", 
                            after_changes = 300, force_alpha_1 = True, random_seed = 1):


        '''
        The default hLouvain community detection algorithm version

        Parameters
        ----------

        HG : hypernetx.hypergraph.

        alphas : array of floats, optional, default= [1] (corresponding to the option of maximizing hypergraph modularity from the beginning to the end of the process)
            The consecutive values of alpha parameter used during the algorithm.
            The array is not-empty.
            Hyperparameter for hypergraph modularity 
        
        change_mode : string, optional, default = "communities" 
            other options: "iter", "modifications"
            Parameter controling  when the change of alpha is made.
            The default value "communities" is based on observing the current number of communities
            Other options: 
            "modification" - based on the number of modifications (node moves) made.
            "iter" - the option forsing alpha change after each iteration of modularity maximization.

        change_frequency : float, optional, default = 0.5
            The parameter used change_mode = "communities" controlling the frequency of changing alpha.
            E.g. value 0.5 forces the change when number of communities is halved when comparing with the value from the previous change.
        
        after_changes: int, optional, default = 300
            The parameter used change_mode = "modifications" controlling the frequency of changing alpha.
            E.g. value 300 forces the change every 300 modifications (node moves).

        force_alpha_1 : boolean, optional, default = True
            The parameter used when alpha is not equal to 1 at the end of algorithm, 
            e.i., the final maximization function is still not the pure hypergraph modularity (alpha < 1)
            When set to False the algorithm finishes despide the fact that alpha is smaller than 1.
            When set to True the algorithm back to its previous pass and do the additional optimization phase with alpha set to 1.

        random_seed: int

        Returns
        -------
        partition: list of sets
        h-modularity of the partition: float
        the sequence of alphas used: list of floats



        Note
        ----
            alpha array is not empty and can be of any length. 
            When algorithm needs to change alpha it always takes the next value from this array.
            If there is no next value then the last value is repeated. 
            E.g. default array [1] corresponds to infinite sequence [1,1,1,...]
        '''

        A1 = self._setHGDictAndA()
        if random_seed > 1:
            self.random_seed = random_seed
        self.changes = 0
        self.new_changes = 0
        self.communities = len(A1)
        self.phase = 1
        self.change_mode = change_mode
        self.after_changes = after_changes
        self.alphas = alphas
        self.level = 0
        self.forsing_alpha_1 = force_alpha_1
        self.iteration = 1
        self.community_factor = 1/change_frequency
        self.community_threshold = self.communities/self.community_factor
        self.back = False
        
        self.neighbors_dict = copy.deepcopy(self._hg_neigh_dict())
        
        while True:
            if self.back == True:
                A1 = self.LPrevious
            A2  = self._next_maximization_phase(A1)
           
            
            A1 = A2
            if self.modified == False:
                if self.alphas[self.level] == 1 or self.forsing_alpha_1 == False:  
                    qnew = self.combined_modularity(A2, self.alphas[self.level], self.hmod_tau, self.resolution) 
                    return A2, qnew, self.alphas[0:self.level+1]
                else:
                    # added in order to force checking alpha = 1 at the last step of the process
                    if self.level > len(self.alphas) - 2:
                        self.alphas.append(1)
                    else:
                        for i  in np.arange(self.level+1, len(self.alphas)):
                            self.alphas[i] = 1 
                    self.back = True
                    self.level+=1
            self.phase+=1
            self.iteration = 1
            if self.change_mode == "phase":
                self.level+=1

                if self.level > len(self.alphas) - 1:
                    self.alphas.append(self.alphas[-1])


    def h_louvain_community_plus_last_step(self,alphas = [1],  change_frequency = 0.5, 
                                           change_mode = "communities", after_changes = 300,  
                                           force_alpha_1 = True, random_seed = 1):
        
        '''
        The enhanced hLouvain community detection algorithm by the additional step based on hypernetx last_step algorithm
        
        Parameters
        ----------

        HG : hypernetx.hypergraph.

        alphas : array of floats, optional, default= [1] (corresponding to the option of maximizing hypergraph modularity from the beginning to the end of the process)
            The consecutive values of alpha parameter used during the algorithm.
            The array is not-empty.
            Hyperparameter for hypergraph modularity 
        
        change_mode : string, optional, default = "communities" 
            other options: "iter", "modifications"
            Parameter controling  when the change of alpha is made.
            The default value "communities" is based on observing the current number of communities
            Other options: 
            "modification" - based on the number of modifications (node moves) made.
            "iter" - the option forsing alpha change after each iteration of modularity maximization.

        change_frequency : float, optional, default = 0.5
            The parameter used change_mode = "communities" controlling the frequency of changing alpha.
            E.g. value 0.5 forces the change when number of communities is halved when comparing with the value from the previous change.
        
        after_changes: int, optional, default = 300
            The parameter used change_mode = "modifications" controlling the frequency of changing alpha.
            E.g. value 300 forces the change every 300 modifications (node moves).

        force_alpha_1 : boolean, optional, default = True
            The parameter used when alpha is not equal to 1 at the end of algorithm, 
            e.i., the final maximization function is still not the pure hypergraph modularity (alpha < 1)
            When set to False the algorithm finishes despide the fact that alpha is smaller than 1.
            When set to True the algorithm back to its previous pass and do the additional optimization phase with alpha set to 1.

        random_seed: int

        Returns
        -------
        partition after last step: list of sets
        partition before the last step: list of sets
        h-modularity of partition after the last step: float
        h-modularity of partition before the last step: float
        the sequence of alphas used: list of floats
        
        '''
        A_basic, qH_basic, alphas_out = self.h_louvain_community(alphas, change_frequency=change_frequency, change_mode=change_mode, after_changes=after_changes,
                                                                   force_alpha_1 = force_alpha_1, random_seed=random_seed)
        A_new, qH_new= self.last_step(A_basic, tau = self.hmod_tau, resolution = self.resolution, delta = 0.01) 
        
        return A_new, A_basic, qH_new, qH_basic, alphas_out
    


    def _last_step_unweighted(self, A,  tau = "infinity", resolution = 1, delta=0.01):

        qH = self.h_modularity(A, tau, resolution)
        

        ## calculating wdc
        Ctr = Counter(self.edge_sizes)
        wdc = {}
        for d in Ctr.keys():
            wdc[d]={}
            for c in np.arange(int(np.floor(d / 2 + 1)), d + 1):
                if tau == "infinity":
                    wdc[d][c] = 1 if c == d else 0
                else:
                    wdc[d][c] = (c/d)**tau if c > d / 2 else 0

        
    

        ## initialize
        ctr_sizes = Counter(self.edge_sizes)
        VolA = [sum([self.HG.degree(i) for i in k]) for k in A]
        
        VolV = np.sum(VolA)
        
        dct_A = hmod.part2dict(A)

        while(True):
            
            n_moves = 0
            for v in list(np.random.permutation(list(self.h_nodes))):
                
                dct_A_v = dct_A[v]
                H_id = [self.HG.incidence_dict[x] for x in self.HG.nodes[v]]
                L = [ [dct_A[i] for i in x] for x in H_id ]
                deg_v = self.HG.degree(v)
                
                ## assume unweighted - EC portion before
                _ctr = Counter([ (Counter(l).most_common(1)[0][1],len(l)) for l in L])
                ec = sum([wdc[k[1]][k[0]]*_ctr[k] for k in _ctr.keys() if k[0] > k[1]/2])

                ## DT portion before
                dt = 0
                for d in ctr_sizes.keys():
                    Cnt = ctr_sizes[d]
                    for c in np.arange(int(np.floor(d/2+1)),d+1):
                        dt += (Cnt*wdc[d][c]*binom.pmf(c,d,VolA[dct_A_v]/VolV)) 
                
                ## move it?
                best = dct_A_v
                best_del_q = 0
                best_dt = 0
                for m in set([i for x in L for i in x])-{dct_A_v}:
                    dct_A[v] = m
                    L = [ [dct_A[i] for i in x] for x in H_id ]
                    ## assume unweighted - EC
                    _ctr = Counter([ (Counter(l).most_common(1)[0][1],len(l)) for l in L])
                    ecp = sum([wdc[k[1]][k[0]]*_ctr[k] for k in _ctr.keys() if k[0] > k[1]/2])
                    ## DT
                    del_dt = -dt
                    for d in ctr_sizes.keys():
                        Cnt = ctr_sizes[d]
                        for c in np.arange(int(np.floor(d/2+1)),d+1):
                            del_dt -= (Cnt*wdc[d][c]*binom.pmf(c,d,VolA[m]/VolV))
                            del_dt += (Cnt*wdc[d][c]*binom.pmf(c,d,(VolA[m]+deg_v)/VolV))
                            del_dt += (Cnt*wdc[d][c]*binom.pmf(c,d,(VolA[dct_A_v]-deg_v)/VolV)) 
                    del_q = ecp-ec-(resolution*del_dt)
                    if del_q > best_del_q:
                        best_del_q = del_q
                        best = m
                        best_dt = del_dt
                if best_del_q > 0.1: ## this avoids some numerical precision issues
                    n_moves += 1
                    dct_A[v] = best
                    VolA[m] += deg_v
                    VolA[dct_A_v] -= deg_v
                    VolV = np.sum(VolA)
                else:
                    dct_A[v] = dct_A_v
            new_qH = self.h_modularity(hmod.dict2part(dct_A), tau,resolution)    
            #print(n_moves,'moves, new qH:',new_qH)
            if (new_qH-qH) < delta:
                break
            else:
                qH = new_qH
        return hmod.dict2part(dct_A), new_qH

    ## THIS ASSUMES WEIGHTED H
    def _last_step_weighted(self, A, tau = "infinity", resolution = 1, delta=0.01):

        qH = self.h_modularity(A,tau,resolution)
        Ctr = Counter(self.edge_sizes)
        wdc = {}
        for d in Ctr.keys():
            wdc[d]={}
            for c in np.arange(int(np.floor(d / 2 + 1)), d + 1):
                if tau == "infinity":
                    wdc[d][c] = 1 if c == d else 0
                else:
                    wdc[d][c] = (c/d)**tau if c > d / 2 else 0


        d = hmod.part2dict(A)

        ## initialize
        ## this is the bottleneck
        VolA = np.repeat(0,1+np.max(list(d.values())))
        m = np.max(self.edge_sizes)
        ctr_sizes = np.repeat(0,1+m)
        S = 0
        for e in self.HG.edges:
            w = self.HG.edges[e].weight
            ctr_sizes[self.HG.size(e)] += w  
            S += w
            for v in self.HG.edges[e]:
                VolA[d[v]] += w 
        VolV = np.sum(VolA)
        dct_A = hmod.part2dict(A)

        ## loop
        while(True):
            n_moves = 0
            for v in list(np.random.permutation(list(self.h_nodes))):

                dct_A_v = dct_A[v]
                H_id = [self.HG.incidence_dict[x] for x in self.HG.nodes[v]]
                L = [ [dct_A[i] for i in x] for x in H_id ]

                ## ec portion before move
                _keys = [ (Counter(l).most_common(1)[0][1],len(l)) for l in L]
                _vals = [self.HG.edges.properties["weight"][x] for x in self.HG.nodes[v]]
                _df = pd.DataFrame(zip(_keys,_vals), columns=['key','val'])
                _df = _df.groupby(by='key').sum()
                ec = sum([ wdc[k[1]][k[0]]*val[0] for (k,val) in _df.iterrows() if k[0]>k[1]/2 ])
                str_v = np.sum(_vals) ## weighted degree

                ## DT portion before move
                dt = 0
                for d in range(len(ctr_sizes)):
                    Cnt = ctr_sizes[d]
                    for c in np.arange(int(np.floor(d/2+1)),d+1):
                        if Cnt > 0:
                            dt += (Cnt*wdc[d][c]*binom.pmf(c,d,VolA[dct_A_v]/VolV)) 


                ## move it?
                best = dct_A_v
                best_del_q = 0
                best_dt = 0
                for m in set([i for x in L for i in x])-{dct_A_v}:
                    dct_A[v] = m
                    L = [ [dct_A[i] for i in x] for x in H_id ]
                    ## EC
                    _keys = [ (Counter(l).most_common(1)[0][1],len(l)) for l in L]
                    _vals = [self.HG.edges.properties["weight"][x] for x in self.HG.nodes[v]]
                    _df = pd.DataFrame(zip(_keys,_vals), columns=['key','val'])
                    _df = _df.groupby(by='key').sum()
                    ecp = sum([ wdc[k[1]][k[0]]*val[0] for (k,val) in _df.iterrows() if k[0]>k[1]/2 ])    


                    ## DT
                    del_dt = -dt
                    for d in range(len(ctr_sizes)):
                        Cnt = ctr_sizes[d]
                        for c in np.arange(int(np.floor(d/2+1)),d+1):
                            if Cnt > 0:
                                del_dt -= (Cnt*wdc[d][c]*binom.pmf(c,d,VolA[m]/VolV))
                                del_dt += (Cnt*wdc[d][c]*binom.pmf(c,d,(VolA[m]+str_v)/VolV))
                                del_dt += (Cnt*wdc[d][c]*binom.pmf(c,d,(VolA[dct_A_v]-str_v)/VolV))
                    del_q = ecp-ec-resolution*del_dt
                    if del_q > best_del_q:
                        best_del_q = del_q
                        best = m
                        best_dt = del_dt

                if best_del_q > 0.1: ## this avoids some precision issues
                    n_moves += 1
                    dct_A[v] = best
                    VolA[m] += str_v
                    VolA[dct_A_v] -= str_v
                    VolV = np.sum(VolA)
                else:
                    dct_A[v] = dct_A_v

            new_qH = self.h_modularity(hmod.dict2part(dct_A), tau,resolution)
    #     print(n_moves,'moves, new qH:',new_qH)
            if (new_qH-qH) < delta:
                break
            else:
                qH = new_qH
        return hmod.dict2part(dct_A), new_qH

    def last_step(self, A,  tau = "infinity", resolution = 1, delta=0.01):

        '''
        Parameters
        ----------
        A : list of sets
        some initial partition of the vertices in HG

        tau : str | float, optional, default= "infinity"
            "infinity" or float greater or equal to 0
            Hyperparameter for hypergraph modularity 

        resolution : float, optional, default = 1 (corresponding to the standard settings for modularity calculation) 
            float greater or equal than 0
            Hyperparameter for hypergraph modularity, controlling the weight of degree tax. 
            The classical definition of modularity is retrieved when the resolution parameter is set to 1.

        delta : float, optional
                convergence stopping criterion

        Returns
        -------
        : list of sets
        A new partition for the vertices in HG
        '''
        
        ## all same edge weights?
        #uniq = (len(Counter(self.HG.edges.properties["weight"]))==1)
        uniq = self.uniq

        if uniq:
            nls, qH = self._last_step_unweighted(A, tau , resolution, delta)
        else:
            nls, qH = self._last_step_weighted(A, tau , resolution, delta)
        return nls, qH
     
            


def load_ABCDH_from_file(filename):
    with open(filename,"r") as f:
        rd = csv.reader(f)
        lines = list(rd)

    print("File loaded")

    Edges = []
    for line in lines:
        Edges.append(list(line))

    HG = hnx.Hypergraph(dict(enumerate(Edges)))

    print("HG created")

    print(len(Edges), "edges")
 
    return HG



def main():


    HG = load_ABCDH_from_file("datasets/results_300_he.txt")


    hL = hLouvain(HG,hmod_tau="infinity")

    A, q2, alphas_out = hL.h_louvain_community(alphas = [0.5,1,0.5,0.7,0.8,1], change_mode="communities", change_frequency= 0.5)
        
    print("")
    print("FINAL ANSWER:")
    print("partition:",A)
    print("qC:", q2)
    print("alphas", alphas_out)



if __name__ == "__main__":
    main()
