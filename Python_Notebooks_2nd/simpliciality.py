import numpy as np
import random
from math import comb
from collections import Counter
from collections.abc import Sequence
import operator
from functools import reduce 
from itertools import combinations as combs

################################
## max subsets from: https://stackoverflow.com/questions/14106121
def is_power_of_two(n):
    return (n & (n - 1)) == 0
def max_subsets(E):
    '''return list of max subsets of length 3+'''
    if len(E) <= 1:
        return list(E)
    if not isinstance(E, Sequence):
        E = list(E)
    incidence = {}
    for i, s in enumerate(E):
        for element in s:
            try:
                incidence[element] |= 1 << i
            except KeyError:
                incidence[element] = 1 << i
    out = [s for s in E if s and len(s)>2 and 
           is_power_of_two(reduce(operator.and_, (incidence[x] for x in s)))]
    return out
################################


################################
### Landry's measures
# These three functions are all from "The simpliciality of higher order networks":
# simplicial_fraction, edit_simpliciality and face_edit_simpliciality
################################


def simplicial_fraction_jordan(V,E):
    # They ignore 1-edges in their paper, so we throw away 1-edges here
    E = [set(e) for e in E if len(e) > 1]
    edge_sets = get_edge_sets(V,E)

    # For each edge, we add 1 to bottom and add 1 to top if the edge is a simplicial complex
    top = 0
    bottom = 0
    for e in E:
        #In their paper, the only check edges of size at least 3
        if len(e) > 2:
            bottom += 1
            #relevant_edges are edges containing v for each v in e
            relevant_edges = set()
            for v in e:
                relevant_edges = relevant_edges.union([frozenset(sorted(edge)) for edge in edge_sets[v]])
            #po_set is the simplicial closure of e. We don't include sets of size 0 or 1, nor do we include e
            po_set = set()
            for k in range(2,len(e)):
                po_set = po_set.union([frozenset(sorted(edge)) for edge in combs(e,k)])

            #Check if e is a simplicial complex
            if (po_set <= relevant_edges):
                top += 1
                
    return top/bottom
            
def get_simplicial_fraction(V,E):
    '''
    Get the simplicial fraction of a hypergraph
    The simplicial fraction is the fraction of edges (of size at least 3) that satisfy downward closure
    V: vertices
    E: edges
    '''
    # They ignore 1-edges in their paper, so we throw away 1-edges here
    E = [set(e) for e in E if len(e) > 1]
    edge_sets = get_edge_sets(V,E)

    # For each edge, we add 1 to bottom and add 1 to top if the edge is a simplicial complex
    top = 0
    bottom = 0
    for e in E:
        # In their paper, the only check edges of size at least 3
        if len(e) > 2:
            bottom += 1
            # relevant_edges are edges containing v for each v in e
            relevant_edges = set()
            for v in e:
                relevant_edges = relevant_edges.union([tuple(sorted(edge)) for edge in edge_sets[v]])
            # po_set is the simplicial closure of e. We don't include sets of size 0 or 1, nor do we include e
            po_set = set()
            for k in range(2,len(e)):
                po_set = po_set.union([tuple(sorted(edge)) for edge in combs(e,k)])

            # Check if e is a simplicial complex
            if (po_set <= relevant_edges):
                top += 1                
    return top/bottom

def edit_simpliciality_jordan(V,E,return_max=False):
    # They ignore 1-edges in their paper, so we throw away 1-edges here
    E = [set(e) for e in E if len(e) > 1]    

    # Finding the set of maximal edges of size at least 3
    E_check = [e for e in E if len(e) > 2]
    edge_sets = get_edge_sets(V,E_check)
    # We build a dictionary of e -> T/F
    is_max = {frozenset(sorted(e)) : True for e in E_check}
    for e in E_check:
        # Some edges will point to false before they are hit in the loop
        if is_max[frozenset(sorted(e))]:
            # Any edge containing e will contain v for all v in e
            # Hence, we find v in e with the fewest edges to check for edges containing e
            v = next(iter(e))
            for u in e:
                if len(edge_sets[u]) < len(edge_sets[v]):
                    v = u
            # For each edge, we check containment in both directions and update accordingly
            for edge in edge_sets[v]:
                if set(e) < set(edge):
                    is_max[frozenset(sorted(e))] == False
                    break
                elif set(edge) < set(e):
                    is_max[frozenset(sorted(edge))] == False
            
    E_max = [e for e in E_check if is_max[frozenset(sorted(e))]]

    # C is the set of edges in the simplicial closure.
    # We start by adding all 2-edges, then add all PoSets of maximal edges
    C = {frozenset(sorted(e)) for e in E if len(e) == 2}
    for e in E_max:
        for k in range(2,len(e)+1):
            C = C.union({frozenset(sorted(edge)) for edge in combs(e,k)})
    if return_max:
        return E_max, len(E)/len(C)
    else:
        return len(E)/len(C)
     
def get_edit_simpliciality(V,E,return_max=False):
    '''
    Get the edit simpliciality of a hypergraph
    The edit simpliciality is |E|/|C| where C is the smallest simplicial complex containing E
    V: vertices
    E: edges
    '''
    # throw away 1-edges
    E = [set(e) for e in E if len(e) > 1]    
    # Finding the set of maximal edges of size at least 3
    E_max = max_subsets(E)
    C = [tuple(sorted(e)) for e in E if len(e)==2]
    for e in E_max:
        for k in range(2,len(e)+1):
            C.extend([tuple(sorted(x)) for x in combs(e,k)])
    C = set(C)
    if return_max:
        return E_max, len(E)/len(C)
    else:
        return len(E)/len(C)
            

#### Warning - E_max returns too many sets ####
def face_edit_simpliciality_jordan(V,E,return_max=False):
    # They ignore 1-edges in their paper, so we throw away 1-edges here
    E = [set(e) for e in E if len(e) > 1]    

    #Finding the set of maximal edges of size at least 3
    E_check = [e for e in E if len(e) > 2]
    edge_sets = get_edge_sets(V,E_check)
    #We build a dictionary of e -> T/F
    is_max = {frozenset(sorted(e)) : True for e in E_check}
    for e in E_check:
        #Some edges will point to false before they are hit in the loop
        if is_max[frozenset(sorted(e))]:
            #Any edge containing e will contain v for all v in e
            #Hence, we find v in e with the fewest edges to check for edges containing e
            v = next(iter(e))
            for u in e:
                if len(edge_sets[u]) < len(edge_sets[v]):
                    v = u
            #For each edge, we check containment in both directions and update accordingly
            for edge in edge_sets[v]:
                if set(e) < set(edge):
                    is_max[frozenset(sorted(e))] == False
                    break
                elif set(edge) < set(e):
                    is_max[frozenset(sorted(edge))] == False
            
    E_max = [e for e in E_check if is_max[frozenset(sorted(e))]]
    
    #We iterate over E_max and compute the edit simpliciality of each maximal face
    FES = 0
    #We now need edge_sets for all edges
    edge_sets = get_edge_sets(V,E)

    #We loop through E_max and compute the edit simpliciality for each edge
    for e in E_max:
        #C_face is the simplicial closure of e
        C_face = {frozenset(sorted(e))}
        for k in range(2,len(e)):
            C_face = C_face.union({frozenset(sorted(edge)) for edge in combs(e,k)})
        #Before computing E_face, we restrict to only edges containing v for each v in e
        relevant_edges = {frozenset(sorted(e))}
        for v in e:
            relevant_edges = relevant_edges.union({frozenset(sorted(edge)) for edge in edge_sets[v]})
        #E_face is the set of edges in C_face that are also in the graph
        E_face = C_face.intersection(relevant_edges)
        FES += len(E_face)/len(C_face)        
    FES = FES/len(E_max)
    if return_max:
        return E_max, FES
    else:
        return FES

def get_face_edit_simpliciality(V,E,return_max=False):
    '''
    Get the face edit simpliciality of a hypergraph
    The face edit simpliciality is the average edit simpliciality across the maximal edges
    V: vertices
    E: edges
    '''
    # Throw away 1-edges here
    E = [set(e) for e in E if len(e) > 1]    
    E_max = max_subsets(E)
    # We iterate over E_max and compute the edit simpliciality of each maximal face
    FES = 0
    # We now need edge_sets for all edges
    edge_sets = get_edge_sets(V,E)

    # We loop through E_max and compute the edit simpliciality for each edge
    for e in E_max:
        C_face = [tuple(sorted(e))]
        for k in range(2,len(e)):
            C_face.extend([tuple(sorted(x)) for x in combs(e,k)])
        C_face = set(C_face)                       
        # Before computing E_face, we restrict to only edges containing v for each v in e
        relevant_edges = {tuple(sorted(e))}
        for v in e:
            relevant_edges = relevant_edges.union({tuple(sorted(edge)) for edge in edge_sets[v]})
        # E_face is the set of edges in C_face that are also in the graph
        E_face = C_face.intersection(relevant_edges)
        FES += len(E_face)/len(C_face)       
    FES = FES/len(E_max)
    if return_max:
        return E_max, FES
    else:
        return FES

######################################################################################################

# Get a dictionary of vertex -> edges containing vertex.
def get_edge_sets(V,E):
    # Initalize edge_set
    edge_sets = {v : [] for v in V}
    #iterate through E and update edge_sets
    for e in E:
        for v in e:
            edge_sets[v].append(e)

    return edge_sets

def get_simplicial_pairs(vertices, edges, as_matrix=False, edge_order=False):
    V = vertices.copy()
    E = edges.copy()
    max_size = max({len(e) for e in E})
    E = [set(e) for e in E]

    ## Getting relevant data
    edge_sets = get_edge_sets(V,E)
    degrees = {v : len(edge_sets[v]) for v in V}

    # Initializing
    if as_matrix:
        M = np.zeros((max_size,max_size))
    elif edge_order:
        M = [0,0]
    else:
        M = 0

    # We iterate through E and compute the number of pairs with e as the smaller edge
    for e in E:
        # Any edge f containing e will contain v for all v in e
        # Thus, we find v with the smallest degree for our search
        v = list(e)[0]
        for u in e:
            if degrees[u] < degrees[v]:
                v = u
        relevant_edges = edge_sets[v]
        # Now we check for f containing e
        # If order matters, we compare the index of e and f in edge_sets[v]
        # edge_sets[v] preserves edge order, so all is good
        if edge_order:
            # Note: we assume there are no multi-edges in E
            # This is precisely why we can't order edges in the Chung-Lu model
            i = relevant_edges.index(e)
            for j,f in enumerate(relevant_edges):
                if e<f:
                    if as_matrix:
                        if i<j:
                            M[len(e)-1][len(f)-1] += 1
                        else:
                            M[len(f)-1][len(e)-1] += 1
                    else:
                        if i<j:
                            M[0] += 1
                        else:
                            M[1] += 1

        # If order doesn't matter, we just count
        else:
            for f in relevant_edges:
                if e<f:
                    if as_matrix:
                        M[len(e)-1][len(f)-1] += 1
                    else:
                        M += 1

    return M

def get_simplicial_ratio(V,E,samples = 1, edge_order = False, multisets = False):
    '''
    Get the simplicial ratio of a hypergraph:
        number of simplicial pairs divided by the expected number of pairs
    Expected number estimated via sampling Chung-Lu hypergraphs (with or without repeated edges)
    V: vertices
    E: edges
    '''
    top = get_simplicial_pairs(V,E,edge_order=edge_order)
    bottom = 0
    for i in range(samples):
        E_cl = chung_lu(V,E,multisets=multisets)
        bottom += get_simplicial_pairs(V,E_cl,edge_order=False)

    if edge_order:
        bottom = [(1+bottom)/(2*samples),(1+bottom)/(2*samples)]
        return [top[0]/bottom[0],top[1]/bottom[1]]
    else:
        bottom = (1+bottom)/samples
        return top/bottom


def get_simplicial_matrix(V,E,samples = 1, edge_order = False, multisets=False):
    '''
    Get the simplicial ratio of a hypergraph for each pair of edge sizes:
        number of simplicial pairs divided by the expected number of pairs
    Expected number estimated via sampling Chung-Lu hypergraphs (with or without repeated edges)
    V: vertices
    E: edges
    '''
    # We first get the matrix of simplicial pairs
    M_top = get_simplicial_pairs(V,E,as_matrix=True,edge_order=edge_order)

    #Next, we get the same matrix for the re-sample, depending on the model specified
    E_resampled = chung_lu(V,E,multisets=multisets)
    #We never order edges in the Chung-Lu model
    #Instead, we split the results in two equal parts
    M_bottom = get_simplicial_pairs(V,E_resampled,as_matrix=True,edge_order=False)
    for i in range(samples-1):
        E_resampled = chung_lu(V,E,multisets=multisets)
        temp = get_simplicial_pairs(V,E_resampled,as_matrix=True,edge_order=False)
        M_bottom = np.add(M_bottom, temp)
    M_bottom = (M_bottom + 1)/samples

    #Now we check for edge_order
    if edge_order:
        for i in range(len(M_bottom)):
            for j in range(i+1,len(M_bottom)):
                M_bottom[i][j] = M_bottom[i][j]/2
                M_bottom[j][i] = M_bottom[i][j]

    return M_top/M_bottom

def chung_lu(vertices,m,degrees = None, multisets = False):
    # If V is a number, we convert to a list
    if type(vertices) is int:
        V = list(range(vertices))
    else:
        V = list(vertices.copy())

    #If m is a list of edges, we get the degrees and the number of edges of each size
    if (type(m[0]) is set) or (type(m[0]) is list):
        degrees = get_degrees(V,m.copy())

        #sort edges by size
        edge_dict = sort_edges(m.copy())

        #convert m to a list of num edges per size
        sizes = edge_dict.keys()
        m = [0]*max(sizes)
        for k in sizes:
            m[k-1] = len(edge_dict[k])

    #If degrees is a list, we convert to a dictionary
    #Note: If V was given as a set, and degrees as a list of degrees, then the degrees might get shuffled
    if type(degrees) is list:
        degrees = dict(zip(V,degrees))
    L = get_volume(V,degrees)

    #choices is a dictionary with degrees[v] keys pointing to v
    #I've tested, and this is much faster than making a list
    choices = dict.fromkeys(set(range(L)))
    counter = 0
    current_vertex = 0
    #We need L keys in total
    while (counter<L):
        for i in range(degrees[V[current_vertex]]):
            choices[counter] = V[current_vertex]
            counter += 1
        current_vertex += 1

    #E is the set of edges to be returned
    E = []
    if multisets:
        for k in range(len(m)):
            #Adding all edges of size k+1
            for i in range(m[k]):
                e = []
                for j in range(k+1):
                    e.append(choices[random.randint(0,L-1)])
                E.append(e)
    else:
        for k in range(len(m)):
            #Adding all edges of size k+1
            for i in range(m[k]):
                e = []
                while len(e)<k+1:
                    v = (choices[random.randint(0,L-1)])
                    if v not in e:
                        e.append(v)
                E.append(e)
    return E

def get_degrees(V,E):
    #Initalize degrees
    degrees = {v : 0 for v in V}

    #Update degrees by iterating through E
    for e in E:
        for v in e:
            degrees[v] = degrees[v]+1

    return degrees

def get_volume(S,X):
    #initialize volume
    volume = 0

    #get degrees if X is a list of edges
    if type(X) is list:
        #We need all vertices in the graph to call get_degrees
        V = S.copy()
        for e in X:
            V.update(e)
        X = get_degrees(V,X)

    #iterate through V and sum degrees
    for v in S:
        volume += X[v]

    return volume

def sort_edges(E):
    #initialize partition
    edge_sizes = {len(e) for e in E}
    partition = {k : [] for k in edge_sizes}

    #iterate through E and update partition
    for e in E:
        partition[len(e)].append(e)

    return partition
