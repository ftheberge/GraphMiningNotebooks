from itertools import combinations as combs
import random
import numpy as np

# Get a dictionary of vertex -> degree.
def degrees(V,E):
    # Initalize degrees
    degrees = {v : 0 for v in V}

    # Update degrees by iterating through E
    for e in E:
        for v in e:
            degrees[v] = degrees[v]+1
    
    return degrees

# Get dictionary of size -> edges of that size
def sort_edges(E):
    # initialize partition
    edge_sizes = {len(e) for e in E}
    partition = {k : [] for k in edge_sizes}
    
    # iterate through E and update partition
    for e in E:
        partition[len(e)].append(e)
        
    return partition

# Get the connected components as a partition of vertices (list of sets)
# The graph is converted to a union-find structure for the sake of speed
def components(vertices,edges):

    V = vertices.copy()
    E = edges.copy()
    #We start by converting to a union-find dataset
    parents = {v : v for v in V}
    sizes = {v : 1 for v in V}
    
    #Inernal function: Finds the current root of a vertex
    def root(v):
        if parents[v] == v:
            return v
        else:
            return root(parents[v])

    #Internal function: Merges sets
    def merge(S):  
        #replace elements by their roots
        S = {root(v) for v in S}

        #If all vertices were in the same component, we don't need to merge
        if len(S) > 1:
            #When we merge, we root using the largest tree
            v = next(iter(S))
            for u in S:
                if sizes[u] > sizes[v]:
                    v = u

            #We now point all of S to v
            for u in S:
                parents[u] = v

            #Lastly, we update the size of v
            #Note: We only care about the sizes of roots, so there is no need to update the size of S\{v}
            #Note: We converted S to a set of roots, so this sum is always correct
            sizes[v] = sum([sizes[u] for u in S])

    #Iterate through E and merge
    for e in E:
        merge(e)

    #Build a dictionary of root -> component
    roots = {root(v) for v in V}
    components = {v : {v} for v in roots}
    for v in V:
        components[root(v)].add(v)

    return list(components.values())

## vertices and edges in the giant component
def giant_component(vertices, edges):
    V = vertices.copy()
    E = edges.copy()
    comps = components(V,E)
    C_max = comps[0]
    for C in comps:
        if len(C) > len(C_max):
            C_max = C
    #We pick a vertex v in e and keep e iff v is in C_max
    #This does not distrupt the order of E
    E = [e for e in E if list(e)[0] in C_max]
    V = C_max
    V = list(V)
    return V, E
 
def giant_component_growth(V, E, cutoff=None, shuffle_edges=True):
    #We start by converting to a union-find dataset
    parents = {v : v for v in V}
    sizes = {v : 1 for v in V}
    
    #Inernal function: Finds the current root of a vertex
    def root(v):
        if parents[v] == v:
            return v
        else:
            return root(parents[v])

    #Internal function: Merges sets
    def merge(S):  
        #replace elements by their roots
        S = {root(v) for v in S}

        #If all vertices were in the same component, we don't need to merge
        if len(S) > 1:
            #When we merge, we root using the largest tree
            v = next(iter(S))
            for u in S:
                if sizes[u] > sizes[v]:
                    v = u

            #We now point all of S to v
            for u in S:
                parents[u] = v

            #Lastly, we update the size of v
            #Note: We only care about the sizes of roots, so there is no need to update the size of S\{v}
            #Note: We converted S to a set of roots, so this sum is always correct
            sizes[v] = sum([sizes[u] for u in S])

    #We shuffle E (unless shuffle_edges=False), add e in E one by one, and keep track of the size of the largest component along the way
    if shuffle_edges:
        random.shuffle(E)
    size_of_giant = {0:1}
    counter = 0
    if cutoff is not None:
        for e in E[:cutoff]:
            counter += 1
            merge(e)
            size_of_giant[counter] = max(size_of_giant[counter-1],sizes[root(next(iter(e)))])
    else:
        for e in E:
            counter += 1
            merge(e)
            size_of_giant[counter] = max(size_of_giant[counter-1],sizes[root(next(iter(e)))])
            
    return size_of_giant

def proc_variables(d, m, q):
    
    #basic fixing
    if type(d) is list:
        d = {v:d[v] for v in range(len(d))}
    #If d is an int, it's assumed to be the desired number of vertices
    #In this case, we give all vertices the same degree
    elif type(d) is int:
        d = {v:1 for v in range(d)}
    V = list(d.keys())
    if type(m) is list:
        m = {(k+1):m[k] for k in range(len(m)) if m[k]>0}
    if (type(q) is float) or (type(q) is int):
        q = {k:q for k in m.keys()}

    #if d and m don't line up, we scale d 
    m_tot = sum(k*m[k] for k in m.keys())
    L = sum(d.values())
    if L != m_tot:
        scale = m_tot/L
        d = {v:round(scale*d[v]) for v in V}
        while sum(d.values()) < m_tot:
            d[np.random.choice(V)] += 1
        while sum(d.values()) > m_tot:
            v = np.random.choice(V)
            if d[v] > 0:
                d[v] -= 1

    #getting the rest of the variables
    L = sum(d.values())
    p = [d[v]/L for v in V]
    E = {k:[] for k in m.keys()}
    
    return d, m, q, V, E, L, p

def check_vertex_list(vertex_list, max_edge_size, V, L, p):
    if len(vertex_list) < max_edge_size:
        return get_vertex_list(V, L, p, False)
    else:
        return vertex_list

def make_normal_edge(k, vertex_list, multisets):
    if multisets:
        e = vertex_list[-k:]
        vertex_list = vertex_list[:-k]
    else:
        e = []
        while len(e)<k:
            v = vertex_list.pop()
            if v not in e:
                e.append(v)
    return e, vertex_list

def make_simplicial_edge(k, E, m, vertex_list, multisets):
    #sample_list is the list of edges NOT of size k
    sample_list = []
    for j in m.keys():
        if j != k:
            sample_list.extend(E[j])
    #if there are no edges to pair with, we build a normal Chung-Lu edge instead
    if len(sample_list) == 0:
        return make_normal_edge(k, vertex_list, multisets)
    else:
        sample = sample_list[np.random.randint(len(sample_list))]
        r = k-len(sample)
        #If the sampled edge is smaller, we build e by joining the sampled edge with vertices from vertex_list
        if r>0:
            e = sample.copy()
            if multisets:
                e.extend(vertex_list[-r:])
                vertex_list = vertex_list[:-r]
            else:
                while len(e)<k:
                    v = vertex_list.pop()
                    if v not in e:
                        e.append(v)
        #Otherwise, we build e by taking a uniform k-subset of e
        else:
            e = list(random.sample(sample,k))
        return e, vertex_list

def get_vertex_list(V, L, p, multisets):
    if multisets:
        return list(np.random.choice(V, size = L, p = p))
    #If multisets are not allowed, we often throw away choices
    #To compensate, we make vertex_list way bigger
    else:
        return list(np.random.choice(V, size = 10*L, p = p))

def get_size_list(m):
    size_list = []
    for k in m.keys():
        size_list.extend([k]*m[k])
    random.shuffle(size_list)
    return size_list

def root(v, parents):
    if parents[v] == v:
        return v
    else:
        return root(parents[v], parents)

def merge(S, parents, sizes):  
    S = {root(v, parents) for v in S}
    if len(S) > 1:
        v = next(iter(S))
        for u in S:
            if sizes[u] > sizes[v]:
                v = u
        for u in S:
            parents[u] = v
        sizes[v] = sum([sizes[u] for u in S])
    return parents, sizes

def normalize(L):
    if type(L) is list:
        norm = sum(L)
        return [x/norm for x in L]
    elif type(L) is dict:
        norm = sum(L.values())
        return {k:v/norm for k,v in L.items()}
    print('invalid input for normalization')
    return

def get_components(V, parents):
    components = {}
    for v in V:
        rv = root(v, parents)
        if rv in components.keys():
            components[rv].append(v)
        else:
            components[rv] = [v]
    return components

def connect_some(k, d, p, components, sizes, multisets):
    e = []
    root_weights = normalize([sizes[v] for v in components.keys()])
    merge_roots = np.random.choice(list(components.keys()), size = k, p = root_weights, replace = False)
    for v in merge_roots:
        weights = normalize({u:d[u] for u in components[v]})
        u = np.random.choice(components[v], p = list(weights.values()))
        e.append(u)
    return e

def connect_all(V, k, d, p, components, multisets):
    e = []
    for v in components.keys():
        weights = normalize({u:d[u] for u in components[v]})
        u = np.random.choice(components[v], p = list(weights.values()))
        e.append(u)
    excess = k-len(components)
    if excess > 0:
        if multisets:
            e_remaining = list(np.random.choice(V, size = excess, p = p, replace = True))
            e.extend(e_remaining)
        else:
            while len(e) < k:
                v_try = np.random.choice(V, p = p)
                if v_try not in e:
                    e.append(v_try)
    return e
 
def make_skeleton(V, d, m, p, size_list, multisets):
    
    parents = {v : v for v in V}
    sizes = {v:p[i] for i, v in enumerate(V)}
    E = []
    components = {v:[v] for v in V}
    
    while True:
        k = size_list.pop()
        if len(components) < k:
            e = connect_all(V, k, d, p, components, multisets)
            E.append(e)
            break
        else:            
            e = connect_some(k, d, p, components, sizes, multisets)
            E.append(e)
            parents, sizes = merge(e, parents, sizes)            
            components = get_components(V, parents)
    return E, size_list

# d: degree sequence
# m: dict of edge sizes
# q: simplicial edge probability
def simplicial_chung_lu(d, m, q, multisets = True, skeleton = False):
    
    d,m,q,V,E,L,p = proc_variables(d,m,q)
    vertex_list = get_vertex_list(V,L,p,multisets)
    size_list = get_size_list(m)

    # We can optionally create a connected skeleton before continuing
    if skeleton:
        E_all, size_list = make_skeleton(V, d, m, p, size_list, multisets = multisets)
    else:
        E_all = []
    
    # Iterate through size_list and build edges
    for k in size_list:
        make_simplicial = (np.random.uniform() < q[k])
        if make_simplicial:
            e, vertex_list = make_simplicial_edge(k, E, m, vertex_list, multisets)
        # Normal chung-lu edges are built here
        else:
            e, vertex_list = make_normal_edge(k, vertex_list, multisets)
        E[k].append(e)
        if not multisets:
            vertex_list = check_vertex_list(vertex_list, max(m.keys()), V, L, p)
        
    for E_k in E.values():
        E_all.extend(E_k)
    return E_all
