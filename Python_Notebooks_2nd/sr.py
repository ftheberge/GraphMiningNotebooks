import numpy as np
import random
import math
from itertools import combinations
from collections import Counter
from time import time
import xgi
TRACK_TIME = False

def simplicial_ratio(V, E, num_pre_process_samples = 1000, num_edge_samples = 1000):

    t0 = time()
    up, down, M, M_ordered = simplicial_pairs(V, E)
    t1 = time()
    if TRACK_TIME:
        print('time to count simplicial pairs = {0}'.format(t1-t0))
    
    M_cl = expected_pairs(V, E, num_pre_process_samples, num_edge_samples)
    pairs_cl = np.sum(M_cl)

    ratio = (up+down)/pairs_cl
    up_ratio = 2*up/pairs_cl
    down_ratio = 2*down/pairs_cl
    M_ratio = M.copy()
    M_ordered_ratio = M_ordered.copy()

    sizes = {len(e) for e in E}
    for k in range(len(M)):
        for l in range(len(M[0])):
            if M_cl[k][l] != 0:
                M_ratio[k][l] = M[k][l]/M_cl[k][l]
                M_ordered_ratio[k][l] = 2*M_ordered[k][l]/M_cl[k][l]
                M_ordered_ratio[l][k] = 2*M_ordered[l][k]/M_cl[k][l]

    return ratio, up_ratio, down_ratio, M_ratio, M_ordered_ratio

def probability_of_multiset(deg, num, edge_size):
    top = 0
    samples = get_samples(num*edge_size, deg)
    for i in range(num):
        edge = samples[i*edge_size: (i+1)*edge_size]
        if len(edge) != len(set(edge)):
            top += 1
    return top/num

def pre_process(edge_sizes, deg, num):
    p_multiset = {}
    for k in edge_sizes:
        p = probability_of_multiset(deg, num, k)
        p_multiset[k] = p
    return p_multiset

def refine_pre_process(deg, p_old, num_old, num_new):
    edge_sizes = list(p_old.keys())
    p_new = pre_process(edge_sizes, deg, num_new)
    w_old = num_old/(num_old + num_new)
    w_new = num_new/(num_old + num_new)
    p_multiset = {k: w_old*p_old[k] + w_new*p_new[k] for k in edge_sizes}
    return p_multiset

def probability_of_pair(edge, k, deg, volume, p_multiset):
    edge_degs = [deg[v] for v in edge]
    k_set_weights = [math.prod(t) for t in combinations(edge_degs, k)]
    edge_weight = math.factorial(k)*sum(k_set_weights)
    total_weight = (1-p_multiset[k])*(volume**k)
    return edge_weight/total_weight

def expected_pairs(V, E, num_pre_process_samples, num_edge_samples):

    t0 = time()
    E_sorted = sort_edges(E)
    m = {k: len(E_sorted[k]) for k in E_sorted}
    
    E = {e: E[e] for e in range(len(E))}
    es = edge_sets(V, E)
    deg = {v: len(es[v]) for v in V}
    
    edge_sizes = list(m.keys())
    k_min = min(edge_sizes)
    k_max = max(edge_sizes)
    volume = sum(deg.values())
    t1 = time()
    if TRACK_TIME:
        print('time to build variables = {0}'.format(t1-t0))
    
    p_multiset = pre_process(edge_sizes, deg, num_pre_process_samples)
    t2 = time()
    if TRACK_TIME:
        print('time to pre-process probabilities = {0}'.format(t2-t1))
    
    M = np.zeros((k_max+1, k_max+1))

    sampled_edges = {l: get_chung_lu_edges(num_edge_samples, l, deg) for l in edge_sizes if l != k_min}
    t3 = time()
    if TRACK_TIME:
        print('time to sample chung-lu edges = {0}'.format(t3-t2))
    for k in edge_sizes:
        for l in edge_sizes:
            if k<l:
                p_pair = 0
                for edge in sampled_edges[l]:
                    p_pair += probability_of_pair(edge, k, deg, volume, p_multiset)
                p_pair = p_pair/len(sampled_edges[l])
                num_expected = m[k]*m[l]*p_pair
                M[k][l] = num_expected
    t4 = time()
    if TRACK_TIME:
        print('time to compute expected matrix = {0}'.format(t4-t3))
    return M

def get_samples(num, deg):
    p = normalize(list(deg.values()))
    return np.random.choice(a=list(deg.keys()), size=num, p=p)

def get_next_edge(num, deg, size, samples, pointer, max_pointer):
    edge = samples[pointer: pointer+size]
    while len(edge) != len(set(edge)):
        pointer += size
        if pointer >= max_pointer:
            samples = get_samples(num*size*size, deg)
            pointer = 0
        edge = samples[pointer: pointer+size]
    pointer += size
    if pointer >= max_pointer:
        samples = get_samples(num*size*size, deg)
        pointer = 0
    return edge, samples, pointer

def get_chung_lu_edges(num, size, deg):
    samples = get_samples(num*size*size, deg)
    edges = [[] for i in range(num)]
    pointer = 0
    max_pointer = num*size
    for i in range(num):
        edge, samples, pointer = get_next_edge(num, deg, size, samples, pointer, max_pointer)
        edges[i].extend(edge)
    return edges

def update_pairs(E, e, E_relevant, up, down, M, M_ordered):
    for f in E_relevant:
        if set(E[e]) < set(E[f]):
            M[len(E[e])][len(E[f])] += 1
            if e < f:
                up += 1
                M_ordered[len(E[e])][len(E[f])] += 1
            else:
                down += 1
                M_ordered[len(E[f])][len(E[e])] += 1   
    return up, down, M, M_ordered

def update_list_of_pairs(E, e, E_relevant, pairs):
    for f in E_relevant:
        if set(E[e]) < set(E[f]):
            pairs.append([e, f]) 
    return pairs


def simplicial_pairs(V, E):
    if type(E) is list:
        E = {e: E[e] for e in range(len(E))}
    E_with = edge_sets(V, E)
    deg = {v : len(E_with[v]) for v in V}
    max_size = max(len(E[e]) for e in E)

    up = 0
    down = 0
    M = np.zeros((max_size+1, max_size+1))
    M_ordered = np.zeros((max_size+1, max_size+1))

    for e in E:
        if len(E[e]) < max_size:
            v = min_deg_in_edge(E[e], deg)
            E_relevant = E_with[v]
            up, down, M, M_ordered = update_pairs(E, e, E_relevant, up, down, M, M_ordered)
    return up, down, M, M_ordered


def list_of_simplicial_pairs(V, E):
    if type(E) is list:
        E = {e: E[e] for e in range(len(E))}
    E_with = edge_sets(V, E)
    deg = {v : len(E_with[v]) for v in V}
    max_size = max(len(E[e]) for e in E)
    
    pairs = []

    for e in E:
        if len(E[e]) < max_size:
            v = min_deg_in_edge(E[e], deg)
            E_relevant = E_with[v]
            pairs = update_list_of_pairs(E, e, E_relevant, pairs)
    return pairs

def edge_sets(V, E):
    E_with = {v: [] for v in V}
    for e in E:
        for v in E[e]:
            E_with[v].append(e)
    return E_with


def min_deg_in_edge(edge, deg):
    edge_deg = {v: deg[v] for v in edge}
    v = min(edge_deg, key = edge_deg.get)
    return v


def normalize(X):
    norm = sum(X)
    return [x/norm for x in X]


def rand_round(x):
    r = x - math.floor(x)
    if random.uniform(0,1) < r:
        return x+1
    else:
        return x


def sort_edges(E):
    sizes = {len(e) for e in E}
    E_sorted = {k: [] for k in sizes}
    for e in E:
        E_sorted[len(e)].append(e)
    return E_sorted


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


def components(V, E):
    parents = {v : v for v in V}
    sizes = {v : 1 for v in V}
    for e in E:
        parents, sizes = merge(e, parents, sizes)
    roots = {root(v, parents) for v in V}
    components = {v : {v} for v in roots}
    for v in V:
        components[root(v, parents)].add(v)
    return list(components.values())


def giant_component(V, E):
    comps = components(V, E)
    C_max = comps[0]
    for C in comps:
        if len(C) > len(C_max):
            C_max = C
    E = [e for e in E if list(e)[0] in C_max]
    V = C_max
    V = list(V)
    return V, E


def make_matrix_pretty(M):
    M = np.delete(M, 0, axis=0)
    M = np.delete(M, 0, axis=0)
    M = np.delete(M, 0, axis=1)
    M = np.delete(M, 0, axis=1)
    for i in range(len(M)):
        for j in range(len(M[i])):
            if M[i][j] > 1000:
                M[i][j] = -1
            elif M[i][j] >= 100:
                M[i][j] = round(M[i][j])
            else:
                M[i][j] = round(M[i][j], 1)
    return M

