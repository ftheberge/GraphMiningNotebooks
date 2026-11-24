## unit square torus distance
def torusDist(U, V):
    dx = np.abs(V[0] - U[0])
    dy = np.abs(V[1] - U[1])
    if dx > 0.5:
        dx = 1 - dx
    if dy > 0.5:
        dy = 1 - dy
    return np.sqrt(dx * dx + dy * dy)
## functions for the SPA and XBS models

import numpy as np
import igraph as ig

# unit square (non-torus) distance
def squareDist(U, V):
    dx = np.abs(V[0] - U[0])
    dy = np.abs(V[1] - U[1])
    return np.sqrt(dx * dx + dy * dy)


# build a SPA graph (in 2-d) given the parameters described above
def buildSPA(n, p, A1, A2, seed=12345):
    ## random points
    rng = np.random.default_rng(seed=seed)
    X = rng.random(n)
    Y = rng.random(n)
    pts = [(X[i], Y[i]) for i in range(n)]

    ## initialize new graph
    t = 0
    G = ig.Graph(directed=True)

    ## add vertices
    for v in range(n):
        G.add_vertex(str(v))
        G.vs[v]["loc"] = pts[v]
        for i in range(v):
            if torusDist(G.vs[v]["loc"], G.vs[i]["loc"]) < G.vs[i]["r"]:
                if rng.random() < p:
                    G.add_edge(v, i)
        t += 1
        ## re-compute radii of the spheres of influence (here in 2-d)
        G.vs["r"] = [np.sqrt(min(1, (A1 * d + A2) / t) / np.pi) for d in G.degree(mode="in")]
    return G

## Naive Xulvi-Brunet and Sokolov algorithm
def XBS(n, d, q, assortative):
    p = d / n
    g = ig.Graph.Erdos_Renyi(n, p)
    g.es["touched"] = False
    ec = g.ecount()
    while True:
        re = np.random.choice(ec, 2, replace=False)
        nodes = list(g.es[re[0]].tuple + g.es[re[1]].tuple)
        if len(set(nodes)) == 4:
            ## with proba q, wire w.r.t. assortativity, else randomly
            if np.random.random() < q:
                idx = np.argsort(g.degree(nodes))
                if assortative:
                    e1 = (nodes[idx[0]], nodes[idx[1]])
                    e2 = (nodes[idx[2]], nodes[idx[3]])
                else:
                    e1 = (nodes[idx[0]], nodes[idx[3]])
                    e2 = (nodes[idx[1]], nodes[idx[2]])
            else:
                np.random.shuffle(nodes)
                e1 = (nodes[0], nodes[1])
                e2 = (nodes[2], nodes[3])
            if (
                g.get_eid(e1[0], e1[1], directed=False, error=False)
                + g.get_eid(e2[0], e2[1], directed=False, error=False)
                == -2
            ):
                g.delete_edges(re)
                g.add_edge(e1[0], e1[1], touched=True)
                g.add_edge(e2[0], e2[1], touched=True)
            else:
                g.es[re[0]]["touched"] = True
                g.es[re[1]]["touched"] = True
        if sum(g.es["touched"]) == g.ecount():
            break
    return g