## Degree correlation functions
import numpy as np
from sklearn.linear_model import LinearRegression
from collections import Counter

# undirected
def deg_corr(G):
    idx = {v: k for k, v in enumerate([i for i in set(G.degree())])}
    idx_inv = {k: v for k, v in enumerate(idx)}
    deg = G.degree()
    L = [[] for i in range(len(idx))]
    for v in G.vs():
        w = [deg[i] for i in G.neighbors(v)]
        L[idx[v.degree()]].extend(w)
    return {idx_inv[i]: np.mean(L[i]) for i in range(len(L))}


## k_nn^{mode1,mode2}(l) : average mode2-degree of out-neighbours of nodes with mode1-degree = l
def deg_corr_directed(G, mode1="all", mode2="all"):
    idx = {v: k for k, v in enumerate([i for i in set(G.degree(mode=mode1))])}
    idx_inv = {k: v for k, v in enumerate(idx)}
    deg = G.degree(mode=mode2)
    L = [[] for i in range(len(idx))]
    for v in G.vs():
        w = [deg[i] for i in G.neighbors(v, mode="out")]  ## do each link only once
        L[idx[v.degree(mode=mode1)]].extend(w)
    return {idx_inv[i]: np.mean(L[i]) for i in range(len(L)) if len(L[i]) > 0}


## degree correlation for neutral undirected graph with degree distribution in G
def deg_corr_neutral(G, mode="all"):
    x = G.degree(mode=mode)
    return np.mean([i**2 for i in x]) / np.mean(x)

## Correlation exponent via linear regression (taking logs)
def corr_exp(G):
    ## compute knn's
    knn = deg_corr(G)
    # Fit the regression
    regressor = LinearRegression()
    x = [np.log(i) for i in knn.keys()]
    y = [np.log(i) for i in knn.values()]
    regressor.fit(np.array(x).reshape(-1, 1), y)
    return regressor.coef_[0]

def corr_exp_directed(G, mode1="all", mode2="all"):
    ## compute knn's
    knn = deg_corr_directed(G, mode1, mode2)
    # Fit the regression
    regressor = LinearRegression()
    A = np.array(
        [[np.log(k), np.log(v)] for k, v in zip(knn.keys(), knn.values()) if k > 0 and v > 0]
    )
    # x = [np.log(i) for i in knn.keys()]
    # y = [np.log(i) for i in knn.values()]
    regressor.fit(A[:, 0].reshape(-1, 1), A[:, 1])
    return regressor.coef_[0]


## for a fixed value for l
def rich_club(g, l=1):
    g.vs()["degree"] = g.degree()
    l_max = np.max(g.degree())
    c = Counter(g.degree())
    n = g.vcount()
    moment = np.sum([k * c[k] / n for k in c]) ** 2
    S = [k * c[k] / n for k in c if k >= l]
    phi_hat = np.sum([x * y for x in S for y in S]) * g.ecount() / moment
    G = g.subgraph([v for v in g.vs() if v["degree"] >= l])
    phi = G.ecount()
    return phi / phi_hat


## Asssortativity for directed graphs is only defined for the ('out','in') case in igraph
## Here is a general function for all cases
def assortativity_directed(G, mode1="out", mode2="in"):
    d1 = G.degree(mode=mode1)
    d2 = G.degree(mode=mode2)
    D = []
    for e in G.es:
        D.append([d1[e.source], d2[e.target]])
    D = np.array(D)
    avg1 = np.mean(D[:, 0])
    avg2 = np.mean(D[:, 1])
    S1 = sum([(v[0] - avg1) * (v[1] - avg2) for v in D])
    S2 = np.sqrt(sum([(v[0] - avg1) * (v[0] - avg1) for v in D]))
    S3 = np.sqrt(sum([(v[1] - avg2) * (v[1] - avg2) for v in D]))
    return S1 / (S2 * S3)


## same as igraph fnction
def assortativity_undirected(G):
    d1 = G.degree(mode="all")
    D = []
    for e in G.es:
        D.append([d1[e.source], d1[e.target]])
        D.append([d1[e.target], d1[e.source]])
    D = np.array(D)
    avg = np.mean(D[:, 0])
    S1 = sum([(v[0] - avg) * (v[1] - avg) for v in D])
    S2 = np.sqrt(sum([(v[0] - avg) * (v[0] - avg) for v in D]))
    S3 = np.sqrt(sum([(v[1] - avg) * (v[1] - avg) for v in D]))
    return S1 / (S2 * S3)