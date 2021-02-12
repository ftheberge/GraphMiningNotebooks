# Graph Mining Notebooks

Notebooks and datasets to accompany the textbook "Mining Complex Networks" by B. Kaminski, P. Pralat and F. Theberge.

## Software environment

The Python Notebooks were created under the following conda environment:

```
conda create --name graphmining python=3.7 numpy pandas jupyter matplotlib scikit-learn statsmodels seaborn cairo pycairo bokeh cython gensim numba datashader holoviews colorcet

conda activate graphmining

pip install python-igraph
pip install plfit
pip install partition-igraph
pip install umap-learn
pip install graphrole
pip install hypernetx
```

Other software used:

Chapter 5 notebook: 
 * install the julia language from https://julialang.org
 * ABCD generator from https://github.com/bkamins/ABCDGraphGenerator.jl (more details in the notebook)

Chapter 6 notebook: 
 * install node2vec from https://snap.stanford.edu/node2vec/
 * GED code from https://github.com/ftheberge/Comparing_Graph_Embeddings (more details in the notebook)

Complementary material: 
 * install the graph2vec code from https://github.com/benedekrozemberczki/graph2vec
 * LFR benchmark with overlapping communities, file binary_networks.tar, from: https://sites.google.com/site/andrealancichinetti/files
 * overlapping NMI measure from: https://github.com/aaronmcdaid/Overlapping-NMI
 * pip install omega-index-py3
 * conda install -c conda-forge folium
 
## References - datasets

* ABCD dataset was generated using the Julia code found at: https://github.com/bkamins/ABCDGraphGenerator.jl

* Actors dataset is from the accompanying material to the book: "Complex Networks: Principles, Methods and Applications", V. Latora, V. Nicosia, G. Russo, Cambridge University Press (2017). ISBN: 9781107103184 The data can be downloaded from: http://www.complex-networks.net

* Airport dataset is available at: https://www.kaggle.com/flashgordon/usa-airport-dataset#Airports2.csv

* GitHubDeveloppers: data is available at https://snap.stanford.edu/data/github-social.html
The source reference is:
```
@misc{rozemberczki2019multiscale,    
       title = {Multi-scale Attributed Node Embedding},   
       author = {Benedek Rozemberczki and Carl Allen and Rik Sarkar},   
       year = {2019},   
       eprint = {1909.13021},  
       archivePrefix = {arXiv},  
       primaryClass = {cs.LG}   
}
```
which is part of project: https://github.com/benedekrozemberczki/MUSAE

* Power grid network data is available at: https://zenodo.org/record/47317#.XzwNDi0ZNhF

* NCI1 and NCI109 datasets can be downloaded from https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets
The source reference is: N. Wale and G. Karypis. Comparison of descriptor spaces for chemical compound retrieval and 
classification. In Proc. of ICDM, pages 678–689, Hong Kong, 2006.

* The Zachary dataset is available within igraph. The source reference is: Zachary, W. W. (1977). "An Information Flow Model for Conflict and Fission in Small Groups". Journal of Anthropological Research. 33 (4): 452–473. JSTOR 3629752.

* Words data build from: http://w3.usf.edu/FreeAssociation/
Reference: Nelson, D. L., McEvoy, C. L., & Schreiber, T. A. (1998). The University of South Florida word association, rhyme, and word fragment norms. http://www.usf.edu/FreeAssociation/.

* Game of Thrones dataset built from: https://github.com/jeffreylancaster/game-of-thrones

* The source graph for the road network example was prepared using the map taken from: https://github.com/pszufe/OpenStreetMapX.jl/blob/master/test/data/reno_east3.osm



## References -- code

* most examples use python-igraph with is available from: https://igraph.org/python/

* ensemble clustering is available at: https://pypi.org/project/partition-igraph/

* graph2vec code can be obtained from: https://github.com/benedekrozemberczki/graph2vec

* node2vec is included in SNAP and can be obtained at: https://snap.stanford.edu/node2vec/

* plfit code is avaialble from: https://pypi.org/project/plfit/

* graphrole code is available from: https://pypi.org/project/graphrole/

* UMAP code is available at: https://pypi.org/project/umap-learn/

* Some graph embedding functions for igraph we wrote are based on the networkx code found in the GEM package: https://github.com/palash1992/GEM

* omega measure: https://pypi.org/project/omega-index-py3/

* LFR with overlapping communities: https://sites.google.com/site/andrealancichinetti/files/

* Normalized Mutual Information (NMI) measure for sets of overlapping clusters: https://github.com/aaronmcdaid/Overlapping-NMI
Reference: Normalized Mutual Information to evaluate overlapping community finding algorithms" by Aaron F. McDaid, Derek Greene, Neil Hurley http://arxiv.org/abs/1110.2515

 * HyperNetX (HNX) can be found at: https://github.com/pnnl/HyperNetX including several good tutorial notebooks.
 
