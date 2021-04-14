# Graph Mining Notebooks

Notebooks and datasets to accompany the textbook "Mining Complex Networks" by B. Kaminski, P. Pralat and F. Théberge.

## Software environment

The Python Notebooks were created under the following conda environment:

```
conda create --name graphmining python=3.7.9 numpy pandas jupyter matplotlib scikit-learn statsmodels seaborn cairo pycairo bokeh cython gensim numba datashader holoviews colorcet

conda activate graphmining

pip install python-igraph==0.8.2
pip install plfit
pip install partition-igraph
pip install umap-learn
pip install --no-dependencies graphrole
pip install hypernetx
```

Other software used:

Chapter 5 notebook: 
 * install the julia language from https://julialang.org
 * ABCD generator from https://github.com/bkamins/ABCDGraphGenerator.jl (more details in the notebook)

Chapter 6 notebook: 
 * install node2vec from SNAP, see: https://snap.stanford.edu/node2vec/
 * GED code from https://github.com/ftheberge/Comparing_Graph_Embeddings which is included in this repo (more details in the notebook)

Complementary material: 
 * install the graph2vec code from https://github.com/benedekrozemberczki/graph2vec
 * LFR benchmark with overlapping communities, file binary_networks.tar, from: https://sites.google.com/site/andrealancichinetti/files
 * overlapping NMI measure from: https://github.com/aaronmcdaid/Overlapping-NMI
 * pip install omega-index-py3
 * conda install -c conda-forge folium
 
## References - datasets

* **ABCD** dataset was generated using the Julia code found at: https://github.com/bkamins/ABCDGraphGenerator.jl

* **Actors** dataset is from the accompanying material to the book: "Complex Networks: Principles, Methods and Applications", V. Latora, V. Nicosia, G. Russo, Cambridge University Press (2017). ISBN: 9781107103184 The data can be downloaded from: http://www.complex-networks.net. It was built from the Internet Movie Database (IMDb, www.imdb.com) and was first studied in: D.J. Watts and S.H. Strogatz, "Collective dynamics of 'small-world' networks", Nature 393 (1998), 440-442.

* **Airport** dataset is available at: https://www.kaggle.com/flashgordon/usa-airport-dataset#Airports2.csv which is part of the Kaggle Public Datasets: https://www.kaggle.com/datasets

* **Football** dataset: is available from http://www-personal.umich.edu/~mejn/netdata/. The source reference is: M. Girvan and M. E. J. Newman,
"Community structure in social and biological networks", Proc. Natl. Acad. Sci. USA 99, 7821-7826 (2002).

* **GitHub Developpers** data is available at https://snap.stanford.edu/data/github-social.html. The source reference is: "Multi-scale Attributed Node Embedding", Benedek Rozemberczki, Carl Allen, and Rik Sarkar, arXi (2019). https://arxiv.org/abs/1909.13021. This is part of project: https://github.com/benedekrozemberczki/MUSAE. SNAP source reference is: Jure Leskovec and Andrej Krevl, "SNAP Datasets: Stanford Large Network Dataset Collection", http://snap.stanford.edu/data, 2014.
 
* **Game of Thrones** dataset is built from: Jeffrey Lancaster GitHub repository: https://github.com/jeffreylancaster/game-of-thrones

* **Grid** network data is available at http://doi.org/10.5281/zenodo.47317. The source reference is: Wiegmans, B. (2016). GridKit: European and North-American extracts [Data set]. Zenodo. http://doi.org/10.5281/zenodo.47317.

* **NCI1 and NCI109** datasets can be downloaded from https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets. Website reference: Kristian Kersting, Nils M. Kriege, Christopher Morris, Petra Mutzel and Marion Neumann, "Benchmark Data Sets for Graph Kernels", 2016. The source reference is: N. Wale and G. Karypis. Comparison of descriptor spaces for chemical compound retrieval and classification. In Proc. of ICDM, pages 678–689, Hong Kong, 2006.

* The source graph for the **Reno** road network example was prepared using the map taken from: https://github.com/pszufe/OpenStreetMapX.jl/blob/master/test/data/reno_east3.osm

* **Words** free-association data was built from the database: http://w3.usf.edu/FreeAssociation/ Reference: Nelson, D. L., McEvoy, C. L., & Schreiber, T. A. (1998). The University of South Florida word association, rhyme, and word fragment norms. Details for using this database are provided in: Douglas L. Nelson, Cathy L. McEvoy and Thomas A. Schreiber, "The University of South Florida free association, rhyme, and word fragment norms", Behavior Research Methods, Instruments, & Computers 2004, 36 (3), 402–407 https://link.springer.com/content/pdf/10.3758/BF03195588.pdf

* **Zachary** dataset is available within *igraph*. The source reference is: Zachary, W. W. (1977). "An Information Flow Model for Conflict and Fission in Small Groups". Journal of Anthropological Research. 33 (4): 452–473. JSTOR 3629752.

## Links -- code

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
 
