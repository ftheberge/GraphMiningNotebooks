# Graph Mining Notebooks

Notebooks and datasets to accompany the textbook "Mining Complex Networks" (https://www.torontomu.ca/mining-complex-networks) by B. Kaminski, P. Pralat and F. Théberge.

We are currently working on a **second edition** of this book (as of June 2025). The corresponding notebooks for this new edition are currently located in the **Python_Notebooks_Second_Edition** directory.
There is also a **yml** file that can be used to build the conda environment via: *conda env create -f environment.yml*.
The notebooks and yml files are subject to frequent updates.

Tutorials for the Python notebooks can be found here: https://www.youtube.com/@MiningComplexNetworks

## Other software used:

Chapter 6 notebook: 
 * install the **Julia** language and package **CGE** from: https://github.com/KrainskiL/CGE.jl

Complementary material: 
 * Chapter 8: install and compile the overlapping NMI measure from: https://github.com/aaronmcdaid/Overlapping-NMI
 * Chapter 9: install the **NEExT** framework and depdendencies, details in the notebook
 * Chapter 11: install the **folium** package, details in the notebook
 
## References - datasets

* **ABCD** dataset was generated using the Julia code found at: https://github.com/bkamins/ABCDGraphGenerator.jl

* **Actors** dataset is from the accompanying material to the book: "Complex Networks: Principles, Methods and Applications", V. Latora, V. Nicosia, G. Russo, Cambridge University Press (2017). ISBN: 9781107103184 The data can be downloaded from: http://www.complex-networks.net. It was built from the Internet Movie Database (IMDb, www.imdb.com) and was first studied in: D.J. Watts and S.H. Strogatz, "Collective dynamics of 'small-world' networks", Nature 393 (1998), 440-442.

* **Airport** dataset is available at: https://www.kaggle.com/flashgordon/usa-airport-dataset#Airports2.csv which is part of the Kaggle Public Datasets: https://www.kaggle.com/datasets

* **Football** dataset: is available from http://www-personal.umich.edu/~mejn/netdata/. The source reference is: M. Girvan and M. E. J. Newman,
"Community structure in social and biological networks", Proc. Natl. Acad. Sci. USA 99, 7821-7826 (2002).

* **GitHub Developpers** data is available at https://snap.stanford.edu/data/github-social.html. The source reference is: "Multi-scale Attributed Node Embedding", Benedek Rozemberczki, Carl Allen, and Rik Sarkar, arXiv (2019). https://arxiv.org/abs/1909.13021. This is part of project: https://github.com/benedekrozemberczki/MUSAE. SNAP source reference is: Jure Leskovec and Andrej Krevl, "SNAP Datasets: Stanford Large Network Dataset Collection", http://snap.stanford.edu/data, 2014.
 
* **Game of Thrones** dataset is built from: Jeffrey Lancaster GitHub repository: https://github.com/jeffreylancaster/game-of-thrones

* **Grid** network data is available at http://doi.org/10.5281/zenodo.47317. The source reference is: Wiegmans, B. (2016). GridKit: European and North-American extracts [Data set]. Zenodo. http://doi.org/10.5281/zenodo.47317.

* **NCI1 and NCI109** datasets can be downloaded from https://ls11-www.cs.tu-dortmund.de/staff/morris/graphkerneldatasets. Website reference: Kristian Kersting, Nils M. Kriege, Christopher Morris, Petra Mutzel and Marion Neumann, "Benchmark Data Sets for Graph Kernels", 2016. The source reference is: N. Wale and G. Karypis. Comparison of descriptor spaces for chemical compound retrieval and classification. In Proc. of ICDM, pages 678–689, Hong Kong, 2006.

* The source graph for the **Reno** road network example was prepared using the map taken from: https://github.com/pszufe/OpenStreetMapX.jl/blob/master/test/data/reno_east3.osm

* **Words** free-association data was built from the database: http://w3.usf.edu/FreeAssociation/ Reference: Nelson, D. L., McEvoy, C. L., & Schreiber, T. A. (1998). The University of South Florida word association, rhyme, and word fragment norms.  
An XML version is available here: http://rali.iro.umontreal.ca/rali/?q=en/USF-FAN  
Details for using this database are provided in: Douglas L. Nelson, Cathy L. McEvoy and Thomas A. Schreiber, "The University of South Florida free association, rhyme, and word fragment norms", Behavior Research Methods, Instruments, & Computers 2004, 36 (3), 402–407 https://link.springer.com/content/pdf/10.3758/BF03195588.pdf  
Our illustration is based on the work in: G. Palla, I. Derényi, I.J. Farkas, T. Vicsek, "Uncovering the overlapping community structure of complex networks in nature and society", Nature 435(7043):814-818, July 2005.

* **Zachary** dataset is available within *igraph*. The source reference is: Zachary, W. W. (1977). "An Information Flow Model for Conflict and Fission in Small Groups". Journal of Anthropological Research. 33 (4): 452–473. JSTOR 3629752.

## Links -- code

* most examples use python-igraph which is available from: https://igraph.org/python/

* ensemble clustering is available at: https://pypi.org/project/partition-igraph/  
Ref: Valérie Poulin and François Théberge, "Ensemble clustering for graphs: comparisons and applications", Appl Netw Sci 4, 51 (2019). https://doi.org/10.1007/s41109-019-0162-z

* graph2vec code can be obtained from: https://github.com/benedekrozemberczki/graph2vec   
Ref: graph2vec: Learning distributed representations of graphs. Narayanan, Annamalai and Chandramohan, Mahinthan and Venkatesan, Rajasekar and Chen, Lihui and Liu, Yang MLG 2017, 13th International Workshop on Mining and Learning with Graphs (MLGWorkshop 2017).

* node2vec is included in SNAP and can be obtained at: https://snap.stanford.edu/node2vec/   
Ref: node2vec: Scalable Feature Learning for Networks. A. Grover, J. Leskovec. ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD), 2016.

* UMAP code is available at: https://pypi.org/project/umap-learn/  
Ref: McInnes, L., Healy, J. UMAP: uniform manifold approximation and projection for dimension reduction. Preprint at https://arxiv.org/abs/1802.03426 (2018).

* Normalized Mutual Information (NMI) measure for sets of overlapping clusters: https://github.com/aaronmcdaid/Overlapping-NMI  
Ref: Normalized Mutual Information to evaluate overlapping community finding algorithms" by Aaron F. McDaid, Derek Greene, Neil Hurley http://arxiv.org/abs/1110.2515

* Some graph embedding functions for igraph we wrote are based on the networkx code found in the GEM package: https://github.com/palash1992/GEM

* HyperNetX (HNX) can be found at: https://github.com/pnnl/HyperNetX including several good tutorial notebooks.
 
### First edition

The notebooks in the **Python_Notebooks_First_Edition** directory were tested under the following conda environment (Python version 3.10.9, igraph version 0.10.4).
It is no longer maintained.

```
conda create --name graphmining numpy pandas jupyter matplotlib scikit-learn statsmodels seaborn cairo pycairo bokeh cython gensim numba datashader holoviews colorcet networkx

conda activate graphmining

pip install igraph
pip install plfit
pip install partition-igraph
pip install umap-learn
pip install --no-dependencies graphrole
pip install celluloid
pip install --no-dependencies hypernetx
```
