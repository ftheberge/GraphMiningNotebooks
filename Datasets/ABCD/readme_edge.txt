# Common parameters

n = "1500"                    # number of vertices in graph
t1 = "2"                      # power-law exponent for degree distribution
d_min = "1"                   # minimum degree
d_max = "50"                  # maximum degree
t2 = "2"                      # power-law exponent for cluster size distribution
c_min = "50"                  # minimum cluster size
c_max = "1000"                # maximum cluster size
isCL = "false"                # if "false" use configuration model, if "true" use Chung-Lu

# fraction of edges to fall in background graph

edge01.dat: xi=0.1
edge04.dat: xi=0.4
edge07.dat: xi=0.7
edge10.dat: xi=1.0
