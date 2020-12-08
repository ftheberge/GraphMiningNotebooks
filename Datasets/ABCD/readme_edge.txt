# Common parameters

seed = "42"                   # RNG seed, use "" for no seeding
n = "2000"                    # number of vertices in graph
t1 = "2.5"                    # power-law exponent for degree distribution
d_min = "2"                   # minimum degree
d_max = "100"                 # maximum degree
t2 = "2.5"                    # power-law exponent for cluster size distribution
c_min = "50"                  # minimum cluster size
c_max = "1000"                # maximum cluster size
isCL = "false"                # if "false" use configuration model, if "true" use Chung-Lu
seed = "1234"

# fraction of edges to fall in background graph

edge00001.dat: xi=0.0001
edge001.dat: xi=0.01
edge01.dat: xi=0.1
edge1.dat: xi=1.0
