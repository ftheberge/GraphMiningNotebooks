Data-Set: Movie actors network
Description: This data set (see Box 2.2 on page 81 of the book)
             contains the co-starring network studied for the first
             time by Duncan Watts and Steven Strogratz in their
             seminal paper "Collective dynamics of small-world
             networks". Each node is an actor, and a link between two
             actors  exists if they have acted togeter in at least one movie. 
Author: Duncan Watts, Steven Strogatz
Reference: D. Watts, S. Strogatz, "Collective behaviour of small-world
           networks", Nature 393, 440-442 (1998)
Additional-Info: This data set is part of the accompanying material of the book
                 "Complex Networks: Principles, Methods and Applications", 
                 V. Latora, V. Nicosia, G. Russo, Cambridge University Press (2017)
                 ISBN: 9781107103184
                 If you use the data set, please add a reference to the book above.
Book-Reference: Chapter 2, Box 2.2

Filename: movie_actors.net
Type: edge_list
Nodes: 248243 
Edges: 8302734
Directed: no
Weighted: no
File-Format: "actor_1 actor__2"

Filename: movie_actors_names_dists
Type: node labels and distances
Rows: 248243
File-Format: "node_label: actor_name d0 d1 d2 d3..."
             - 'node_label' is the label of the node
             - 'actor_name' is the name of the actor
             - 'd0', 'd1', 'd2' etc. are the number of nodes at
               distance 0,1,2, etc. from 'node_label'
