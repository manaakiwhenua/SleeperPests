library(igraph) #Garrett used this package in the INA package so lets use it also for the social networks.

###IMPORT OR GENERATE A SOCIAL NETWORK

#OPTION 1: import social network data

#Norman's code includes the following line:
#SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread

#convert the matrix to a graph object
mtx <- SEAM
SEAM1 <-graph_from_adjacency_matrix(mtx) #convert the matrix to an iGraph object

#OPTION 2: construct a random social network (Erdős–Rényi model) based on a given probability of land-owners being connected to each other.

n <- #number of nodes in the network, i.e. how many land-owners (nodes) you want to include
p <- #the probability for drawing a link (connection) between two random nodes  
#OR
m <- #the number of links in the network
# in the "type", choose type = c("gnp" or "gnm") based on the choice between p and m  

SEAM2 <- erdos.renyi.game( 
  n, 
  p.or.m, 
  type = c("gnp", "gnm"), 
  directed = TRUE,
  loops = FALSE
)

#example:
SEAM2 <- erdos.renyi.game( #igraph package
  30, #number of nodes in the network, i.e. how many land-owners (nodes) you want to include
  0.05, #either the probability for drawing a link (connection) between two random nodes, or the number of links in the network
  type = c("gnp"), #see  line above
  directed = TRUE,
  loops = FALSE
)


#OPTION 3: small-world model for a social network
dim <- #the dimension of the starting lattice, use 1
size <- #number of land-owners in the network  
nei <- #the neighborhood within which the links will be connected
p <- #the rewiring probability
  
SEAM3 <- sample_smallworld(dim, size, nei, p, loops = FALSE, multiple = FALSE)

#example: 
SEAM3 <- sample_smallworld(1, 10, 3, 0.1)


#### NETWORK ROLES FOR TARGETING SPECIFIC LAND-OWNERS

graph <- SEAM2

#land-owners with most direct connections to other land-owners:
V(graph)$degree <- degree(graph) #calculate network degree for each land-owner (how many other land-owners each land-owner is connected to) and store the degree as node attribute.
#Those with degree = 0 are so called network isolates, i.e. not connected to any other land-owners.

#land-owners with strongest connections
E(graph)$weight <- seq(0, 1, length = ecount(graph)) #give random weights for the network since we don't have data for it
V(graph)$w_degree <- strength(graph, mode = "out") #mode = out refers to out-going connections, i.e. land-owners influence or share knowledge to others. Replace with mode = "in" to
#find land-owners with strongest incoming connections

#land-owners in brokerage roles, i.e. connect otherwise unconnected land-owners
V(graph)$betweenness <- betweenness(graph)






