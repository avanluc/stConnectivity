# S2M
Master thesis project


The presented algorithm pretend to resolve the ST-CONNECTIVITY problem on undirected graphs performing multiple BFS and then connecting the sources.
The basic idea is to perform several BFS from different sources and continue until a treshold of edges is reached. During the visits a matrix of adjacence is created and it represents the graph formed by all the BFS sources. Whenever is found a connection between two BFS the matrix is updated. 
At the and of this first phase we have visited a part of the input graph corrisponding to the Treshold. So a parallel visit of the matrix is performed looking for a connectivity of original source and target of the problem. If the connectivity result is true than the algorithm stops and print the computation time, otherwise the algorithm starts phase two.
The second phase goal is to complete the visit of the graph running backward, so it is a Bottom-Up visit. Every non visited vertex looks for a neighboor already visited and, if it exist, the vertex connects to it updating the adjacence matrix. At the end of this visit every vertex of the graph is visited and we can perform a second parallel visit of the matrix looking for connectivity. If the result is still false it means that there isn't a path from source to target.



