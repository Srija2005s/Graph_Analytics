# Graph_Analytics
## Efficient Graph Algorithm analysis for Betweenness Centrality
### Introduction:
Betweenness Centrality (BC) is a fundamental graph metric that measures how often a node lies on shortest paths between other nodes. This helps quantify how important or influential a node is in the graph. Importance may relate to shortest travel routes, fraud detection, communication bottlenecks, information flow, point of contact, or popularity ranking. It is widely used in social networks, communication networks, biological systems and web ranking. Various algorithms are designed to calculate the betweenness centrality on large graphs, the algorithms explored are:
1. Brandes' algorithm by Ulrik Brande
2. Modularity-based clustering(E1C-FastBC algorithm - modified brandes)

This project includes:
1. Understanding theoretical BC and Brandes algorithm
2. Understanding clustering-based FastBC concepts
3. Designing and analysing C++ implementations of the algorithms

The motivation is that Brandes’ BC is accurate but extremely slow for large graphs
(time complexity O(V × E)). Hence, E1C-FastBC inspired workflow is considered, where clustering, boundary nodes, and local/global dependency decomposition help scale BC computation. In this project, gradual versions of cluster-based BC are designed, each improving accuracy over earlier attempts.

### Betweenness Centrality:
                          BC(v) = ∑(s!=v!=t) (σst​(v)​​/σst)
σst​(v) = number of those shortest paths that pass through node v
​​σst = number of shortest paths from s to t

### Naive Approach:
1. For every pair of nodes (s,t), computes all shortest paths between s and t.
2. store all shortest paths and count of the paths that include v
3. for every node v!=s,t, we check whether v lies on each of the shortest paths
4. then update the BC by adding the ratio of number of paths with v in it to total number of paths(σst​(v)​​/σst).
5. this is repeated for all the pairs

#### Time complexity:
Unweighted graphs: O(V^3)
Weighted graphs: O(nm+n^2logn)

This is not feasible for larger graphs.
### Brandes workflow(Exact Betweenness Centrality):
Standard betweenness centrality is calculated on various graphs depending on kind of graphs -- Unweighted graphs use BFS, weighted graphs use Dijkstra's -- the steps includes:
1. Runs BFS/Dijksta's from every source
2. Computes shortest-path counts (σ) and predecessors
3. Performs backward dependency accumulation (δ) and the dependecies of every traversal are added with their contribution to global betweenness centrality to get overall betweeness centrality of the nodes.

The overall time complexity of this algorithm is O(V*E) for unweighted graphs, which is still infeasible for large networks.
To address the scalability, cluster-based BC computation is addressed.

### Cluster-based BC workflow:
Clustering is a technique used to divide a large graph into smaller, densely connected groups of nodes called clusters. In the context of Betweenness Centrality computation, clustering is a fundamental idea used in algorithms like E1C-FastBC. 

Clustering is essential to E1C-FastBC because it:

1. Reduces the size of subgraphs where BC is computed
2. Separates local and global contributions
3. Identifies "boundary nodes" that lie between communities
4. Enables building a compressed cluster graph

Louvain-like clustering is adapted in this analysis.

1.Initial Clustering(Louvain-like):
    Clustering is done by starting all nodes in their own community with population of one each. Using modularity gain:
   
                   gain = (k_i_in/2m) - (dv*total)/2m^2
   
        k_i_in = number of edges from node i to all other nodes in same community
        dv = degree of node i in whole graph
        total = sum of degrees of all nodes in community
        m = total number of edges in the graph
this is repeated until there is no more change, usually it's not repeated more than 10 times
2.Compress cluster IDs:
    1.Assign compact IDs to each discovered cluster.
    2.Create cluster → member node list.
3.Identify Boundary Nodes:
    Boundary nodes are nodes whose neighbors belong to different clusters.
    These nodes are critical because:
    1. All inter-cluster shortest paths must pass through them
    2. They influence global betweenness centrality heavily
4.Mapping Members & Boundary Nodes:
    Store:
      1.Internal nodes of each cluster
      2.Boundary nodes shared across clusters
    These are used as sources for local BC.
5.Local Betweenness Centrality Computation: (Local_BC)
    Local BC handles paths fully inside its own cluster.
    Depending on the implementation:
        1.Version 1 (merged_bc.cpp):
            Run Brandes-from-subset using all internal nodes + boundary nodes.
        2.Version 2 (boundary-only):
            Run Brandes using only boundary nodes.
        3.Version 3 (boundary + important internal nodes):
            Choose top-K high-degree internal nodes + boundary nodes.
    Local BC captures intra-cluster importance.
6.Build the Cluster Graph:
    1.Each cluster becomes a super-node.
    2.Edges between clusters appear if any original edge connected them.
    3.Duplicate edges are removed.
  This compressed graph preserves inter-cluster connectivity.
7. Compute Global Betweenness Centrality:(Global_BC)
    Run Brandes on the cluster graph (which is very small).
    This yields BC for each cluster (global scores).
8. Distribute Global BC to Nodes
    Distribute cluster-level BC to individual nodes using degree-proportional assignment:
          BC(v)+=BC(cluster(v))×(degree(v)/∑u∈cluster(v)degree(u))
    Nodes with higher degree receive more global importance.
9. Combine Local & Global BC:
    Final BC is:
          BC(v) = Local_BC(v) + Global_BC(v)
    This provides a fast and accurate approximation of exact BC.

Compare cluster-based BC output with exact BC by comparing top nodes. 

### Local BC Algorithm Variants:
In this project, three different Local_BC computation strategies were implemented, each trading off accuracy and efficiency differently. Since the Local BC component heavily influences the final BC score, understanding these variants is crucial.

Version 1: Full Local BC (merged_bc.cpp)
  1.Uses all internal nodes + all boundary nodes as Brandes sources.
  2.Highest accuracy but slowest.
  3.Best for dense clusters or social networks.

Version 2: Boundary-Only Local BC (merged_12th_boundary.cpp)
  1.Uses only boundary nodes as sources.
  2.Fastest but can miss internal cluster structure.
  3.Best for sparse or tree-like clusters.

Version 3: Hybrid Local BC (merged_boundaryAndImp.cpp)
  1.Uses boundary nodes + top-K important internal nodes (high-degree nodes).
  2.Balanced accuracy and speed.
  3.Best for real-world networks with hubs.

### Experimental Evaluation:
To evaluate the effectiveness of cluster-based BC, the outputs of the three Local BC variants were compared with exact Brandes BC using top-K ranked nodes. Results show that Version 1 closely matches exact BC on dense clusters, Version 2 performs well on sparse and tree-like graphs, and Version 3 achieves the best balance between accuracy and computation cost on real-world graphs. This confirms that clustering preserves important nodes while significantly improving scalability.

### Conclusion:
This project shows that while Brandes’ algorithm provides exact Betweenness Centrality, it does not scale to large graphs. By using Louvain clustering, boundary node detection, and local–global BC decomposition, the cluster-based approach significantly reduces computation while preserving the ranking of important nodes. Among the three Local BC variations, the hybrid method (boundary + important internal nodes) offers the best balance between speed and accuracy. Overall, clustering proves to be an effective strategy for scalable BC analysis on large real-world networks.
