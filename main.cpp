#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <random>
#include <limits>
using namespace std;
using Graph = vector<vector<int>>;

// -------------------- Read Graph --------------------
void read_graph(Graph &g, int &V, int &E) {
    if (!(cin >> V >> E)) {
        cerr << "Failed to read V E\n";
        exit(1);
    }
    g.assign(V, {});
    for (int i = 0; i < E; i++) {
        int u, v;
        cin >> u >> v;
        if (u>=0 && v>=0 && u<V && v<V) {
            g[u].push_back(v);
            g[v].push_back(u);
        }
    }
}
int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int V,E; Graph G;
    read_graph(G,V,E);

    auto bc_exact = brandes_full(G);
    auto bc_cluster = cluster_based_bc_louvain(G);

    cout<<"=== Exact Brandes BC ===\n";
    print_vector(bc_exact);

    cout<<"=== Cluster-based BC (Improved) ===\n";
    print_vector(bc_cluster);

    auto R1 = rank_nodes(bc_exact);
    auto R2 = rank_nodes(bc_cluster);

    int K=min(10,V);
    cout<<"=== Top-"<<K<<" Exact ===\n";
    for(int i=0;i<K;i++) cout<<R1[i].second<<": "<<R1[i].first<<"\n";

    cout<<"=== Top-"<<K<<" Cluster-based ===\n";
    for(int i=0;i<K;i++) cout<<R2[i].second<<": "<<R2[i].first<<"\n";

    return 0;
}
