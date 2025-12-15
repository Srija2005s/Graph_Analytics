
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

// -------------------- Brandes (Exact) --------------------
vector<double> brandes_full(const Graph &G) {
    int V = G.size();
    vector<double> BC(V,0.0);

    for(int s=0; s<V; s++){
        stack<int> S;
        vector<vector<int>> P(V);
        vector<int> d(V,-1);
        vector<double> sigma(V,0.0);

        d[s]=0; sigma[s]=1.0;
        queue<int>q; q.push(s);

        while(!q.empty()){
            int v=q.front(); q.pop();
            S.push(v);
            for(int w:G[v]){
                if(d[w]<0){
                    d[w]=d[v]+1;
                    q.push(w);
                }
                if(d[w]==d[v]+1){
                    sigma[w]+=sigma[v];
                    P[w].push_back(v);
                }
            }
        }

        vector<double> delta(V,0.0);
        while(!S.empty()){
            int w=S.top(); S.pop();
            for(int v:P[w]){
                if(sigma[w]!=0)
                    delta[v] += (sigma[v]/sigma[w])*(1.0+delta[w]);
            }
            if(w!=s) BC[w]+=delta[w];
        }
    }
    return BC;
}
