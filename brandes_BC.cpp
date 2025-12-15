
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

// Brandes but restricted source set
vector<double> brandes_from_sources_subset(const Graph &G, const vector<int> &sources) {
    int V = G.size();
    vector<double> BC(V,0.0);
    for(int s : sources){
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
                    delta[v]+=(sigma[v]/sigma[w])*(1.0+delta[w]);
            }
            if(w!=s) BC[w]+=delta[w];
        }
    }
    return BC;
}

// -------------------- Louvain-Like Clustering --------------------
struct Louvain {
    const Graph &G;
    int V;
    long long m2;
    vector<int> com, deg, tot;

    Louvain(const Graph &g):G(g),V(g.size()) {
        deg.assign(V,0);
        for(int i=0;i<V;i++) deg[i]=G[i].size();
        m2=0; for(int d:deg) m2+=d;
        com.resize(V);
        for(int i=0;i<V;i++) com[i]=i;
        tot=deg;
    }

    long long node_ki_in(int i, int c){
        long long cnt=0;
        for(int nb:G[i]) if(com[nb]==c) cnt++;
        return cnt;
    }

    bool one_level_move(){
        bool moved_any=false;
        vector<int> order(V);
        iota(order.begin(),order.end(),0);
        for(int v:order){
            int orig=com[v];
            long long dv=deg[v];
            tot[orig]-=dv;

            unordered_map<int,int> neigh;
            for(int nb:G[v]) neigh[ com[nb] ]++;

            int best = orig;
            double bestGain=0;

            for(auto &kv:neigh){
                int c=kv.first;
                long long k_i_in=kv.second;
                double gain = (double)k_i_in - (double)dv*(double)tot[c]/(double)m2;
                if(gain>bestGain){
                    bestGain=gain;
                    best=c;
                }
            }
            com[v]=best;
            tot[best]+=dv;
            if(best!=orig) moved_any=true;
        }
        return moved_any;
    }

    vector<int> run(int maxlv=5){
        bool moved=one_level_move();
        if(!moved) return com;
        // compress IDs
        unordered_map<int,int> mp; int nid=0;
        vector<int> out(V);
        for(int i=0;i<V;i++){
            int c=com[i];
            if(!mp.count(c)) mp[c]=nid++;
            out[i]=mp[c];
        }
        return out;
    }
};

// -------------------- Boundary Node Detection --------------------
vector<int> find_boundary_nodes(const Graph &G, const vector<int> &cluster) {
    int V = G.size();
    vector<int> boundary;
    vector<bool> mark(V,false);
    for(int u=0;u<V;u++){
        for(int v:G[u]) if(cluster[u]!=cluster[v]){
            mark[u]=true;
            break;
        }
    }
    for(int i=0;i<V;i++) if(mark[i]) boundary.push_back(i);
    return boundary;
}

// -------------------- Cluster Graph Construction --------------------
unordered_map<int, vector<int>> group_members(const vector<int> &cluster){
    unordered_map<int,vector<int>> mem;
    for(int i=0;i<cluster.size();i++)
        mem[ cluster[i] ].push_back(i);
    return mem;
}

Graph build_cluster_graph(const Graph &G, const vector<int> &cluster, int &K) {
    int V=G.size();
    unordered_map<int,int> mp; K=0;
    for(int c:cluster) if(!mp.count(c)) mp[c]=K++;
    Graph cg(K);
    unordered_set<long long> seen;
    for(int u=0;u<V;u++){
        for(int v:G[u]){
            int cu=mp[cluster[u]];
            int cv=mp[cluster[v]];
            if(cu==cv) continue;
            long long key=((long long)min(cu,cv)<<32)|max(cu,cv);
            if(!seen.count(key)){
                cg[cu].push_back(cv);
                cg[cv].push_back(cu);
                seen.insert(key);
            }
        }
    }
    return cg;
}

vector<double> distribute_cluster_bc(const Graph &G,const vector<int> &cluster,const vector<double> &bc_c){
    int V=G.size();
    unordered_map<int,vector<int>> mem=group_members(cluster);
    vector<double> out(V,0.0);

    unordered_map<int,int> mp; int k=0;
    for(auto &p:mem) mp[p.first]=k++;

    for(auto &p:mem){
        int cid=p.first;
        int mcid=mp[cid];
        double total=0;
        for(int v:p.second) total += max(1,(int)G[v].size());
        if(total==0){
            double share=bc_c[mcid]/p.second.size();
            for(int v:p.second) out[v]+=share;
        } else {
            for(int v:p.second)
                out[v]+= bc_c[mcid]*((double)max(1,(int)G[v].size())/total);
        }
    }
    return out;
}

// -------------------- Improved Cluster Based BC --------------------
vector<double> cluster_based_bc_louvain(const Graph &G){
    int V=G.size();
    Louvain LV(G);
    vector<int> cluster = LV.run();

    unordered_map<int,vector<int>> members = group_members(cluster);
    vector<double> delta_local(V,0.0);

    // NEW: boundary nodes fix
    vector<int> boundary = find_boundary_nodes(G, cluster);

    for(auto &kv:members){
        vector<int> sources = kv.second;
        for(int b:boundary) sources.push_back(b);

        vector<double> part = brandes_from_sources_subset(G, sources);
        for(int i=0;i<V;i++) delta_local[i]+=part[i];
    }

    int K;
    Graph cg = build_cluster_graph(G,cluster,K);
    vector<double> bc_c = (K>0? brandes_full(cg): vector<double>(1,0));

    vector<double> delta_global = distribute_cluster_bc(G, cluster, bc_c);

    vector<double> BC(V);
    for(int i=0;i<V;i++)
        BC[i]=delta_local[i]+delta_global[i];

    return BC;
}

// -------------------- Helpers --------------------
vector<pair<double,int>> rank_nodes(const vector<double>&BC){
    vector<pair<double,int>> r;
    for(int i=0;i<BC.size();i++) r.push_back({BC[i],i});
    sort(r.begin(),r.end(),[](auto &a,auto &b){
        if(a.first!=b.first) return a.first>b.first;
        return a.second<b.second;
    });
    return r;
}

void print_vector(const vector<double>&A){
    cout<<fixed<<setprecision(6);
    for(int i=0;i<A.size();i++){
        cout<<A[i]<<(i+1<A.size()?' ':'\n');
    }
}
