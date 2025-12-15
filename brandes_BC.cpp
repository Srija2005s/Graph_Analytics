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
