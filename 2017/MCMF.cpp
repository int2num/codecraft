//
// Created by int2num on 3/24/17.
//

#include <cstdio>
#include "MCMF.h"

using namespace std;
MCMF::MCMF(char *topo[MAX_EDGE_NUM], int line_num) {
    int linei = 0;
    sscanf(topo[linei++], "%d%d%d", &n, &m, &consumers);
    linei++;
    danci=0;
    int init=0;
    vector<int>serv(3,-1);
    while(true)
    {
        if(sscanf(topo[linei++], "%d%d%d",&serv[0],&serv[1],&serv[2])<3)break;
        danci++;
        servinfo.push_back(serv);
        delta.push_back(serv[1]-init);
        init=serv[1];
        upper.push_back(serv[1]);
        fixedcharge.push_back(serv[2]);
    }
    cout<<"danci is "<<danci<<endl;
    lower.push_back(0);
    for(int i=0;i<danci-1;i++)
        lower.push_back(servinfo[i][1]);
    maxbandw=servinfo[danci-1][1];
    for(int i=0;i<n;i++)
    {
        int a=0,b=0;
        sscanf(topo[linei++], "%d%d",&a,&b);
        nodecharge.push_back(b);
    }
    linei++;
    g.resize(2*n + consumers + 3);
    s = n;
    function<void(int, int, int, int)> addedge = [&](int a, int b, int c, int d) {
        edges.emplace_back(a, b, c, d);
        g[a].emplace_back(edges.size() - 1);
        edges.emplace_back(b, a, 0, -d);
        g[b].emplace_back(edges.size() - 1);
    };
    for (int i = 0; i < m; i++) {
        int f, t, b, c;
        sscanf(topo[linei + i], "%d%d%d%d", &f, &t, &b, &c);
        addedge(f, t, b, c);
        addedge(t, f, b, c);
    }
    linei += m + 1;
    s = n, t = n + 1;
    jibie=1;
    //if(s>800)
       // jibie=1;
    n += consumers + 2;
    total_band = 0;
    for (int i = 0; i < consumers; i++) {
        int cid, adjid, bandwith;
        sscanf(topo[linei + i], "%d%d%d", &cid, &adjid, &bandwith);
        addedge(adjid, n - 1 - cid, bandwith, 0);
        total_band += bandwith;
        addedge(n - 1 - cid, t, bandwith, 0);
    }

    for(int j=0;j<s;j++) {
        int inidelta=0;
        vector<double>tvdelta;
        for (int i = 0; i < danci; i++) {
            tvdelta.push_back(fixedcharge[i] + nodecharge[j] - inidelta);
            inidelta = fixedcharge[i] + nodecharge[j];
        }
        vdelta.push_back(tvdelta);
    }
    cout<<"before init"<<endl;
    for(int j=0;j<s;j++)
    {
        vector<pair<int,int>>trx;
        vector<pair<int,int>> trxrange;
        vector<pair<int,int>> tsxrange;
        vector<pair<int,int>>tsx;
        int step=vdelta[j][0];
        trx.push_back(make_pair(step,0));
        trxrange.push_back(make_pair(0,1));
        trx.push_back(make_pair(0,step));
        trxrange.push_back(make_pair(1,upper[0]+1));
        tsx.push_back(make_pair(0,0));
        tsxrange.push_back(make_pair(lower[0],upper[0]));
        int slope=0;
        int bottom=0;
        for(int i=1;i<danci;i++)
        {   slope-=vdelta[j][i];
            bottom+=tsx[i-1].first*(tsxrange[i-1].second-tsxrange[i-1].first);
            trx.push_back(make_pair(slope,step));
            trxrange.push_back(make_pair(lower[i]+1,upper[i]+1));
            tsx.push_back(make_pair(-slope,bottom));
            tsxrange.push_back(make_pair(lower[i],upper[i]));
            step+=slope*(upper[i]-lower[i]);
        }
        rx.push_back(trx);
        rxrange.push_back(trxrange);
        sxrange.push_back(tsxrange);
        sx.push_back(tsx);
    }
    bakn=n,rn=bakn;
    for (int i = 0; i <s; i++) {
        addedge(s,i+n,maxbandw,0);
        for(int j=0;j<sx[i].size();j++)
        {
            addedge(n+i,i,sxrange[i][j].second-sxrange[i][j].first,sx[i][j].first);
            //cout<<sxrange[i][j].second-sxrange[i][j].first<<endl;
        }
    }
    offset=g[s][0];
    n=g.size();
}
int MCMF::getcostfromrx(double x,int index)
{
    int cost=0;
    for(int i=0;i<rx[index].size();i++)
    {
        if(x>=rxrange[index][i].first&&x<rxrange[index][i].second)
            cost+=rx[index][i].first*(x-rxrange[index][i].first)+rx[index][i].second;
        /*if(x>sxrange[index][i].first&&x<=sxrange[index][i].second)
        {   cout<<(x-sxrange[index][i].first)<<" "<<sx[index][i].first*(x-sxrange[index][i].first)<<" "<<sx[index][i].second<<endl;
            cost+=sx[index][i].first*(x-sxrange[index][i].first)+sx[index][i].second;
        }*/
    }

    return cost;

}
vector<pair<int,int>>MCMF::purelab(int dst)
{
    vector<double> dis(n, INT_MAX / 2);
    vector<int> vis(n, 0);
    int piS = 0, flow = 0;
    double cost = 0;
    dis = vector<double>(n, INT_MAX / 2);
    vector<int> inq(n, 0);
    dis[dst] = 0;
    deque<int> q;
    q.push_back(dst);
    inq[dst] = 1;
    while (!q.empty()) {
        int u = q.front();
        q.pop_front();
        inq[u]--;
        for (unsigned i = 0; i < g[u].size(); i++) {
            if (edges[g[u][i] ^ 1].bandwith == 0) continue;
            Edge &edge = edges[g[u][i]];
            if (dis[u] - edge.unit_price+eps< dis[edge.to]) {
                dis[edge.to] = dis[u] - edge.unit_price;
                if (inq[edge.to]) continue;
                if (!q.empty() && dis[edge.to]+eps< dis[q.front()]) q.push_front(edge.to);
                else q.push_back(edge.to);
                inq[edge.to]++;
            }
        }
    }
    vector<pair<int,int>>result;
    for(int i=0;i<s;i++)
    {
        int nodeid=i+rn;
        if(dis[nodeid]<0&&edges[g[s][i]^1].bandwith<=0)
        {
            result.push_back(make_pair(i,dis[nodeid]));
        }
    }
    sort(result.begin(),result.end(),quegreater());
    //for(int i=0;i<result.size();i++)
        //cout<<result[i].first<<":"<<result[i].second<<" ";
    //cout<<endl;
    return  result;
}
pair<double, int> MCMF::mcmf(int src, int dst, enum SPWAY cache = NORMAL) {
    vector<double> dis(n, 0);
    vector<int> vis(n, 0);
    int piS = 0, flow = 0;
    double cost = 0;
    path_result.clear();

    function<bool()> relable2 = [&]() -> bool {
        double d = INT_MAX / 2;
        for (int i = 0; i < n; i++) {
            if (vis[i]) {
                for (int j = 0; j < g[i].size(); j++) {
                    Edge &edges1 = edges[g[i][j]];
                    if (edges1.bandwith && !vis[edges1.to]) {
                        d = min(d * 1.0, dis[edges1.to] + edges1.unit_price - dis[edges1.from]);
                    }
                }
            }
        }
        if (d == INT_MAX / 2) return false;
        for (int i = 0; i < n; i++)
            if (vis[i]) {
                dis[i] += d;
            }
        piS = dis[src];
        if (cache == ROTATE_DELETE) {
            return dis[src] < INT_MAX / 2;
        }
//        cout<<"dis is "<<dis[src]<<endl;
        if (cache == ROTATE)
            return dis[src] < 0;
        return dis[src] < INT_MAX / 2;
    };
    function<bool()> relable = [&]() -> bool {
        dis = vector<double>(n, INT_MAX / 2);
        vector<int> inq(n, 0);
        dis[dst] = 0;
        deque<int> q;
        q.push_back(dst);
        inq[dst] = 1;
        while (!q.empty()) {
            int u = q.front();
            q.pop_front();
            inq[u]--;
            for (unsigned i = 0; i < g[u].size(); i++) {
                if (edges[g[u][i] ^ 1].bandwith == 0) continue;
                Edge &edge = edges[g[u][i]];
                if (dis[u] - edge.unit_price + eps < dis[edge.to]) {
                    dis[edge.to] = dis[u] - edge.unit_price;
                    if (inq[edge.to]) continue;
                    if (!q.empty() && dis[edge.to] + eps < dis[q.front()]) q.push_front(edge.to);
                    else q.push_back(edge.to);
                    inq[edge.to]++;
                }
            }
        }
        piS += dis[src];
        if (cache == ROTATE_DELETE) {
            return dis[src] < INT_MAX / 2;
        }
        //cout<<"dis is "<<dis[src]<<endl;
        if (cache == ROTATE)
            return dis[src] < 0;
        return dis[src] < INT_MAX / 2;
    };
    vector<int> curpath;
    function<int(int src, int u)> aug = [&](int augsrc, int u) -> int {
        if (augsrc == dst) {
            if (cache == ROUTE) {
                path_result.push_back(curpath);
                path_result.back().push_back(u);
            }
            flow += u;
            cost += piS * u;
            return u;
        }
        vis[augsrc] = 1;
        int l = u;
        for (unsigned i = 0; i < g[augsrc].size(); i++) {
            Edge &edge = edges[g[augsrc][i]];
            if (edge.bandwith > 0 && !vis[edge.to] && abs(dis[edge.from] - dis[edge.to] - edge.unit_price) < eps) {
                if (cache == ROUTE) curpath.push_back(g[augsrc][i]);
                int d = aug(edge.to, l < edge.bandwith ? l : edge.bandwith);
                if (cache == ROUTE) curpath.pop_back();
                edge.bandwith -= d, l -= d;
                if (cache != ROUTE)
                    edges[g[augsrc][i] ^ 1].bandwith += d;
                if (!l) return u;
            }
        }
        return u - l;
    };
    if(cache==ROTATE || cache==ROTATE_DELETE) relable();
    do {
        do {
            vis = vector<int>(n, 0);
        }while (aug(src, total_band));
    } while (relable2());
    return make_pair(cost, flow);
}
void MCMF::reset() {
//    reset the bandwith
    for (unsigned i = 0; i < edges.size(); i += 2) {
        edges[i].bandwith += edges[i ^ 1].bandwith;
        edges[i ^ 1].bandwith = 0;
    }
    vector<double> mj(g[s].size(), 0);
    setmj(mj);
}
void MCMF::setmj(const vector<double> &mj) {
    for (unsigned i = 0; i < g[s].size(); i++) {
        edges[g[s][i]].unit_price = mj[i];
        edges[g[s][i] ^ 1].unit_price = -mj[i];
    }
}
void MCMF::setmji(int index,double value) {
        edges[g[s][index]].unit_price=value;
        edges[g[s][index] ^ 1].unit_price=-value;
}
void MCMF::setcap(int index,int cap) {
        edges[g[s][index]].bandwith=cap;
        edges[g[s][index] ^ 1].bandwith=0;
}
void MCMF::setcoste(int index,double vv){
    edges[index*2].unit_price=vv;
    edges[index*2+1].unit_price=-vv;
}
void MCMF::setcape(int index,int cap) {
    edges[index].bandwith=cap;
}
int MCMF::rotate_delete(int a,vector<Edge>&back)
{
    vector<Edge>bak=edges;
    edges=back;
    int store=edges[g[s][a]^1].bandwith;
    edges[g[s][a]].bandwith=0;
    edges[g[s][a]^1].bandwith=0;
    edges.push_back(Edge(a+rn,n-1,store,0));
    g[a+rn].push_back(edges.size()-1);
    edges.push_back(Edge(n-1,a+rn,0,0));
    g[n-1].push_back(edges.size()-1);
    pair<double,int> pp=mcmf(s,n-1,ROTATE_DELETE);
    edges.pop_back(),edges.pop_back();
    g[a+rn].pop_back(),g[n-1].pop_back();
    if(pp.second>=store)
    {return pp.first;}
    else
    {edges=bak;return INFCOST;}
}
int MCMF::rotate(int a,vector<Edge>&back)
{
    vector<Edge>bbk=edges;
    edges.assign(back.begin(),back.end());
    edges.push_back(Edge(s,n-1,maxbandw,0));
    g[s].push_back(edges.size()-1);
    edges.push_back(Edge(n-1,s,0,0));
    g[n-1].push_back(edges.size()-1);
    pair<double,int> pp=mcmf(a+rn,n-1,ROTATE);
    edges.pop_back(),edges.pop_back();
    g[s].pop_back(),g[n-1].pop_back();
    if(pp.first<0) {
        //cout<<"add success"<<endl;
        edges[g[s][a] ^ 1].bandwith = pp.second;
        edges[g[s][a]].bandwith = maxbandw - pp.second;
        //cout<<a<<" "<<pp.first<<" "<<pp.second<<" "<<get_cost_allocation().first<<endl;
    }
    else
        edges=bbk;
    return pp.first;
}
pair<int,vector<Edge>> MCMF::push(int a,int band,vector<int>flowx,vector<Edge>bak)
{
    edges=bak;
    for(int i=0;i<s;i++)
    {   if(flowx[i]>0){
            int dan=getlevelofb(flowx[i]).first;
            setcape(g[s][i],upper[dan]-flowx[i]);
        }
    }
    edges[g[s][a]].bandwith=0;
    edges[g[s][a]^1].bandwith=0;
    edges.push_back(Edge(a+bakn,n-1,band,0));
    g[a+bakn].push_back(edges.size()-1);
    edges.push_back(Edge(n-1,a+bakn,0,0));
    g[n-1].push_back(edges.size()-1);
    int rp=mcmf(s,n-1).second;
    edges.pop_back(),edges.pop_back();
    g[a+bakn].pop_back(),g[n-1].pop_back();
    if(rp>=band)
    {
        edges[g[s][a]^1].bandwith=flowx[a]-band;
        return make_pair(get_cost_allocation().first,edges);
    }
    else {
        edges = bak;
        return make_pair(-1,edges);
    }

}
int MCMF::compress(string key) {
    int bestv=TaoMCMf(key);
    cout<<"copress init is "<<bestv<<endl;
    vector<int>flows = getflowx();
    bestflows=getallflow();
    bestflowx=getflowx();
    minvalue=bestv;
    vector<Edge> bak = edges;
    vector<Edge> bbk = edges;
    int fealure=0;
    for(int k=0;k<s;k++)
    {
        vector<pair<int, int>> over;
        vector<pair<int, int>> remain;
        for (int i = 0; i < s; i++) {
        int dan = getlevelofb(flows[i]).first;
        if (flows[i] > lower[dan]) {
            over.push_back(make_pair(i, flows[i] - lower[dan]));
            remain.push_back(make_pair(i, upper[dan] - flows[i]));
            }
        }
        sort(over.begin(), over.end(), quppcmp());
        pair<int, vector<Edge>> pp = push(over[0].first, over[0].second, flows, bak);
        if(pp.first<0)
            fealure++;
        if(fealure>10)
            break;
        edges = pp.second;
        if(pp.first<bestv&&pp.first>0)
        {
            bestv=pp.first;
            bestflows=getallflow();
            bestflowx=getflowx();
            minvalue=pp.first;
        }
        flows = getflowx();
        bak = edges;
    }
    edges=bbk;
    return bestv;
}
int MCMF::MAXCMCMF(double maxd)
{
    edges.push_back(Edge(t,n-1,maxd,0));
    g[t].push_back(edges.size() - 1);
    edges.push_back(Edge(n-1,t, 0,0));
    g[n-1].push_back(edges.size() - 1);
    pair<double,int>pp=mcmf(s,n-1);
    edges.pop_back(),edges.pop_back();
    g[t].pop_back(),g[n-1].pop_back();
    if(pp.second>=maxd)
        return pp.second;
    else
        return INFCOST;
}

vector<int> MCMF::getflowx(){
    vector<int> flowx;
    for(int i=0;i<g[s].size();i++)
        flowx.push_back(edges[g[s][i]^1].bandwith);
    return flowx;
}
vector<int> MCMF::getallflow(){
    vector<int> flows;
    for(int i=0;i<edges.size();i++)
        flows.push_back(edges[i^1].bandwith);
    return flows;
}
string MCMF::getstring()
{
    string ss;
    for(int i=0;i<g[s].size();i++)
        if(edges[g[s][i]^1].bandwith>0)
            ss.push_back('1');
        else
            ss.push_back('0');
    return ss;
}
pair<int,int> MCMF::getlevelofb(int bandwith)
{
    if(bandwith<0)
    {   cout<<"erro - flow"<<endl;
        return make_pair(-2,-1);
    }
    for(int i=0;i<servinfo.size();i++)
    {
        if(servinfo[i][1]>=bandwith)
            return make_pair(i,servinfo[i][2]);
    }
    cout<<"erro no range"<<endl;
    return make_pair(-1,-1);
}
pair<int, vector<int>> MCMF::get_cost_allocation() {
    int cost = 0;
    vector<int> xj(g[s].size(), 0);
    for (unsigned i = 0; i < edges.size(); i += 2) {
        if (edges[i ^ 1].bandwith == 0) continue;
        if(edges[i].from!=s)
            cost += edges[i ^ 1].bandwith * edges[i].unit_price;
    }
    for(int i=0;i<g[s].size();i++)
    {
        if(edges[g[s][i]^1].bandwith>0)
        {
            xj[i]=edges[g[s][i]^1].bandwith;
            if(bakn>0)
                cost+=getcostfromrx(xj[i],i);
            else {
                int dan=getlevelofb(xj[i]).first;
                cost += servinfo[dan][2]+nodecharge[i];
            }
        }
    }
    return make_pair(cost,xj);
}
string MCMF::print_result(vector<int> bestflow) {
    reset();
    for (unsigned i = 0; i < edges.size(); i++) {
        edges[i].bandwith = bestflow[i];
        edges[i ^ 1].bandwith = 0;
        i++;
    }
    mcmf(s, t, ROUTE);
    ostringstream out;
    out << path_result.size() << endl;
    int cost=0;
    for (unsigned i = 0; i < path_result.size(); i++) {
        out << endl;
        string seq ="";
        int node=edges[path_result[i][2]].from;
        int dist=0;
        for (unsigned j = 2; j < path_result[i].size() - 1; j++) {
            dist+=edges[path_result[i][j]].unit_price;
            if (j == path_result[i].size() - 2)
                out << seq << rn-edges[path_result[i][j]].from-1;
            else
                out << seq << edges[path_result[i][j]].from;
            seq = " ";
        }
        dist+=edges[path_result[i][1]].unit_price;
        cost+=dist*path_result[i].back();
        out << seq << path_result[i].back()<<seq<<getlevelofb(bestflowx[node]).first;
    }
    for(int i=0;i<bestflowx.size();i++)
        if(bestflowx[i]>0)
        {
            int vv=0;
            if(bakn>0)
                 vv=getcostfromrx(bestflowx[i],i);
            else
            {   int dan=getlevelofb(bestflowx[i]).first;
                vv= servinfo[dan][2]+nodecharge[i];
            }
            cost+=vv;
        }
    cout<<"really cost is "<<cost<<endl;
    return out.str();
}
int MCMF::TaoMCMf(string choose)
{
    reset();
    for(int i=0;i<choose.size();i++)
        if(choose[i]=='1')
            setcap(i,maxbandw);
        else
            setcap(i,0);
    pair<double,int> pp=mcmf(s,t);
    if(pp.second>=total_band)
        return get_cost_allocation().first;
    return INFCOST;
}
int MCMF::COMMCMf()
{
    mcmf(s,t);
    return get_cost_allocation().first;
}
void MCMF::GenMaxMinGraph()
{
    for(int i=0;i<s;i++)
    {
        maxnode.push_back(0);
        same.push_back(vector<pair<int,int>>());
    }
    for(int i=0;i<s;i++)
        mark.push_back(i);
    for(int i=0;i<s;i++) {
        int capacity=0;
        maxflow.push_back(0);
        for (int j = 0; j <g[i].size();j++)
            capacity+=edges[g[i][j]].bandwith;
        capacities.push_back(capacity);
    }
    reset();
    for(int i=0;i<s;i++)
        setmji(i,rx[i][0].first*1000);
    minvalue=COMMCMf();
    bestkey=getstring();
    cout<<bestkey<<endl;
    /*for(int i=0;i<bestkey.size();i++)
        if(bestkey[i]=='1')
            maxnode[i]++;*/
    cout<<"init value "<<minvalue<<endl;
    bestX=getallflow();
}
/*int MCMF::op_delete(string key,int vv)
{
    if(vv<0)
        return 1;
    double value=vv;
    reset();
    TaoMCMf(key);
    vector<Edge>bak=edges;
    reset();
    for(int i=0;i<key.size();i++){
        if (key[i] == '1') {
            key[i] = '0';
            if (Visited.find((key)) == Visited.end()) {
                double delta = rotate_delete(i, bak);
                double newv = get_cost_allocation().first;
                string kk = getstring();
                //candv[i]++;
                if (delta < INFCOST && newv < value + value / decrease) {
                    if (Visited.find(kk) == Visited.end()) {
                        bottomque.push(make_pair(kk, newv));
                    }
                    if (newv < minvalue) {
                        minvalue = newv;
                        bestkey = kk;
                        bestX = getallflow();
                        bestflowx = getflowx();
                        bestflows = getallflow();
                        cout << "change in delete search: " << minvalue << endl;
                        if(jibie > 0)
                        {cout<<"i is "<<i<<endl;return 1;}
                    }
                }
                Visited.insert(key);
                Visited.insert(kk);
            }
            key[i] = '1';
            gettimeofday(&endtime, NULL);
            if (endtime.tv_sec - starttime.tv_sec>500)
                return -1;
        }
    }
    return 1;
}*/
int MCMF::op_delete(string key,int vv)
{
    if(vv<0)
        return 1;
    double value=vv;
    reset();
    TaoMCMf(key);
    vector<Edge>bak=edges;
    reset();
    int i=delstore;
    while(true) {
        if (key[i] == '1') {
            key[i] = '0';
            if (Visited.find((key)) == Visited.end()) {
                double delta = rotate_delete(i, bak);
                double newv = get_cost_allocation().first;
                string kk = getstring();
                if (delta < INFCOST && newv < value + value / decrease) {
                    if (Visited.find(kk) == Visited.end()) {
                        bottomque.push(make_pair(kk, newv));
                    }
                    if (newv < minvalue) {
                        minvalue = newv;
                        bestkey = kk;
                        bestX = getallflow();
                        bestflowx = getflowx();
                        bestflows = getallflow();
                        if (jibie>0)
                        {
                            cout<<"-x "<<i<<" ";
                            delstore=(++i)%key.size();
                            return 1;
                        }
                    }
                }
                Visited.insert(key);
                Visited.insert(kk);
                if (jibie > 0&&state==0)
                {
                    cout<<"-x "<<i<<" ";
                    delstore=(++i)%key.size();
                    return 1;
                }
            }
            key[i] = '1';
            gettimeofday(&endtime, NULL);
            if (endtime.tv_sec - starttime.tv_sec >200)
                return -1;

        }
        i++;
        i=i%(key.size());
        if(i==delstore)
        {
            delstore=i;
            break;
        }
    }
    cout<<"no delete"<<endl;
    return 1;
}
/*int MCMF::op_add(string key,int value)
{
    if(value>0)
    {
        reset();
        TaoMCMf(key);
        vector<Edge>bak=edges;
        cout<<"candy size is "<<candy.size()<<endl;
        for(int j=0;j<candy.size();j++) {
            if(j>20)break;
            int i=candy[j];
            if (key[i] == '0' && maxnode[i] >0) {
                key[i] = '1';
                if (Visited.find((key)) == Visited.end()) {
                    int v = rotate(i, bak);
                    if (v < 0) {
                        int newv = get_cost_allocation().first;
                        string kk = getstring();
                        if (Visited.find(kk) == Visited.end()) {
                            bottomque.push(make_pair(kk, newv));
                        }
                        if (newv < minvalue) {
                            minvalue = newv;
                            bestkey = kk;
                            bestflows = getallflow();
                            bestflowx = getflowx();
                            cout << "change in add search: " << minvalue << endl;
                            if (jibie > 0) {
                                return 1;
                            }
                        }
                        Visited.insert(key);
                        Visited.insert(kk);
                    }
                }
                key[i] = '0';
                gettimeofday(&endtime, NULL);
                if (endtime.tv_sec - starttime.tv_sec >87)
                    return -1;
            }
            i++;
        }
    }
}*/
int MCMF::op_add(string key,int value)
{
    if(value>0)
    {
        reset();
        TaoMCMf(key);
        vector<Edge>bak=edges;
        int i=addstore;
        int addcount=0;
        while(true) {
            if (key[i] == '0' && maxnode[i] >0) {
                key[i] = '1';
                if (Visited.find((key)) == Visited.end()) {
                    int v = rotate(i, bak);
                    if (v < 0) {
                        int newv = get_cost_allocation().first;
                        string kk = getstring();
                        if (Visited.find(kk) == Visited.end()) {
                            bottomque.push(make_pair(kk, newv));
                        }
                        if (newv < minvalue) {
                            minvalue = newv;
                            bestkey = kk;
                            bestflows = getallflow();
                            bestflowx = getflowx();
                            cout << "change in add search: " << minvalue << endl;
                            if (jibie > 0) {
                                addstore = (++i) % key.size();
                                return 1;
                            }
                        }
                        Visited.insert(key);
                        Visited.insert(kk);
                    }
                }
                key[mark[i]] = '0';
                addcount++;
                if(addcount>5&&state==0)
                    break;
                gettimeofday(&endtime, NULL);
                if (endtime.tv_sec - starttime.tv_sec >87)
                    return -1;
            }
            i++;
            i=i%(key.size());
            if(i==addstore)
            {
                addstore=i;
                break;
            }
        }
        return 1;
    }
}
vector<int> MCMF::findcandy(int s,vector<Edge>bak){
    vector<Edge>bbk=edges;
    edges=bak;
    vector<pair<int,int>> vpp=purelab(s);
    vector<int>res;
    int k=0;
    for(int i=0;i<vpp.size();i++)
            res.push_back(vpp[i].first);
    edges=bbk;
    return res;
};
/*int MCMF::op_add(string key,int value)
{
    if(value>0)
    {
        reset();
        TaoMCMf(key);
        vector<Edge>bak=edges;
        int j=0;
        set<int>haset;
        for(int j=0;j<30;j++) {
            int ran=-1;
            while(true){
                ran=rand()%candynode.size();
                if(haset.find(ran)==haset.end()){
                    haset.insert(ran);break;
                }
            }
            int i=candynode[ran];
            if (key[i] == '0') {
                key[i] = '1';
                if (Visited.find((key)) == Visited.end()) {
                    int v = rotate(i, bak);
                    if (v < 0) {
                        int newv = get_cost_allocation().first;
                        string kk = getstring();
                        if (Visited.find(kk) == Visited.end()) {
                            bottomque.push(make_pair(kk, newv));
                        }
                        if (newv < minvalue) {
                            minvalue = newv;
                            bestkey = kk;
                            bestflows = getallflow();
                            bestflowx = getflowx();
                            //cout << "change in add search: " << minvalue << endl;
                            if (jibie > 0)
                            {
                                cout<<"+"<<i<<" m:"<<maxnode[i]<<endl;
                                return 1;
                            }

                        }
                        Visited.insert(key);
                        Visited.insert(kk);
                    }
                }
                key[i] = '0';
                gettimeofday(&endtime, NULL);
                if (endtime.tv_sec - starttime.tv_sec >200)
                    return -1;
            }
        }
        cout<<"n +"<<endl;
        return 1;
    }
}*/
void updatemax(vector<double >&max,vector<double >&now)
{
    for(int i=0;i<max.size();i++)
        if(max[i]<now[i])
            max[i]=now[i];
}
void updatefor(vector<double >&max,vector<double >&now)
{
    for(int i=0;i<max.size();i++)
        if(now[i]>0)
            max[i]=now[i];
}
void MCMF::nearsearch()
{
    int j=0;
    for(int i=0;i<s;i++)
        if(maxnode[i]>0)
        {   j++;
            candv.push_back(0);
            candynode.push_back(i);
        } else
            candv.push_back(INFCOST);
    cout<<endl;
    cout<<"candidate set num: "<<j<<endl;
    bottomque.push(make_pair(bestkey,minvalue));
    cout<<"begin value: "<<minvalue<<endl;
    int count=0;
    decrease=4;
    addstore=0;
    delstore=0;
    state=0;
    while(!bottomque.empty())
    {
        int j=0;

        sort(candies.begin(),candies.end(),quppcmp());
        pair<string,int>pp=bottomque.top();
        bottomque.pop();
        cout<<"pop "<<pp.second<<endl;
        decrease*=decrease*1.2;
        if(pp.second==minvalue)
        {
            changeque.push(pp);
            decrease=4;count=0;
        }
        count++;
        if(op_delete(pp.first,pp.second)<0||op_add(pp.first,pp.second)<0)break;
    }
    while(!bottomque.empty())bottomque.pop();
    cout<<"final min value is: "<<minvalue<<endl;
}
int MCMF::getseg(double band,int index)
{
    for(int i=0;i<rxrange[index].size();i++)
        if(band>=rxrange[index][i].first&&band<=rxrange[index][i].second)
            return i;
    return -1;
}
pair<vector<int>,vector<int>>MCMF::solve_dcup(vector<double>epsilon,vector<int>ya)
{
    cout<<"solve dcup"<<endl;
    vector<int>initX(s,0);
    //ya.assign(initX.begin(),initX.end());
    vector<int> pre(s,-1);
    for(int k=0;k<10000;k++) {
        reset();
        for (int i=0;i<s;i++)
            if (ya[i]==0)
                setmji(i, (double)nodecharge[i]/epsilon[i]);
        double value= COMMCMf();
        initX=getflowx();
        for(int i=0;i<initX.size();i++)
            if(maxflow[i]<initX[i])
                maxflow[i]=initX[i];
        cout<<"nw it:"<<value<<endl;
        if(value<minvalue) {
            bestkey=getstring();
            minvalue = value;
            bestflows=getallflow();
            bestflowx=getflowx();
        }
        pre.assign(ya.begin(),ya.end());
        for (int i=0;i<s;i++)
        {
            if (initX[i]>epsilon[i])
            {ya[i]=1;}
            else
                ya[i]=0;
            if(initX[i]>0)
                maxnode[i]++;
        }
        int flag=0;
        for(int i=0;i<ya.size();i++)
            if(ya[i]!=pre[i])
            {flag=1;break;}
        if(flag==0)
        {
            break;
        }
    }
    return make_pair(initX,ya);
}
void MCMF::solve_adcup(vector<int>ya)
{
    double alpha=0.5;
    cout<<"solve ssss"<<endl;
    niceX=ya;
    vector<double> epsilon;
    for(int i=0;i<s;i++)
        if(ya[i]>0)
            ya[i]=1;
    epsilon.clear();
    for(int i=0;i<capacities.size();i++)
        epsilon.push_back(capacities[i]);
    int flag=0;
    vector<int> preya;
    for(int i=0;i<20000;i++) {
        flag=0;
        for(int j=0;j<capacities.size();j++)
            if(niceX[j]<epsilon[j]&&epsilon[j]>1)
            {epsilon[j]=epsilon[j]*alpha;
                flag=1;
            }
        if(flag==0)
        {cout<<"early break"<<endl;break;}
        preya=niceX;
        pair<vector<int>,vector<int>>pp=solve_dcup(epsilon,ya);
        ya=pp.second;
        niceX=pp.first;
    }
    solvekey=getstring();
}
pair<vector<Edge>,int> MCMF::changeGraph()
{
    reset();
    vector<Edge> bakedges=edges;
    for(int i=0;i<s;i++)
    {
        edges[g[rn+i][1]].bandwith=maxbandw;
        edges[g[rn+i][1]].unit_price=0;
        edges[g[rn+i][1]^1].unit_price=0;
        for(int j=2;j<g[rn+i].size();j++)
            edges[g[rn+i][j]].bandwith=0,edges[g[rn+i][j]^1].bandwith=0;
    }
    int bbn=bakn;bakn=0;
    return make_pair(bakedges,bbn);
}
//useless
/*pair<vector<int>,vector<int>>MCMF::solve_dcup(vector<int>initx)
{
    cout<<"solve dcup"<<endl;
    vector<int>y;
    for(int i=0;i<s;i++)
        y.push_back(getseg(initx[i],i));
    vector<int> pre(s,-1);
    for(int k=0;k<100;k++) {
        reset();
        for(int i=0;i<s;i++)
        {
            setmji(i,rx[i][y[i]].first);
            cout<<rx[i][y[i]].first<<" ";
        }
        cout<<endl;
        int newv=COMMCMf();
        cout<<"newv is"<<newv<<endl;
        vector<int>flowx=getflowx();
        for(int i=0;i<s;i++)
            y[i]=getseg(flowx[i],i);
        for(int i=0;i<s;i++)
            cout<<y[i]<<" ";
        cout<<endl;
    }
    return make_pair(initx,pre);
}*/
pair<int,string> MCMF::localsearch()
{
    int count=0;
    vector<int>bkey;
    for(int i=0;i<localkey.size();i++)
        if(localkey[i]=='1')
            count++,bkey.push_back(i);
    int choose=rand()%count;
    choose=bkey[choose];
    string key=localkey;
    key[choose]='0';
    int init=INFCOST;
    string bbk=localkey;
    for(int i=0;i<localkey.size();i++)
        if(localkey[i]=='0')
        {
            key[i]='1';
            int gg=TaoMCMf(key);
            if(gg<init)
                init=gg,bbk=getstring();
            key[i]='0';
        }
    return make_pair(init,bbk);
};
void MCMF::changebak(pair<vector<Edge>,int> pp)
{
    edges=pp.first,bakn=pp.second;
}

int MCMF::compresearch()
{
    int k=0;
    map<int,int> mv;
    totalmaxnode=0;
    for(int i=0;i<maxnode.size();i++)
        if(maxnode[i]>0)
        {
            nodepriority.insert(make_pair(i,maxnode[i])),k++;
            candynode.push_back(i);
            candies.push_back(make_pair(i,maxnode[i]));
        }
    cout<<"candidate set is "<<k<<endl;
    vector<Edge>allbak=edges;
    TaoMCMf(bestkey);
    vector<Edge>bak=edges;
    int i=-1;
    int j=-1;
    int cont=0;
    while(true)
    {
        cont++;
        i=(++i)%bestkey.size();
        gettimeofday(&endtime, NULL);
        if (endtime.tv_sec - starttime.tv_sec >87)
        { edges=allbak;return -1;}
        if(bestkey[i]=='1')
        {
            int delsucess = rotate_delete(i, bak);
            if (delsucess < INFCOST);
            {
                mv[i]=cont;
                double newv = get_cost_allocation().first;
                if (newv < minvalue) {
                    cout << "change in delete " << minvalue << endl;
                    minvalue = newv;
                    bestkey = getstring();
                    bak = edges;
                }
                int count=0;
                //cout<<"b a"<<endl;
                while(true)
                {
                    j=(++j)%candynode.size();
                    gettimeofday(&endtime, NULL);
                    if (endtime.tv_sec - starttime.tv_sec >87)
                    { edges=allbak;return -1;}
                    count++;
                    if(count>30)
                        break;
                    int addnode=candynode[j];
                    set<int>haset;
                    while(true){
                        addnode=rand()%candynode.size();
                        if(haset.find(addnode)==haset.end()){
                            haset.insert(addnode);break;
                        }
                    }
                    addnode=candynode[addnode];
                    if(bestkey[addnode]=='0') {
                        int delsucess=rotate(addnode,edges);
                        nodepriority[addnode]-= 1;
                        int newv=get_cost_allocation().first;
                        //if(delsucess<0)
                            //cout<<newv<<endl;
                        if (delsucess<0&&newv < minvalue) {
                            cout << "change in add " << minvalue << endl;
                            minvalue = newv;
                            bak = edges;
                            bestkey=getstring();
                            mv[addnode]=cont;
                            break;
                        }
                    }
                }
                //cout<<"af add"<<endl;
            }
        }
    }
    /*while(true)
    {
        cont++;
        i=(++i)%bestkey.size();
        gettimeofday(&endtime, NULL);
        if (endtime.tv_sec - starttime.tv_sec >87)
        { edges=allbak;return -1;}
        if(bestkey[i]=='1')
        {
            //cout<<"before dd"<<endl;
            int delsucess = rotate_delete(i, bak);
            if (delsucess < INFCOST);
            {
                double newv = get_cost_allocation().first;
                if (newv < minvalue) {
                    cout << "change in delete " << minvalue << endl;
                    minvalue = newv;
                    bestkey = getstring();
                    bak = edges;
                }
                int count=0;
                j=-1;
                while(true)
                {
                    j=(++j)%candynode.size();
                    gettimeofday(&endtime, NULL);
                    if (endtime.tv_sec - starttime.tv_sec >87)
                    { edges=allbak;return -1;}
                    count++;
                    if(count>30)
                        break;
                    int addnode=candynode[j];
                    /*set<int>haset;
                    while(true){
                        addnode=rand()%candynode.size();
                        if(haset.find(addnode)==haset.end()&&cont-mv[addnode]>5){
                            haset.insert(addnode);break;
                        }
                    }
                    addnode=candynode[addnode];
                    if(bestkey[addnode]=='0') {
                        int rosucess=rotate(addnode,edges);
                        nodepriority[addnode]-= 1;
                        newv = get_cost_allocation().first;
                        if (rosucess<0&&newv < minvalue) {
                            cout << "change in add " << minvalue << endl;
                            minvalue = newv;
                            bak = edges;
                            bestkey=getstring();
                            //cout<<bestkey<<endl;
                            break;
                        }
                    }
                }
            }
        }
    }*/
    return 0;
}
string MCMF::NewSolve()
{
    gettimeofday(&starttime,NULL);
    srand(0);
    //string bbk="000000000000000000000000000000010001000000000000000000000000110000000000000000000000000000100000000000000000000001000000000000000000000000000000000000001000000000100010000000000000000000000000000100000000000000000000100000000001000000000000001100000000000000001000000000000000100000000010000000000000000000100000000001100000000000000001010000000010000000011010000000000000000010100000000000000000101100000000000110000000000000000000000000000010000000000000010000000000000000000000000000010000000100000001010001001000000000000000000000000000000000000000100000000000000000100000000100000000000000010000";
    //cout<<"gold check for me "<<TaoMCMf(bbk)<<endl;
    GenMaxMinGraph();
    pair<vector<Edge>,int>pp=changeGraph();
    vector<int>y(s,0);
    solve_adcup(y);
    changebak(pp);
    if(TaoMCMf(bestkey)>minvalue)
    {
        cout<<"first kind "<<endl;
        changeGraph();
    }
    compresearch();
    cout<<"after compress"<<compress(bestkey)<<endl;
    return print_result(bestflows);
}


