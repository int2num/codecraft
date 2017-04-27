#ifndef CDN_SKW_H_H
#define CDN_SKW_H_H
#endif //CDN_SKW_H_H
/*autor peng*/
#include <cstdio>
#include <cstring>
#include <iostream>
#include "string.h"
#include <sstream>
#include <float.h>
#define N 3000
#define M 10000
#define inf  1000000
#define eps 1e-11
struct zkw {
    int C;
    double D;
    int s, t,n,m;
    int fst[N],nxt[M],vv[M],cap[M],flow[M],e;
    int backcap[M];
    double cost[M],d[N];
    bool vis[N];
    int F,offset;
    void init() {
        memset(fst, -1, sizeof fst);
        e = 0;
        C = D = 0;
    }
    void add(int u, int v, int f, double c) {
        vv[e] = v, nxt[e] = fst[u], cost[e] = c, cap[e] = f,backcap[e]=f,flow[e] = 0, fst[u] = e++;
        vv[e] = u, nxt[e] = fst[v], cost[e] = -c, cap[e] =0,backcap[e]=0,flow[e] = 0, fst[v] = e++;
    }
    double aug(int u, double f) {
        if(u == t) {
            C += D * f;
            F += f;
            return f;
        }
        vis[u] = 1;
        double tmp = f;
        for(int i = fst[u]; ~i; i = nxt[i]) {
            int v = vv[i];
            if(cap[i] > flow[i] && fabs(cost[i])<eps && !vis[v]) {
                double d = aug(v, tmp < cap[i] - flow[i]? tmp: cap[i] - flow[i]);
                flow[i] += d;
                flow[i^1] -= d;
                tmp -= d;
                if(!tmp) return f;
            }
        }
        return f - tmp;
    }
    bool modLabel() {
        for(int i=0;i<n;i++)
            d[i]=inf;
        d[t] = 0;
        deque<int> q; q.push_back(t);
        while(!q.empty()) {
            double dt;
            int u = q.front(); q.pop_front();
            //cout<<u<<endl;
            for(int i = fst[u]; ~i; i = nxt[i]) {

                int v = vv[i];
                if(cap[i^1] > flow[i^1] &&(dt=d[u]-cost[i])<d[v]) {
                    d[v] = dt;
                    double tmp = inf;
                    if(q.size() > 0) tmp = d[q.front()];
                    if((dt-tmp)<eps) q.push_front(v);
                    else q.push_back(v);
                }
            }
        }
        for(int i = 0; i <= t; ++i) {
            for(int j = fst[i]; ~j; j = nxt[j]) {
                int v = vv[j];
                cost[j] += d[v] - d[i];
            }
        }
        D += d[s];
        return d[s] < inf;
    }
    int MCMF(int s, int t) {
        this -> s = s, this -> t = t;
        F = 0;
        memset(flow,0,sizeof flow);
        while(modLabel()) {
            do memset(vis, 0, sizeof vis);
            while (fabs(aug(s, inf))>eps);
        }
        for(int i=0;i<m;i++)
            if(cap[i]<0)
                cout<<cap[i]<<endl;
        //cout<<F<<endl;
        return F;
    }
    //allay functions
    void getflows(int num,vector<int>&result)
    {
        int cc=0;
        while(cc<num)
        {
            result[cc]=flow[offset+2*cc];
            cc++;
        }
    }
    string getstring(int num)
    {

        stringstream st;
        int cc=0;
        while(cc<num)
        {
            if(flow[offset+2*cc]>0)
                st<<"1";
            else
                st<<"0";
            cc++;
        }
        return st.str();
    }
    void setoffset(int off,int node,int medge)
    {
        n=node;
        m=medge;
        offset=off;
    }
    void setcost(int index,double value){
        cost[index*2+offset]=value;
        cost[index*2+offset+1]=0;
    }
    void setcapacity(int index,int value) {
        cap[2 * index + offset] = value;
        cap[2 * index + offset + 1] = 0;
    }
    void clearcost()
    {
        memset(cost+offset,0,sizeof(int)*(M-offset));
    }
    void loadcapacity(){
        memcpy(cap,backcap,sizeof cap);
    }
};