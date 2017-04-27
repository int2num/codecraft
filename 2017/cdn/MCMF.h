//
// Created by int2num on 3/24/2017.
//

#ifndef CDN_MCMF_H
#define CDN_MCMF_H

#include <lib_io.h>
#include <bits/stdc++.h>
#include<sys/time.h>
#define MAXNODE 5000
#define INFCAPACITY 100000000
#define INFCOST INT_MAX / 2
#define eps 1e-8
#define MAXTYPE 1
#define FORTYPE 2
#define SERTYPE 3
using namespace std;
enum SPWAY {NORMAL,ROUTE,ROTATE,ROTATE_DELETE,PUSH};
struct Edge {
    int from, to, bandwith;
    double unit_price;

    Edge(int f, int t, int b, double up) : from(f), to(t), bandwith(b), unit_price(up) {};
};
class paircmp
{
public:
    bool operator()(pair<string,int>a,pair<string,int>b)
    {if(fabs(a.second)>fabs(b.second))
            return true;
        return false;};
};
class quppcmp
{
public:
bool operator()(pair<int,int>a,pair<int ,int>b)
{
    if(a.second<b.second)
        return true;
    return false;};
};
class quegreater
{
    public:
    bool operator()(pair<int,int>a,pair<int ,int>b)
    {
        if(a.second<b.second)
            return true;
        return false;};
};
class ppcmp{
public:
    bool operator()(pair<int,double>a,pair<int,double>b){
        if(a.second>b.second)
            return true;
        return false;
    }
};
class greatcmp{
public:
    bool operator()(pair<int,int> a,pair<int,int> b)
    {return a.second>b.second;}
};

class MCMF {
public:
    int n, m, s, t;
    int consumers;
    int cdn_cost;
    int total_band;
    vector<vector<int>> path_result;
    vector<Edge> edges;
    vector<vector<int>> g;
    vector<int>maxnode;
    vector<int>maxflow;
    priority_queue<pair<int,int>,vector<pair<int,int>>,quppcmp> delset,addset;
    vector<pair<int,double>>candidateset;
    string solvekey,servkey,maxkey,localkey;
    vector<int>mark;
    vector<int>capacities;
    vector<int>niceX;
    vector<int>bestX;
    unordered_set<string> Visited;
    int minvalue,localvalue;
    string bestkey;
    vector<int>bestflows;
    vector<int>bestflowx;
    unordered_set<string> blf;
    vector<vector<pair<int,int>>> same;
    double decrease;
    struct timeval starttime, endtime;
    priority_queue<pair<int, string>, vector<pair<int, string>>, greater<pair<int, string>>> q;
    pair<double, int> mcmf(int src, int dst, enum SPWAY cache) ;
    void reset();
    void setmj(const vector<double> &mj);
    void setcap(int index,int cap);
    pair<int, vector<int>> get_cost_allocation();
    string print_result(vector<int>);
    vector<double> mask2mjcdn(const string &);
    void setcoste(int index,double vv);
    void setcape(int index,int cap);
    int MAXCMCMF(double maxd);
    vector<vector<int>> servinfo;
    int danci;
    int maxbandw;
    vector<double> delta;
    vector<double> lower;
    vector<double>upper;
    vector<double>nodecharge;
    int offset;
    int bakn,rn,jibie;
    int delstore,addstore;
    vector<pair<int,int>>candies;
    vector<int>candv;
    vector<int>candy;
    int state;
    vector<int>candynode;
    vector<double>deltaprice;
    vector<double> fixedcharge;
    void solve_sub(vector<int>);
    int compress(string key);
    int compress(vector<int>);
    int push(int a,int band,vector<Edge>&back);
    pair<int,vector<Edge>> push(int a,int band,vector<int>flowx,vector<Edge>bak);
    int sa();
    pair<int,string> localsearch();
    //add by zhang
    priority_queue <pair<string,int>,vector<pair<string,int>>,paircmp> bottomque;
    priority_queue <pair<string,int>,vector<pair<string,int>>,paircmp> changeque;
    int TaoMCMf(string choose);
    int COMMCMf();
    void GenMaxMinGraph();
    void setmji(int index,double value);
    vector<int> getflowx();
    vector<int> getallflow();
    string getstring();
    void solve_adcup(vector<int>ya);
    int op_delete(string key,int vv);
    int op_add(string key,int value);
    void nearsearch();
    pair<string,pair<int,vector<int>>> arrange(pair<string,pair<int,vector<int>>>data,int type,int dd);
    void Lagrange(int dd,int ss);
    int rotate(int a,vector<Edge>&back,int level=0);
    int rotate_delete(int a,vector<Edge>&back);
    int replace(int a,int b,vector<Edge>&bak);
    vector<vector<int>> getroutes(vector<int> bestflow);
    void solve_adcup(vector<int>ya,vector<int>delcap);
    pair<int,int> getlevelofb(int bandwith);
    string xjbh(string key);
    vector<vector<double>>vdelta;
    vector<vector<pair<int,int>>>rx;
    vector<vector<pair<int,int>>>rxrange;
    vector<vector<pair<int,int>>>sx;
    vector<vector<pair<int,int>>>sxrange;
    int getcostfromrx(double x,int index);
    int getseg(double band,int i);
    pair<vector<int>,vector<int>>solve_dcup(vector<double>epsilon,vector<int>ya);
    pair<vector<Edge>,int> changeGraph();
    void changebak(pair<vector<Edge>,int> pp);
    vector<int> findcandy(int s,vector<Edge>bak);
    vector<pair<int,int>>purelab(int dst);
    int compresearch();
    vector<Edge>classic;
    int testdel;
    map<int,int>nodepriority;
    int totalmaxnode;
    vector<pair<double,double>>maxrange;
    int choosenode();
    int zuirisearch();
    map<int,int>del,add;
    int distof(int a,int b,vector<Edge>bak);
    vector<Edge>dbak;
    vector<int>addevalue,delevalue;
    int getdelnode();
    int getaddnode();
    int updatedel(set<int>&index,int v);
    int updateadd(set<int>&index,int v);
    int gcount;
    map<int,int>ifadded;
    map<int,int> ifdeleted;
    int reversearch();
    vector<int> super_rotate(int a,vector<Edge>&back);
    pair<int,int> getcostpair(vector<Edge>&edges);
    void changegraphv();
public:
    MCMF(char *topo[MAX_EDGE_NUM], int line_num);
    string NewSolve();

};


#endif //CDN_MCMF_H
