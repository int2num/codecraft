//
// Created by int2num on 17-3-29.
//
#include "MCMF.h"
#include <set>
#ifndef CDN_TABU_H
#define CDN_TABU_H

#endif //CDN_TABU_H
class fixque{
public:
    int flag;
    vector<int>cap;
    fixque(int size):cap(size,-1),flag(0){};
    void push(int value){
        flag=flag%cap.size();
        cap[flag]=value;
    }
    bool find(int value)
    {
        for(int i=0;i<cap.size();i++)
            if(cap[i]=value)
                return true;
        return false;
    }

};
class Tabu {
public:
    vector<int> z, t, h;
    vector<int> forbiden;
    string tkey;
    vector<Edge> tedges;
    vector<pair<int, double >> candidates;
    fixque aspireset;
    MCMF &mm;
    int k, kesaia, kesaib, kesai_in, kesai_out, bestvalue, candmark;
    int nowv;
    set<int> hasdel, hasadd;
    map<int, priority_queue<pair<int, int>, vector<pair<int, int>>, quegreater>> adjnode;
    map<int, priority_queue<pair<int, int>, vector<pair<int, int>>, quegreater>> dejnode;
    vector<pair<int, double>> eva;

    Tabu(MCMF &mf, int size) : mm(mf), z(mf.s, INFCOST), t(mf.s, -100000), h(mf.s, 0), aspireset(size),
                               forbiden(mf.s, -10000), eva(mm.edges.size() / 2, make_pair(-1, 0)) {
        kesaia = 5, kesaib = 10, candmark = 0;
        candidates = mm.candidateset;
        tkey = mf.bestkey;
        mf.TaoMCMf(tkey);
        tedges = mf.edges;
        k = 0, nowv = mm.minvalue;
    };

    void genkesai() {
        kesai_in = rand() % (kesaib - kesaia);
        kesai_out = rand() % (kesaib - kesaia);
        kesai_in += kesaia;
        kesai_out += kesaia;
    };

    int add_delete(int opnum, set<int> &ford, set<int> &fora, int delta = 0, int thred = 1000) {
        mm.reset();
        int localmv = mm.TaoMCMf(tkey);
        vector<Edge> bak(mm.edges);
        vector<int> bpos(opnum, -1), dpos(opnum, -1);
        string localkey = tkey;
        mm.reset();
        for (int j = 0; j < opnum; j++) {
            for (int i = 0; i < candidates.size(); i++) {
                int pos = candidates[i].first;
                if (candidates[pos].second < thred && tkey[pos] == '0' && fora.find(pos) == fora.end() &&
                    (k - t[pos] > 4)) {
                    int can = mm.rotate(pos, bak);
                    if (mm.get_cost_allocation().first < localmv) {
                        localmv = mm.get_cost_allocation().first;
                        localkey = mm.getstring(), tkey = localkey, bpos[j] = pos;
                        hasadd.insert(pos);
                        ford.insert(pos);
                        break;
                    }
                }
            }
        }
        tkey = localkey;
        mm.TaoMCMf(tkey);
        bak = mm.edges;
        for (int j = 0; j < opnum; j++) {
            for (int i = 0; i < tkey.size(); i++)
                if (tkey[i] == '1' && ford.find(i) == ford.end() && (k - t[i] > 4)) {
                    string kk = tkey;
                    kk[i] = '0';
                    int del = mm.rotate_delete(i, bak);
                    if (del < INFCOST && mm.get_cost_allocation().first < localmv + delta) {
                        localmv = mm.get_cost_allocation().first;
                        localkey = mm.getstring();
                        tkey = localkey;
                        dpos[j] = i;
                        fora.insert(i);
                        hasdel.insert(i);
                        break;
                    }
                }
        }

        for (int i = 0; i < opnum; i++) {
            if (dpos[i] >= 0) { t[bpos[i]] = k; }
            if (bpos[i] >= 0) { t[dpos[i]] = k; }
        }
        if (localmv < mm.minvalue)
            mm.minvalue = localmv, mm.bestkey = tkey;
        k++;
        return localmv;
    }

    int localsearch() {
        int k = 0;
        int prev = mm.minvalue;
        set<int> empty;
        for (int i = 0; i < 100; i++) {
            int v = add_delete(1, empty, empty);
            cout << v << endl;
            if (prev == v)
                k++;
            if (k >= 5)break;
            prev = v;
        }

    }

    int replace(string key, int s, int t) {
        key[s] = '0';
        if (t >= 0)
            key[t] = '1';
        return mm.TaoMCMf(key);
    }

    pair<int, int> lkh(string &key, int value, int init, int kk = 1) {
        if (kk == 1)
            forbiden[init] = k;
        priority_queue<pair<int, int>, vector<pair<int, int>>, quegreater> que = adjnode[init];
        int add = -1, flag = 0;
        int dnewv = replace(key, init, -1);
        if (dnewv < value) {
            cout << "lkh  delet v: " << value << endl;
            value = dnewv;
            key = mm.getstring();
            flag = 1;
        }
        while (!que.empty() && flag == 0) {
            add = que.top().first;
            que.pop();
            if (k - forbiden[add] < 3)
                continue;
            int newv = replace(key, init, add);
            if (kk == 1)
                forbiden[add] = k;
            if (newv < value) {
                cout << "lkh v: " << value << endl;
                value = newv;
                key = mm.getstring();
                flag = 1;
                break;
            }
        }
        if (flag == 0)
            return make_pair(-1, value);
        return make_pair(add, value);
    };


    int choosepoint(string &key) {
        int init = -1;
        for (int i = candidates.size() - 1; i >= 0; i--)
            if (key[candidates[i].first] == '1' && k - forbiden[candidates[i].first] > 3) {
                init = candidates[i].first;
                break;
            }
        if (init < 0)
            cout << "< 0" << endl;
        return init;
    }

    int klkh(string &key, int &value) {
        while (true) {
            k = 0;
            cout << "init v is " << value << endl;
            int add = -1;
            int init = choosepoint(key);
            while (true) {
                pair<int, int> pp = lkh(key, value, init);
                value = pp.second;
                add = pp.first;
                if (add >= 0)
                    break;
                init = choosepoint(key);
                if (init < 0)
                    return 1;
            }
            cout << "choose init" << init << endl;
            for (int i = 0; i < 1000; i++) {
                priority_queue<pair<int, int>, vector<pair<int, int>>, quegreater> dque = dejnode[add];
                add = -1;
                while (!dque.empty()) {
                    int delet = dque.top().first;
                    dque.pop();
                    if (k - forbiden[delet] < 3)
                        add = lkh(key, value, delet).first;
                    if (add >= 0)
                        break;
                }
                k++;
                if (add < 0)
                    break;
            }
        }
        return 0;
    }


    string xjbh(string key) {
        set<int> exin;
        int initv = mm.TaoMCMf(key);
        vector<Edge> bak = mm.edges;
        for (int i = 0; i < key.size(); i++) {
            if (i == 9 && key[i] == '1' && exin.find(i) == exin.end()) {
                vector<Edge> dedge;
                int best = initv;
                int delta = mm.rotate_delete(i, bak);
                int value = mm.get_cost_allocation().first;
                key[i] = '0';
                cout << "really value is" << mm.TaoMCMf(key) << endl;
                key[i] = '1';
                int shouldv = mm.get_cost_allocation().first;
                if (shouldv != value) {
                    cout << "erro when deleting " << i << endl;
                    cout << delta + initv - mm.cdn_cost << endl;
                    cout << shouldv << " " << value << endl;
                }
                dedge=mm.edges;
                for(int j=0;j<mm.s;j++)
                    if(key[j]=='0'&&mm.maxnode[j]>0)
                    {
                        string kk=key;
                        kk[i]='0';
                        kk[j]='1';
                        //int ivv=mm.TaoMCMf(kk);if(ivv<best)key=kk,best=mm.get_cost_allocation().first;
                        int ivv=mm.rotate(j,dedge);

                        if(ivv<INFCOST&&mm.get_cost_allocation().first<best)
                            key=kk,best=ivv+shouldv+mm.cdn_cost;
                    }
                initv=best;
                int check=mm.TaoMCMf(key);
                cout<<"check edge flow"<<endl;

                bak=mm.edges;
                cout<<"a ex:"<<check<<" "<<best<<endl;
            }
        }
        return key;
    }

    int improvesearch(string key) {
        mm.reset();
        vector<Edge> gbak = mm.edges;
        int bakv = mm.TaoMCMf(key);
        vector<Edge> bak = mm.edges;
        for (int i = 0; i < mm.edges.size(); i += 2)
            if (mm.edges[i ^ 1].bandwith > 0)
                eva[i / 2] = make_pair(i / 2, mm.edges[i].unit_price);
            else
                eva[i / 2] = make_pair(i / 2, 0);
        for (int i = 0; i < mm.g[mm.s].size(); i++)
            if (mm.edges[mm.g[mm.s][i] ^ 1].bandwith > 0)
                eva[mm.g[mm.s][i] / 2].second = (double) mm.cdn_cost / (double) mm.edges[mm.g[mm.s][i] ^ 1].bandwith;
        sort(eva.begin(), eva.end(), ppcmp());
        int maxprice = eva[0].second;
        double totalprice = 0;
        int count = 0;
        for (int i = 0; i < eva.size(); i++) {
            totalprice += eva[i].second;
            if (eva[i].second > 0)
                count++;
            cout << eva[i].second << " ";
        }
        cout << endl;
        double average = totalprice / count;
        cout << average << endl;
        for (int i = 0; i < eva.size(); i++)
            if (eva[i].first >= mm.g[mm.s][0] / 2) {
                mm.reset();
                double wb = bak[2 * eva[i].first + 1].bandwith;
                if (wb != 121)
                    continue;
                cout << "first wb is:" << wb << endl;
                for (int j = 0; j < eva.size(); j++) {
                    if (eva[j].first >= mm.g[mm.s][0] / 2) {
                        if (i == j)
                            mm.setcoste(eva[j].first, mm.cdn_cost), mm.setcape(eva[j].first, 0);
                        if (eva[j].second > 0)
                            if (eva[j].second < 0)
                                mm.setcoste(eva[j].first, eva[j].second), mm.setcape(eva[j].first, INFCAPACITY);
                            else
                                mm.setcoste(eva[j].first, mm.cdn_cost), mm.setcape(eva[j].first, 0);

                        else if (mm.maxflow[eva[j].first - mm.g[mm.s][0] / 2] > 0)
                            mm.setcoste(eva[j].first, (double) mm.cdn_cost / (double) mm.maxflow[eva[j].first -
                                                                                                 mm.g[mm.s][0] /
                                                                                                 2]), mm.setcape(
                                    eva[j].first, INFCAPACITY);
                        else
                            mm.setcoste(eva[j].first, mm.cdn_cost), mm.setcape(eva[j].first, 0);

                    } else {
                        double ss = min(bak[eva[j].first * 2].bandwith * 1.0, wb);
                        if (eva[j].second < average)
                            mm.setcoste(eva[j].first, ss * bak[2 * eva[j].first].unit_price),
                                    mm.setcape(eva[j].first, bak[2 * eva[j].first + 1].bandwith);
                        else
                            mm.setcoste(eva[j].first, mm.cdn_cost), mm.setcape(eva[j].first, 0);
                    }
                }
                cout << wb << endl;
                cout << "value is " << mm.MAXCMCMF(wb) << endl;
                string kk = key;
                string newkey = mm.getstring();
                cout << "before key" << endl << key << endl;
                cout << "add key " << endl << newkey << endl;
                for (int k = 0; k < key.size(); k++)
                    if (key[k] == '1')
                        key[k] = '1';
                //key[eva[i].first-mm.g[mm.s][0]/2]='0';
                mm.edges = gbak;
                key[68] = '0';
                //key[241]='1';
                //key[76]='1';
                //key[77]='0';
                //key[168]='1';
                //key[177]='0';
                key[271] = '1';
                //key[190]='0';
                int newv = mm.TaoMCMf(key);
                cout << newv << endl;
                if (newv < bakv) {
                    bakv = newv;
                    cout << "change vlaue" << newv << endl;
                }
                break;
            }

        cout << endl;
    }

    void ppk(string key) {
        for (int i = 0; i < key.size(); i++)
            cout << key[i] << "";
        cout << endl;
    }
};
