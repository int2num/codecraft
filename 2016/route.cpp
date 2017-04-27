
#include <bits/stdc++.h>
#include "route.h"
#include "lib_record.h"



using namespace std;


namespace BAOLI
{
    #define M 5010
    #define N 610
    #define INF 1000000
    struct Edge
    {
        int id;
        int head;
        int tail;
        int cost;
    };

    Edge e[M];
    int vertex[N];


    set<int> v;
    int s, t;
    int edge_num;
    int vertexVisit[N];
    bool _flag;
    bool ltEdge(const Edge &a, const Edge &b)
    {
        return a.head < b.head || (a.head == b.head && a.tail < b.tail);
    }

    void parse(char *topo[5000], int edge_num, char *demand)
    {
        BAOLI::edge_num = edge_num;
        for(int i = 0; i < edge_num; i++)
            sscanf(topo[i], "%d,%d,%d,%d", &e[i].id, &e[i].head, &e[i].tail, &e[i].cost);

        sort(e, e+edge_num, ltEdge);
        int num = 0;
        while(*demand != ',')
            num = num*10 + *(demand++) - '0';
        s = num;
        num = 0;
        demand++;
        while(*demand != ',')
            num = num*10 + *(demand++) - '0';
        t = num;
        num = 0;
        demand++;
        while((*demand <= '9' && *demand >= '0' ) || *demand == '|')
        {
            if(*demand == '|')
            {
                v.insert(num);
                num = 0;
            }
            else
                num = num*10 + *demand - '0';
            demand++;
        }
        v.insert(num);

        for(int i = 0; i < N; i++)
            vertex[i] = edge_num+1;
        vertex[e[0].head] = 0;
        for(int i = 1; i < edge_num; i++)
        {
            if(e[i-1].head != e[i].head)
            {
                 vertex[e[i].head] = i;
                 vertex[e[i].head+1] = edge_num;
            }
        }
        for(int i = N-1; i; i--)
            if(vertex[i] == edge_num+1)
                vertex[i] = vertex[i+1];
    }

    void showPath(vector<int> &path, char *str)
    {
        printf("%s", str);
        int len = 0;
        for(int i = 0; i < path.size(); i++)
        {
            len += e[path[i]].cost;
            printf("%d ", e[path[i]].head);
            if(i == path.size()-1)
                printf("%d ", e[path[i]].tail);
        }
        printf("\nlen:%d\n", len);

    }

    vector<int> pathAns;
    int minLen = INT_MAX;
    time_t tb;
    double TIME_OUT = 0.1;

    bool solv(int cur, vector<int> &path, int len, set<int> &cv)
    {
        time_t te = clock();
        if((te-tb*1.0)/CLOCKS_PER_SEC > TIME_OUT)
        {
            _flag = false;
            return false;
        }

        if(len > minLen) return false;

        if(cur == t)
        {
            if(cv != v) return false;
            if(len < minLen)
            {
                //showPath(path, "");
                minLen = len;
                pathAns = path;
            }
            return true;
        }
        if(vertexVisit[cur])
            return false;
        else
            vertexVisit[cur] = 1;
        if(v.find(cur) != v.end())
            cv.insert(cur);
        bool ans = false;
        for(int i = vertex[cur]; i < vertex[cur+1]; i++)
        {
            if(vertexVisit[e[i].tail]) continue;
            path.push_back(i);
            if(solv(e[i].tail, path, len+e[i].cost, cv))
                ans = true;
            path.erase(--path.end());
        }
        vertexVisit[cur] = 0;
        cv.erase(cur);
        return ans;
    }

    pair<int, vector<int>> solvx(char *topo[5000], int edge_num, char *demand, bool &flag)
    {
        _flag = true;
        tb = clock();
        parse(topo, edge_num, demand);
        vector<int> t;
        memset(vertexVisit, 0, sizeof(vertexVisit));
        int len = 0;
        set<int> cv;
        solv(s, t, len, cv);
        vector<int> ans;
        for(int i = 0; i < pathAns.size(); i++)
            ans.push_back(e[pathAns[i]].id);
        //showPath(pathAns, " ");
        flag = _flag;
        return make_pair(minLen, ans);
    }
}

namespace PATH
{
    vector<pair<int,pair<int,int> > > gra[160005];
    int pre[160005];
    int pid[160005];
    int dis[160005];
    int nn;
    bool vis[160005];
    int inV[160005];
    bool del[160005];
    inline int Dijkstra(int s)
    {
        int mdis,mx,k=1,to,w;
        priority_queue<pair<int,int>,vector<pair<int,int> >,greater<pair<int,int> > > q;
        for(int i=0;i<nn;i++) dis[i]=0x3fffffff,vis[i]=0;
        q.push(make_pair(0,s));dis[s]=0;pre[s]=-1;
        while(!q.empty())
        {
            mdis=q.top().first;
            mx=q.top().second;q.pop();
            if(vis[mx]) continue;
            if(inV[mx]&&mx!=s) return mx;
            vis[mx]=1;
            for(int i=0;i<gra[mx].size();i++)
                if(dis[to=gra[mx][i].first]>(w=gra[mx][i].second.first)+mdis)
                {
                    if(del[to]==0&&vis[to]==0)
                    {
                        pre[to]=mx;
                        pid[to]=gra[mx][i].second.second;
                        dis[to]=w+mdis,q.push(make_pair(dis[to],to));
                    }
                }
        }
        return -1;
    }
    vector<int> solv(int n,int s,int t,vector<int>&id,vector<int>&px,vector<int>&py,vector<int>&pw,vector<int>&V, bool &flagD)
    {
        flagD = false;
        nn=n;
        int m=id.size();
        int k=V.size();
        for(int i=0;i<n;i++)
            gra[i].clear(),inV[i]=del[i]=0;
        for(int i=0;i<k;i++)
            inV[V[i]]=1;
        for(int i=0;i<m;i++)
        {
            if(inV[py[i]])gra[px[i]].push_back(make_pair(py[i],make_pair(pw[i],id[i])));
            else if(py[i]==t)gra[px[i]].push_back(make_pair(py[i],make_pair(pw[i]+1000,id[i])));
            else gra[px[i]].push_back(make_pair(py[i],make_pair(pw[i],id[i])));
        }
        inV[t]=1;
        vector<int>ans;
        vector<int>tem;
        int sum=0;
        while(s!=t)
        {
            tem.clear();
            if(sum==k)
            {
                for(int i=0;i<n;i++)
                    for(int j=0;j<gra[i].size();j++)
                    {
                        if(gra[i][j].first==t)
                            gra[i][j].second.first-=1000;
                    }
            }
            int c=Dijkstra(s);
            if(c==t&&sum!=k)
            {
//                    ans.clear();
                return ans;
            }
            int f=c;
            if(c==-1)
            {
//                    ans.clear();
                return ans;
            }
            while(c != s)
            {
                tem.push_back(pid[c]);
                c=pre[c];
                del[c]=1;
            }
            reverse(tem.begin(),tem.end());
            ans.insert(ans.end(),tem.begin(),tem.end());
            s=f;
            sum++;
        }
        flagD = true;
        return ans;
    }

    void parse(char *topo[5000], int edge_num, char *demand, vector<vector<int> > &data, int &s, int &t, int &n)
    {
        n = 0;
        vector<int> vid, vx, vy, vw, vv;
        for(int i = 0; i < edge_num; i++)
        {
            int a, b, c, d;
            sscanf(topo[i], "%d,%d,%d,%d", &a, &b, &c, &d);
            vid.push_back(a);
            vx.push_back(b);
            vy.push_back(c);
            n = max(b, max(n, c));
            vw.push_back(d);
        }
        n++;
        int num = 0;
        while(*demand != ',')
            num = num*10 + *(demand++) - '0';
        s = num;
        num = 0;
        demand++;
        while(*demand != ',')
            num = num*10 + *(demand++) - '0';
        t = num;
        num = 0;
        demand++;
        while((*demand <= '9' && *demand >= '0' ) || *demand == '|')
        {
            if(*demand == '|')
            {
                vv.push_back(num);
                num = 0;
            }
            else
                num = num*10 + *demand - '0';
            demand++;
        }
        vv.push_back(num);
        data.push_back(vid);
        data.push_back(vx);
        data.push_back(vy);
        data.push_back(vw);
        data.push_back(vv);
    }

    int calLen(vector<int> &path, vector<int> &vc)
    {
        if(path.size() == 0) return INT_MAX;
        int ans = 0;
        for(int i = 0; i < path.size(); i++)
            ans += vc[path[i]];
        return ans;

    }
    pair<int, vector<int>> solvx(char *topo[5000], int edge_num, char *demand)
    {
        srand(time(nullptr));
        time_t ts = clock();
        vector<vector<int> > data;
        int s, t, n;
        parse(topo, edge_num, demand, data, s, t, n);
        vector<int> w = data[3], ans;
        int ansLen = INT_MAX;
        for(int i = 0; i < 1000000; i++)
        {
            if(i%5 == 0) w = data[3];
            if(((double)clock()-ts)/CLOCKS_PER_SEC > 8) break;
            bool flag = false;
            vector<int> anst = PATH::solv(n, s, t, data[0], data[1], data[2], w, data[4], flag);
            if(flag)
            {
                int t = calLen(anst, data[3]);
                if(t < ansLen)
                {
                    printf("Optimal solution: %d\n", t);
                    ansLen = t;
                    ans = anst;
                }
            }
            vector<int> anst2 = anst;
            anst = PATH::solv(n, t, s, data[0], data[2], data[1], w, data[4], flag);
            if(flag)
            {
                int t = calLen(anst, data[3]);
                if(t < ansLen)
                {
                    printf("Optimal solution reverse: %d\n", t);
                    ansLen = t;
                    reverse(anst.begin(), anst.end());
                    ans = anst;
                }
            }
            for(int j = 0; j < anst.size(); j++)
                w[anst[j]] += rand()%3;
            for(int j = 0; j < anst2.size(); j++)
                w[anst2[j]] += rand()%3;
        }
        return make_pair(ansLen, ans);
    }

}

//    pair<int, vector<int>> solvx(char *topo[5000], int edge_num, char *demand)
//    {
//        pair<int, vector<int>> ans1 = BAOLI::solvx(topo, edge_num, demand);
//        pair<int, vector<int>> ans2 = PATH::solvx(topo, edge_num, demand);
//        if(ans1.first < ans2.first)
//            return ans1;
//        else return ans2;
//    }


// wenxiong dij+搜索
namespace pAth
{
    time_t ts;
    vector<pair<int,pair<int,int> > > gra[160005];
    int dis[160005];
    int pre[160005];
    int pid[160005];
    int nn;
    bool vis[160005];
    bool _flag;
    inline int Dijkstra(int s)
    {
        int mdis,mx,k=1,to,w;
        priority_queue<pair<int,int>,vector<pair<int,int> >,greater<pair<int,int> > > q;
        for(int i=0;i<nn;i++) dis[i]=0x3fffffff,vis[i]=0;
        q.push(make_pair(0,s));dis[s]=0;
        while(!q.empty())
        {
            mdis=q.top().first;
            mx=q.top().second;q.pop();
            if(vis[mx])continue;
            vis[mx]=1;
            for(int i=0;i<gra[mx].size();i++)
                if(dis[to=gra[mx][i].first]>(w=gra[mx][i].second.first)+mdis)
                {
                    pre[to]=mx;
                    pid[to]=gra[mx][i].second.second;
                    dis[to]=w+mdis,q.push(make_pair(dis[to],to));
                }
        }
    }
    vector<int>pa[61][61],pp[61][61];
    int tt,kk;
    int len[61][61];
    int fans;
    bool use[605];
    int sum;
    vector<int>po;
    int te[61],lk;
    int sori[61][61];
    int cnt;
    inline bool ok(int x,int y)
    {
        for(int i=0;i<pp[x][y].size();i++)
            if(use[pp[x][y][i]])
                return 0;
        return 1;
    }
    inline void modi(int x,int y,int f)
    {
        for(int i=0;i<pp[x][y].size();i++)
            use[pp[x][y][i]]=f;
    }
    void dfs(int s,int nd)
    {
        cnt++;
        if(s==tt)
        {
            if(fans>nd)
            {
                fans=nd;
                po.clear();
                for(int i=0;i<lk;i++)
                    po.push_back(te[i]);
            }
        }
        else
        {
            for(int i=0;i<kk;i++)
            {
                int x=sori[s][i];
                if(s!=x&&ok(s,x))
                {
                    if(x==tt&&sum<kk-2)continue;
                    modi(s,x,1);
                    sum++;
                    te[lk++]=x;
                    if(nd+len[s][x]+len[x][tt]<fans)
                        dfs(x,nd+len[s][x]);
                    if((clock()-(double)ts)/CLOCKS_PER_SEC > 1.8)
                    {
                        _flag = false;
                        return;
                    }
                    modi(s,x,0);
                    sum--;
                    lk--;
                }
            }
        }
    }
    vector<int> solv(int n,int s,int t,vector<int>&id,vector<int>&px,vector<int>&py,vector<int>&pw,vector<int>&v, bool &flag)
    {
        _flag = true;
        ts =clock();
        nn=n;
        int m=px.size();
        for(int i=0;i<n;i++)
            gra[i].clear();
        for(int i=0;i<m;i++)
            gra[px[i]].push_back(make_pair(py[i],make_pair(pw[i],id[i])));
        v.push_back(t);
        v.push_back(s);
        kk=v.size();
        for(int i=0;i<kk;i++)
            for(int j=0;j<kk;j++)
                pa[i][j].clear(),pp[i][j].clear();
        for(int i=0;i<kk;i++)
        {
            int x=v[i];
            Dijkstra(x);
            for(int j=0;j<kk;j++)
            {
                int y=v[j],t=v[j];
                len[i][j]=dis[y];
                if(dis[y]==0x3fffffff)
                    continue;
                while(t!=x)
                {
                    pa[i][j].push_back(pid[t]);
                    pp[i][j].push_back(t);
                    t=pre[t];
                }
                reverse(pa[i][j].begin(),pa[i][j].end());
            }
        }
        for(int i=0;i<kk;i++)
        {
            pair<int,int>c[61];
            for(int j=0;j<kk;j++)
                c[j]=make_pair(len[i][j],j);
            sort(c,c+kk);
            for(int j=0;j<kk;j++)
                sori[i][j]=c[j].second;
        }
        tt=kk-2;
        memset(use,0,sizeof(use));
        sum=0;
        po.clear();
        use[v[kk-1]]=1;
        lk=0;
        fans=1000000;
        cnt=0;
        te[lk++]=kk-1;
        dfs(kk-1,0);
        vector<int>ans;
        if(po.size()==0) return ans;
        for(int i=0;i<po.size()-1;i++)
        {
            int x=po[i],y=po[i+1];
            for(int j=0;j<pa[x][y].size();j++)
                ans.push_back(pa[x][y][j]);
        }
        flag = _flag;
        return ans;
    }

    void parse(char *topo[5000], int edge_num, char *demand, vector<vector<int> > &data, int &s, int &t, int &n)
    {
        n = 0;
        vector<int> vid, vx, vy, vw, vv;
        for(int i = 0; i < edge_num; i++)
        {
            int a, b, c, d;
            sscanf(topo[i], "%d,%d,%d,%d", &a, &b, &c, &d);
            vid.push_back(a);
            vx.push_back(b);
            vy.push_back(c);
            n = max(b, max(n, c));
            vw.push_back(d);
        }
        n++;
        int num = 0;
        while(*demand != ',')
            num = num*10 + *(demand++) - '0';
        s = num;
        num = 0;
        demand++;
        while(*demand != ',')
            num = num*10 + *(demand++) - '0';
        t = num;
        num = 0;
        demand++;
        while((*demand <= '9' && *demand >= '0' ) || *demand == '|')
        {
            if(*demand == '|')
            {
                vv.push_back(num);
                num = 0;
            }
            else
                num = num*10 + *demand - '0';
            demand++;
        }
        vv.push_back(num);
        data.push_back(vid);
        data.push_back(vx);
        data.push_back(vy);
        data.push_back(vw);
        data.push_back(vv);
    }

    int calLen(vector<int> &path, vector<int> &vc)
    {
        if(path.size() == 0) return INT_MAX;
        int ans = 0;
        for(int i = 0; i < path.size(); i++)
            ans += vc[path[i]];
        return ans;

    }

    void showPath(vector<int> &path, vector<int> &vc, vector<int> &px, vector<int> &py)
    {
        printf("len:%d\n", calLen(path, vc));
        for(int i = 0; i < path.size(); i++)
        {
            printf("%d ", px[path[i]]);
            if(i == path.size()-1)
            {
                printf("%d\n", py[path[i]]);
            }
        }
    }
    pair<int, vector<int>> solvx(char *topo[5000], int edge_num, char *demand, bool &flag)
    {
        vector<vector<int> > data;
        int s, t, n;
        parse(topo, edge_num, demand, data, s, t, n);
        vector<int> ans = pAth::solv(n, s, t, data[0], data[1], data[2], data[3], data[4], flag);
//        showPath(ans, data[3], data[1], data[2]);
        return make_pair(calLen(ans, data[3]), ans);
    }
}

	 namespace ANT
	 {
		#define N 900
		#define M 8000
		#define INF 1000000
		using namespace std;
		typedef struct Ed{
			int head, tail;
			double weight;
			int mark,enable,bacw;
			double Pheromone;
		}Edge;
		typedef struct Gr{
			Edge* incL[M+10];
			int p[N][200];
                        int f[N][200];
		}Graph;
		void Graph_addEdge(Edge* e,int s, int t, double w,int mark,int bacw,int enable,double Pheromone){
			e->tail = t;
			e->head = s;
			e->weight = w;
			e->mark=mark;
			e->enable=enable;
			e->bacw=bacw;
			e->Pheromone=Pheromone;

		}
		int s=-1,t=-1,n_num=0;
		vector<int> V;
                int flag[N+10];
		int pre[N];
		int peg[N];
		int ISinV[N]={0};
		Graph *G = (Graph *)malloc(sizeof(Graph));
		const double ALPHA=3;
		const double BETA=0.05;//0.2 不出14.
		const double GAMA=0.8;//ss
		const double YITA=1;//jin
		const double NAMENDA=500;//s
		const double ROU=50;//sss
		const double PI=-0.0;//jin
		const double CTA=5.0;
		const double YIBUCNO=0.0;
		const bool MaxMinEnable=true;
		const double Low=0.0;
		const double High=30;
		const double InitPheromone=1;
		double OptimalWeight=INF;
                //dijstra
	        class Heap{
		private:
			Edge *h[N+10];
			int post[N+10];
			int nodeNum;
			void fix(int fixID){
				while(fixID > 1){
					int fat = fixID/2;
					if(h[fixID]->weight < h[fat]->weight){
						swap(post[h[fixID]->tail], post[h[fat]->tail]);
						swap(h[fixID], h[fat]);
						fixID /= 2;
					}
					else
						break;
				}
			}
		public:
			Heap(){
				nodeNum = 0;
			}
			void push(int vertID, int w){
				nodeNum++;
				h[nodeNum] = (Edge *)malloc(sizeof(Edge));
				h[nodeNum]->tail = vertID;
				h[nodeNum]->weight = w;
				post[vertID] = nodeNum;
				fix(nodeNum);
			}
			void update(int vertID, int w){
				int p = post[vertID];
				h[p]->weight = w;
				fix(p);
			}
			int pop(){
				int ret = h[1]->tail;
				free(h[1]);
				h[1] = h[nodeNum--];
				if(nodeNum > 0)
					post[h[1]->tail] = 1;

				int cur = 1;
				while(cur * 2 <= nodeNum){
					int son = cur * 2;
					if(son + 1 <= nodeNum && h[son]->weight > h[son+1]->weight)
						son = son + 1;
					if(h[cur]->weight > h[son]->weight){
						swap(post[h[cur]->tail], post[h[son]->tail]);
						swap(h[cur], h[son]);
						cur = son;
					}
					else
						break;
				}

				return ret;
			}
			int empty(){
				if(nodeNum <= 0)
					return 1;
				else
					return 0;
			}
			void display(){
				printf("********************\n");
				for(int i = 1; i <= nodeNum; i++)
					printf("%d %d\n", h[i]->tail, h[i]->weight);
				for(int i = 0; i < N; i++)
					printf("pot_%d: %d\n", i, post[i]);
				printf("********************\n");
			}
			~Heap(){
				while(nodeNum > 0)
					pop();
			}
	};
	void dijkstraWithHeap(Graph *G, int s, int t, int d[],int n_num){

		//printf("inner ters1\n");
		for(int i = 0; i < n_num; i++)
			if(i == s)
				d[i] = 0;
			else
				d[i] = INF;
		for(int i = 0; i < n_num; i++)
			{flag[i] = 0;
		         pre[i]=-1;
		         peg[i]=-1;
		         }
	  //printf("inner ters2\n");
		//t_s = clock();
		int cur=s;
		Heap heap;
		for(int i = 0; i < n_num; i++)
			heap.push(i, d[i]);
		do{

			//printf("poping\n");
			int cur = heap.pop();
		        flag[cur]=1;
		        if(cur == t)
				break;
			//printf("pop end\n");
			//heap.display();

			for(int i = 0;G->p[cur][i]>=0;i++){
		                //printf("p[cur][i]:%d\n",G->p[cur][i]);
				Edge* e = G->incL[G->p[cur][i]];
				if(e->enable>0&&d[e->tail]>d[e->head]+e->weight&&flag[e->tail]==0){
					d[e->tail] = d[e->head] + e->weight;
					//printf("updating %d %d\n", e.head, d[e.head]);
					heap.update(e->tail, d[e->tail]);
		                        pre[e->tail]=e->head;
		                        peg[e->tail]=G->p[cur][i];
					//printf("update end\n");
				}
			}
			//heap.display();
		}while(!heap.empty());
		//t_t = clock();

	   }


		//返回指定范围内的随机整数

		double Adapter(double delta)
		{if(MaxMinEnable)
		  {if(delta<Low)
		    return Low;
		   if(delta>High)
		    return High;}
		 return delta;
		}
		int rnd(int nLow,int nUpper)
		{
		    return nLow+(nUpper-nLow)*rand()/(RAND_MAX+1);
		}

		//返回指定范围内的随机浮点数
		double rnd(double dbLow,double dbUpper)
		{
		    double dbTemp=rand()/((double)RAND_MAX+1.0);
		    return dbLow+dbTemp*(dbUpper-dbLow);
		}

		//返回浮点数四舍五入取整后的浮点数
		double ROUND(double dbA)
		{
		    return (double)((int)(dbA+0.5));
		}

		class Ant
		{
		  public:
		     int Number;
		     int Findt=0;
		     vector<int> Path;
		     vector<int> EdgeID;
		     double PathLen;
		     vector<int> VistedV;
		     vector<int> VistedNode;
		     vector<int> NodeSequence;
		     int BadEdge[N][10];
		     int CurrentPosition;
		     int IsVisted[N]={0};
		  public:
		   Ant(int num)
		   {this->Number=num;
		    Init();
		   }
		   void Init()
		   {Findt=0;
		    PathLen=0;
		    memset(BadEdge,0,sizeof(BadEdge));
		    NodeSequence.push_back(s);
		    VistedNode.push_back(s);
		    CurrentPosition=s;
		    IsVisted[s]=1;
		   }

		   int ChooseNext()
		   {   //printf("in choose next!\n");
		       double TotalPheromone=0;
		       double Prob[200]={0};
		       for(int i=0;G->p[CurrentPosition][i]>=0;i++)
			  {//printf("in 1\n");
			   Edge* edge=G->incL[G->p[CurrentPosition][i]];
			   if(IsVisted[edge->tail]==0&&BadEdge[CurrentPosition][i]==0)//!(VistedV.size()<V.size()&&edge->tail==t)
			     { TotalPheromone=TotalPheromone+pow(edge->Pheromone,ALPHA)*pow(1.0/edge->weight,BETA);
			       Prob[i]=pow(edge->Pheromone,ALPHA)*pow(1.0/edge->weight,BETA);
			     }
			  }

		       if(TotalPheromone==0)
			   return -2;
		       else
			  {
			   double dbTemp=0.0;
			   dbTemp=rnd(0.0,TotalPheromone);
			   for(int i=0;G->p[CurrentPosition][i]>=0;i++)
			      {Edge* edge=G->incL[G->p[CurrentPosition][i]];
			       if(IsVisted[edge->tail]==0&&BadEdge[CurrentPosition][i]==0)
			       {    dbTemp=dbTemp-Prob[i];
				    if (dbTemp < 0.0)
				      {EdgeID.push_back(i);
				       return G->p[CurrentPosition][i];}
				}
			      }
			  }

		    }
		   void move()
		   {  int NextEdge=-1;
		      //printf("in move()!\n");
		      if(CurrentPosition==t)
			  NextEdge=-1;
		      else
			  NextEdge=this->ChooseNext();
		      //printf("in move()2\n");
		      if(this->PathLen>OptimalWeight )
			  {NextEdge=-1;
			   //for(int z=0;z<this->Path.size();z++)
			      // G->incL[this->Path[z]]->Pheromone=Adapter(G->incL[this->Path[z]]->Pheromone-YIBUCNO);
			  }
		      // printf("in move()3\n");
		      if(NextEdge>=0)
			  {
			   int NextNode=G->incL[NextEdge]->tail;
			   PathLen=PathLen+G->incL[NextEdge]->weight;
			   Path.push_back(NextEdge);

			   VistedNode.push_back(NextNode);
			   NodeSequence.push_back(NextNode);
			  // printf("add %d\n",NextNode);
			   if(ISinV[NextNode]==1)
			       VistedV.push_back(NextNode);
			   CurrentPosition=NextNode;
			   IsVisted[NextNode]=1;
			   if(NextNode==t)
			       Findt=1;

			  // G->incL[NextEdge]->Pheromone=Adapter(G->incL[NextEdge]->Pheromone-YITA);
			  }
		    //  printf("in move()4\n");
		      if(NextEdge==-1)//成功得到解或者超重，回到s
			 { PathLen=0;
			   Path.clear();
			   VistedNode.clear();
			   NodeSequence.clear();
			   VistedV.clear();
			   memset(IsVisted,0,N*sizeof(int));
			   this->Init();
			 }
		      // printf("in move()5\n");
		      if(NextEdge==-2)//回朔。
			{
			//printf("ant%d in back up procedue\n",this->Number);
			//printf("set good edge down %d\n",NodeSequence[NodeSequence.size()-1]);
			for(int i=0;i<10;i++)
			   BadEdge[NodeSequence[NodeSequence.size()-1]][i]=0;
			Path.pop_back();
			IsVisted[VistedNode[VistedNode.size()-1]]=0;
			if(ISinV[VistedNode[VistedNode.size()-1]])
			    VistedV.pop_back();
			NodeSequence.pop_back();
			VistedNode.pop_back();
			CurrentPosition=NodeSequence[NodeSequence.size()-1];
			//printf("Current position is :%d\n",CurrentPosition);
			BadEdge[CurrentPosition][EdgeID[EdgeID.size()-1]]=1;
			EdgeID.pop_back();
		       // printf("set %d,%d to be :",CurrentPosition,EdgeID[EdgeID.size()-1)
			//for(int r=0;r<N;r++)
			     //if(IsVisted[r]>0)
				  // printf("%d-",r);
			}
			 //printf("in move()6\n");
		    }
		};

               void Check(vector<int> OptimalSequence)
                {int Isv[N]={0};
                 int count=0;
                 for(int i=0;i<OptimalSequence.size();i++)
                 {if(Isv[G->incL[OptimalSequence[i]]->tail]<=0)
                      {Isv[G->incL[OptimalSequence[i]]->tail]=1;
                       if(ISinV[G->incL[OptimalSequence[i]]->tail]>0)
                           count++;
                       }
                  else
                    {printf("\nerro 1overlap node!%d\n",G->incL[OptimalSequence[i]]->tail);
                     break;
                    }

                  }
                  if(count<V.size())
                      printf("erro unreach node!\n");
                  else
                    printf("correct!\n");

                 }


            vector<int> DjiDirect2(double &weight ,int edge_num,vector<int> EdgeSequence)
                {vector<int> NodeSquence;
                 weight=0;
                 for(int i=0;i<EdgeSequence.size();i++)
                      NodeSquence.push_back(G->incL[EdgeSequence[i]]->head);
                   NodeSquence.push_back(t);
                 reverse(NodeSquence.begin(),NodeSquence.end());
                 reverse(EdgeSequence.begin(),EdgeSequence.end());
                 vector<int> OptimalSequence;
                 vector<int> SubOptimalSequence;
                 vector<int>SubNodeSequence;
                // printf("hrea 1\n");
                 int ISinNode[N]={0};
                 for(int i=0;i<NodeSquence.size();i++)
                     ISinNode[NodeSquence[i]]=1;
                 for(int i=0;i<NodeSquence.size();i++)
                     if(ISinV[NodeSquence[i]]>0||NodeSquence[i]==s)
                          for(int k=0;G->p[NodeSquence[i]][k]>=0;k++)
                                     G->incL[G->p[NodeSquence[i]][k]]->enable=0;
                 int start=0;
                 int end;
                 int d[n_num];
                 //double weight=0.0;
                  //printf("hrea 1\n");
                 for(int i=1;i<NodeSquence.size();i++)
                  {

                    if(ISinV[NodeSquence[i]]||NodeSquence[i]==s)
                      { end=i;
                        for(int k=0;G->p[NodeSquence[end]][k]>=0;k++)
                                     G->incL[G->p[NodeSquence[end]][k]]->enable=1;
                        //printf("\n%d< ----------------------%d\n",NodeSquence[start],NodeSquence[end]);
                        dijkstraWithHeap(G,NodeSquence[end],NodeSquence[start],d,n_num);
                        int f=pre[NodeSquence[start]];
                        SubOptimalSequence.clear();
                        SubNodeSequence.clear();
                        SubOptimalSequence.push_back(peg[NodeSquence[start]]);

                       while(f!=NodeSquence[end]&&f!=-1)
                       {
                        SubNodeSequence.push_back(f);
                       // printf("%d<-",f);
                        SubOptimalSequence.push_back(peg[f]);
                        f=pre[f];
                       }
                       if(f!=-1)
                         { SubNodeSequence.push_back(f);
                            //printf("%d<-",f);
                         }
                        if(f==-1)
                           { //printf("dao gua le!!!!!!\n");
                             OptimalSequence.clear();
                             weight=INF;
                            return OptimalSequence;

                            }
                        else
                            {//printf("mei gua!\n");
                             for(int j=0;j<SubNodeSequence.size();j++)
                                 {//printf("sail %d\n",SubNodeSequence[j]);
                                  for(int k=0;G->p[SubNodeSequence[j]][k]>=0;k++)
                                       G->incL[G->p[SubNodeSequence[j]][k]]->enable=-2;
                                  }

                             for(int k=0;G->p[NodeSquence[start]][k]>=0;k++)
                                     G->incL[G->p[NodeSquence[start]][k]]->enable=-2;
                             for(int h=0; h<SubOptimalSequence.size();h++)
                                 {OptimalSequence.push_back(SubOptimalSequence[h]);
                                  weight=weight+G->incL[SubOptimalSequence[h]]->weight;
                                  }

                            }
                       int fe=peg[NodeSquence[end]];
                       start=end;
                       //printf("\n");
                      }

                  }
                  vector<int> OptimalNode;
                  for(int i=0;i<OptimalSequence.size();i++)
                     OptimalNode.push_back(G->incL[OptimalSequence[i]]->head);
                  OptimalNode.push_back(t);
              //    printf("in Direct2 optimal w is:%f\n",weight);
                 // Check(OptimalSequence);
                  // for(int y=0;y<OptimalSequence.size();y++)
                     //printf("%d->",OptimalSequence[y]);
              //    printf("\n");
                    //for(int y=0;y<OptimalSequence.size();y++)
                     // record_result(OptimalSequence[y]);
                   for(int i=0;i<edge_num;i++)
                      G->incL[i]->enable=1;
                   reverse(OptimalSequence.begin(),OptimalSequence.end());
                   reverse(OptimalNode.begin(),OptimalNode.end());
                   return OptimalSequence;
                }


               vector<int> DjiDirect(double& weight ,int edge_num,vector<int> EdgeSequence)
                {vector<int> NodeSquence;
                 weight=0;
                 for(int i=0;i<EdgeSequence.size();i++)
                     NodeSquence.push_back(G->incL[EdgeSequence[i]]->head);
                   NodeSquence.push_back(t);
                 //for(int i=0;i<NodeSquence.size();i++)
                   //  printf("->%d",NodeSquence[i]);
                 vector<int> OptimalSequence;
                 vector<int> SubOptimalSequence;
                 vector<int>SubNodeSequence;
                // printf("hrea 1\n");
                 int ISinNode[N]={0};
                 for(int i=0;i<NodeSquence.size();i++)
                     ISinNode[NodeSquence[i]]=1;
                 for(int i=0;i<NodeSquence.size();i++)
                     if(ISinV[NodeSquence[i]]>0||NodeSquence[i]==t)
                          for(int k=0;G->f[NodeSquence[i]][k]>=0;k++)
                                     G->incL[G->f[NodeSquence[i]][k]]->enable=0;
                 int start=0;
                 int end;
                 int d[n_num];
                 //double weight=0.0;
                  //printf("hrea 1\n");
                 for(int i=1;i<NodeSquence.size();i++)
                  {

                    if(ISinV[NodeSquence[i]]||NodeSquence[i]==t)
                      { end=i;
                        for(int k=0;G->f[NodeSquence[end]][k]>=0;k++)
                                     G->incL[G->f[NodeSquence[end]][k]]->enable=1;
                        //printf("\n%d ---------------------->%d\n",NodeSquence[start],NodeSquence[end]);
                        dijkstraWithHeap(G,NodeSquence[start],NodeSquence[end],d,n_num);
                        int f=pre[NodeSquence[end]];
                        SubOptimalSequence.clear();
                        SubNodeSequence.clear();
                        SubOptimalSequence.push_back(peg[NodeSquence[end]]);

                       while(f!=NodeSquence[start]&&f!=-1)
                       {
                        SubNodeSequence.push_back(f);
                        //printf("%d<-",f);
                        SubOptimalSequence.push_back(peg[f]);
                        f=pre[f];
                       }
                       if(f!=-1)
                         { SubNodeSequence.push_back(f);
                            //printf("%d<-",f);
                         }
                        if(f==-1)
                           {
                           // printf("gua le!!!!!!\n");
                            OptimalSequence.clear();
                            weight=INF;
                            return OptimalSequence;
                            }
                        else
                            {//printf("mei gua!\n");
                             for(int j=0;j<SubNodeSequence.size();j++)
                                 {//printf("sail %d\n",SubNodeSequence[j]);
                                  for(int k=0;G->f[SubNodeSequence[j]][k]>=0;k++)
                                       G->incL[G->f[SubNodeSequence[j]][k]]->enable=-2;
                                  }

                             for(int k=0;G->f[NodeSquence[end]][k]>=0;k++)
                                     G->incL[G->f[NodeSquence[end]][k]]->enable=-2;
                             for(int h=SubOptimalSequence.size()-1;h>=0;h--)
                                 {OptimalSequence.push_back(SubOptimalSequence[h]);
                                  weight=weight+G->incL[SubOptimalSequence[h]]->weight;
                                  }

                            }
                       int fe=peg[NodeSquence[end]];
                       start=end;
                       //printf("\n");
                      }

                  }
                  vector<int> OptimalNode;
                  for(int i=0;i<OptimalSequence.size();i++)
                     OptimalNode.push_back(G->incL[OptimalSequence[i]]->head);
                  OptimalNode.push_back(t);
                  //for(int i=0;i<OptimalNode.size();i++)
                    //printf("%d->",OptimalNode[i]);
                  //printf("int DIdir optimal w is:%f\n",weight);
                  //for(int y=0;y<OptimalSequence.size();y++)
                     //printf("%d->",OptimalSequence[y]);
                  //printf("\n");
                   //Check(OptimalSequence);
                  for(int i=0;i<edge_num;i++)
                      G->incL[i]->enable=1;
                  return OptimalSequence;


                }


               vector<int> DjiManage2(double &weight ,int edge_num,vector<int> EdgeSequence)
                {vector<int> NodeSquence;
                 weight=0;
                 for(int i=0;i<EdgeSequence.size();i++)
                      NodeSquence.push_back(G->incL[EdgeSequence[i]]->head);
                   NodeSquence.push_back(t);
                 reverse(NodeSquence.begin(),NodeSquence.end());
                 reverse(EdgeSequence.begin(),EdgeSequence.end());
                 vector<int> OptimalSequence;
                 vector<int> SubOptimalSequence;
                 vector<int>SubNodeSequence;
                // printf("hrea 1\n");
                 int ISinNode[N]={0};
                 for(int i=0;i<NodeSquence.size();i++)
                     ISinNode[NodeSquence[i]]=1;
                 for(int i=0;i<edge_num;i++)
                     if(ISinNode[G->incL[i]->head]>0||G->incL[i]->head==s||G->incL[i]->head==t)
                         G->incL[i]->enable=0;
                 int start=0;
                 int end;
                 int d[n_num];
                 //double weight=0.0;
                  //printf("hrea 1\n");
                 for(int i=1;i<NodeSquence.size();i++)
                  {
                    for(int k=0;G->p[NodeSquence[i]][k]>=0;k++)
                      if(G->incL[G->p[NodeSquence[i]][k]]->enable>=0)
                         G->incL[G->p[NodeSquence[i]][k]]->enable=1;
                    if(ISinV[NodeSquence[i]]||NodeSquence[i]==s)
                      { end=i;
                        //printf("\n%d< ----------------------%d\n",NodeSquence[start],NodeSquence[end]);
                        dijkstraWithHeap(G,NodeSquence[end],NodeSquence[start],d,n_num);
                        int f=pre[NodeSquence[start]];
                        SubOptimalSequence.clear();
                        SubNodeSequence.clear();
                        SubOptimalSequence.push_back(peg[NodeSquence[start]]);

                       while(f!=NodeSquence[end]&&f!=-1)
                       {
                        SubNodeSequence.push_back(f);
                       // printf("%d<-",f);
                        SubOptimalSequence.push_back(peg[f]);
                        f=pre[f];
                       }
                       if(f!=-1)
                         { SubNodeSequence.push_back(f);
                            //printf("%d<-",f);
                         }
                        if(f==-1)
                           {
                           for(int j=start;j<=end;j++)
                            {
                             for(int k=0;G->p[NodeSquence[j]][k]>=0;k++)
                                {
                                 G->incL[G->p[NodeSquence[j]][k]]->enable=-2;
                                }
                            }
                           for(int h=start;h<end;h++)
                                {OptimalSequence.push_back(EdgeSequence[h]);
                                 //printf("%d->",EdgeSequence[h]);
                                 weight=weight+G->incL[EdgeSequence[h]]->weight;
                                }
                            }
                        else
                            {//printf("mei gua!\n");
                             for(int j=0;j<SubNodeSequence.size();j++)
                                 {//printf("sail %d\n",SubNodeSequence[j]);
                                  for(int k=0;G->p[SubNodeSequence[j]][k]>=0;k++)
                                       G->incL[G->p[SubNodeSequence[j]][k]]->enable=-2;
                                  }

                             for(int k=0;G->p[NodeSquence[start]][k]>=0;k++)
                                     G->incL[G->p[NodeSquence[start]][k]]->enable=-2;
                             for(int h=0; h<SubOptimalSequence.size();h++)
                                 {OptimalSequence.push_back(SubOptimalSequence[h]);
                                  weight=weight+G->incL[SubOptimalSequence[h]]->weight;
                                  }

                            }
                       int fe=peg[NodeSquence[end]];
                       start=end;
                       //printf("\n");
                      }

                  }
                  vector<int> OptimalNode;
                  for(int i=0;i<OptimalSequence.size();i++)
                     OptimalNode.push_back(G->incL[OptimalSequence[i]]->head);
                  OptimalNode.push_back(t);
                 // printf("in mana2 optimal w is:%f\n",weight);
                      //Check(OptimalSequence);
                    //for(int y=0;y<OptimalSequence.size();y++)
                     // record_result(OptimalSequence[y]);
                   for(int i=0;i<edge_num;i++)
                      G->incL[i]->enable=1;
                   reverse(OptimalSequence.begin(),OptimalSequence.end());
                   reverse(OptimalNode.begin(),OptimalNode.end());
                   return OptimalSequence;
                }




               vector<int> DjiManage(double& weight ,int edge_num,vector<int> EdgeSequence)
                {vector<int> NodeSquence;
                 weight=0;
                 for(int i=0;i<EdgeSequence.size();i++)
                     NodeSquence.push_back(G->incL[EdgeSequence[i]]->head);
                   NodeSquence.push_back(t);
                 //for(int i=0;i<NodeSquence.size();i++)
                   //  printf("->%d",NodeSquence[i]);
                 vector<int> OptimalSequence;
                 vector<int> SubOptimalSequence;
                 vector<int>SubNodeSequence;
                // printf("hrea 1\n");
                 int ISinNode[N]={0};
                 for(int i=0;i<NodeSquence.size();i++)
                     ISinNode[NodeSquence[i]]=1;
                 for(int i=0;i<edge_num;i++)
                     if(ISinNode[G->incL[i]->tail]>0||G->incL[i]->tail==s||G->incL[i]->tail==t)
                         G->incL[i]->enable=0;
                 int start=0;
                 int end;
                 int d[n_num];
                 //double weight=0.0;
                  //printf("hrea 1\n");
                 for(int i=1;i<NodeSquence.size();i++)
                  {
                    for(int k=0;G->f[NodeSquence[i]][k]>=0;k++)
                      if(G->incL[G->f[NodeSquence[i]][k]]->enable>=0)
                         G->incL[G->f[NodeSquence[i]][k]]->enable=1;
                    if(ISinV[NodeSquence[i]]||NodeSquence[i]==t)
                      { end=i;
                        //printf("\n%d ---------------------->%d\n",NodeSquence[start],NodeSquence[end]);
                        dijkstraWithHeap(G,NodeSquence[start],NodeSquence[end],d,n_num);
                        int f=pre[NodeSquence[end]];
                        SubOptimalSequence.clear();
                        SubNodeSequence.clear();
                        SubOptimalSequence.push_back(peg[NodeSquence[end]]);

                       while(f!=NodeSquence[start]&&f!=-1)
                       {
                        SubNodeSequence.push_back(f);
                        //printf("%d<-",f);
                        SubOptimalSequence.push_back(peg[f]);
                        f=pre[f];
                       }
                       if(f!=-1)
                         { SubNodeSequence.push_back(f);
                            //printf("%d<-",f);
                         }
                        if(f==-1)
                           {
                           for(int j=start;j<=end;j++)
                            {
                             for(int k=0;G->f[NodeSquence[j]][k]>=0;k++)
                                {
                                 G->incL[G->f[NodeSquence[j]][k]]->enable=-2;
                                }
                            }
                           for(int h=start;h<end;h++)
                                {OptimalSequence.push_back(EdgeSequence[h]);
                                 //printf("%d->",EdgeSequence[h]);
                                 weight=weight+G->incL[EdgeSequence[h]]->weight;
                                }
                            }
                        else
                            {//printf("mei gua!\n");
                             for(int j=0;j<SubNodeSequence.size();j++)
                                 {//printf("sail %d\n",SubNodeSequence[j]);
                                  for(int k=0;G->f[SubNodeSequence[j]][k]>=0;k++)
                                       G->incL[G->f[SubNodeSequence[j]][k]]->enable=-2;
                                  }

                             for(int k=0;G->f[NodeSquence[end]][k]>=0;k++)
                                     G->incL[G->f[NodeSquence[end]][k]]->enable=-2;
                             for(int h=SubOptimalSequence.size()-1;h>=0;h--)
                                 {OptimalSequence.push_back(SubOptimalSequence[h]);
                                  weight=weight+G->incL[SubOptimalSequence[h]]->weight;
                                  }

                            }
                       int fe=peg[NodeSquence[end]];
                       start=end;
                       //printf("\n");
                      }

                  }
                  vector<int> OptimalNode;
                  for(int i=0;i<OptimalSequence.size();i++)
                     OptimalNode.push_back(G->incL[OptimalSequence[i]]->head);
                  OptimalNode.push_back(t);
                  //for(int i=0;i<OptimalNode.size();i++)
                    //printf("%d->",OptimalNode[i]);
                 //printf("in mana 1 optimal w is:%f\n",weight);
                   //Check(OptimalSequence);
                  for(int i=0;i<edge_num;i++)
                      G->incL[i]->enable=1;
                  return OptimalSequence;


                }

		void parse(Graph* G,char *topo[5000], int edge_num, char *demand){
		    int head,tail,weight,id,j,i,num,mark;
		    int pmark[N]={0};
		    int fmark[N]={0};
		    for(i=0;i<N;i++)
		      for(j=0;j<200;j++)
			   G->p[i][j]=-1;
                     for(i=0;i<N;i++)
		      for(j=0;j<200;j++)
			   G->f[i][j]=-1;
		    while(*demand != ',')
			num = num*10 + *(demand++) - '0';
		    t = num;
                    //printf("num is :%d",num);
		    num = 0;
		    demand++;
		    while(*demand != ',')
			num = num*10 + *(demand++) - '0';
		    s = num;
		    demand++;
		    sscanf(demand,"%d",&i);
		    V.push_back(i);
		    ISinV[i]=1;
		    demand++;
			 if(i>9)
			     demand++;
			 if(i>99)
			     demand++;
		    while(sscanf(demand,"|%d",&i)==1)
			{ V.push_back(i);
			  ISinV[i]=1;
			 demand=demand+2;
			 if(i>9)
			     demand++;
			 if(i>99)
			     demand++;
			}
		    for(j=0;j<edge_num;j++)
		    { sscanf(topo[j],"%d,%d,%d,%d",&id,&tail,&head,&weight);
		      if(head>n_num)
			  num=head;
		      if(tail>n_num)
			  n_num=tail;
		      mark=0;
		      for(i=0;i<V.size();i++)
			  if(V[i]==tail)
			      mark=1;
		      Edge* e=(Edge*)malloc(sizeof(Edge));
		      if(mark==1)
			Graph_addEdge(e,head,tail,(double)weight,mark,weight,1,InitPheromone+CTA);
		      else
			 if(tail==t)
			     Graph_addEdge(e,head,tail,(double)(weight),mark,weight,1,InitPheromone);
			  else
			      Graph_addEdge(e,head,tail,(double)weight,mark,weight,1,InitPheromone);
		      G->incL[j]=e;
		      //printf("pmark :%d\n",pmark[head]);
		      G->p[head][pmark[head]++]=j;
		        G->f[tail][fmark[tail]++]=j;
		    }
		  // printf("t is :%d\n",t);
		   n_num++;

		}
		int ANT_NUM=10000;
		int Times=10000;
		void SolveByAnt(char *topo[5000], int edge_num, char *demand)
		{   srand(time(0));
		    time_t t_s=clock();
		    vector<vector<int> >Sequences;
		    Ant* Colony[ANT_NUM];
		    vector<int> OptimalSequence;
                    vector<int> EdgeSequence;
		    //printf("here !1\n");
		    for(int i=0;i<ANT_NUM;i++)
		       Colony[i]=new Ant(i);
		    //int Count =10000;
		     //printf("here !2\n");
		    for(int j=0;j<Times;j++)
			{//printf(".................................................\n");
			 time_t t_t = clock();
			 if((double)(t_t - t_s)>9870000.0)
			     break;
			 int update=0;
			 int max=0;
			 //Count=Count+100;
			  //printf("here !3\n");
			 for(int k=0;k<ANT_NUM;k++)
			     { //printf("here !31\n");
			      Colony[k]->move();
			     // printf("here !32\n");
			      //printf("shit sadadd!\n");
			      if(Colony[k]->VistedV.size()>max)
				   max=Colony[k]->VistedV.size();
			      if(Colony[k]->Findt>0)
				  update=1;
			      //printf("shit!\n");
			     }
			 //printf("max is :%d\n",max);
			 // printf("here !4\n");
			 if(update>0)
			     for(int k=0;k<edge_num;k++)
				  G->incL[k]->Pheromone= Adapter(G->incL[k]->Pheromone*GAMA);
			 for(int k=0;k<ANT_NUM;k++)
			    {if(ISinV[Colony[k]->CurrentPosition])
				 {for(int l=0;l<Colony[k]->Path.size();l++)
				    G->incL[Colony[k]->Path[l]]->Pheromone=Adapter(G->incL[Colony[k]->Path[l]]->Pheromone+PI);
				 }

			     if(Colony[k]->Findt>0)//&&Colony[k]->VistedV.size()>=V.size())
				 {
                                 for(int l=0;l<Colony[k]->Path.size();l++)
				    G->incL[Colony[k]->Path[l]]->Pheromone=Adapter(G->incL[Colony[k]->Path[l]]->Pheromone+NAMENDA*pow((double)Colony[k]->VistedV.size()/(double)V.size(),ROU));

				  if(Colony[k]->VistedV.size()>=V.size()&&OptimalWeight>Colony[k]->PathLen)
				      {//printf("haha!\n");

                                       //printf("optimal !!,weight is:%.2f\n",Colony[k]->PathLen);
				        OptimalWeight=Colony[k]->PathLen;
				       OptimalSequence.clear();
                                       EdgeSequence.clear();
				       for(int h=0;h<Colony[k]->NodeSequence.size();h++)
				          OptimalSequence.push_back(Colony[k]->	NodeSequence[h]);
				       for(int h=0;h<Colony[k]->Path.size();h++)
				          EdgeSequence.push_back(Colony[k]->Path[h]);
				       Sequences.push_back(EdgeSequence);

				      }

				 }
			     }
			   // printf("here !5\n");
			 }
                      //Check(OptimalSequence);
                      //for(int i=0;i<OptimalSequence.size();i++)
                        //  printf("%d->",OptimalSequence[i]);
                      //printf("......................................................\n");

                     vector<vector<int> > Opti;
                     double minw=INF;
                     int Optimali;
		     if(Sequences.size()>0)
                    {
                       for(int i=0;i<Sequences.size();i++)
                       {  double w1=0;
                          Opti.push_back(DjiManage(w1,edge_num,DjiManage2(w1,edge_num,Sequences[i])));
                          if(w1<minw)
                          {Optimali=Opti.size()-1;
                           minw=w1;
                          }
                          Opti.push_back(DjiManage2(w1,edge_num,DjiManage(w1,edge_num,Sequences[i])));
                           if(w1<minw)
                          {Optimali=Opti.size()-1;
                           minw=w1;
                          }
                          Opti.push_back(DjiDirect(w1,edge_num,Sequences[i]));
                          if(w1<minw)
                          {Optimali=Opti.size()-1;
                           minw=w1;
                          }
                          Opti.push_back(DjiDirect2(w1,edge_num,Sequences[i]));
                          if(w1<minw)
                          {Optimali=Opti.size()-1;
                           minw=w1;
                          }

                       }

                    double w=0;
                    reverse(Opti[Optimali].begin(),Opti[Optimali].end());
                    for(int i=0;i<Opti[Optimali].size();i++)
                      record_result(Opti[Optimali][i]);
                   // Check(Opti[Optimali]);
                     //printf("opti:%d,minw:%f,actw:%f\n",Optimali,minw,w);
                    }


		}
	  }//namespaceend;
vector<int> solv(char *topo[5000], int edge_num, char *demand)
{
    int n = 0;
    int a, b, c, d;
    ANT::parse(ANT::G,topo,edge_num,demand);
    for(int i = 0; i < edge_num; i++)
    {
        sscanf(topo[i], "%d,%d,%d,%d", &a, &b, &c, &d);
        n = max(n, max(b, c));
    }
    if(n >= 550)
    {
        ANT::SolveByAnt(topo,edge_num,demand);
        vector<int> ans;
        return ans;
    }
    bool flag = false;
    pair<int, vector<int>> ans1 = BAOLI::solvx(topo, edge_num, demand, flag);
    if(flag) return ans1.second;
    pair<int, vector<int>> ans2 = pAth::solvx(topo, edge_num, demand, flag);
    if(flag) return ans2.second;

    pair<int, vector<int>> ans3 = PATH::solvx(topo, edge_num, demand);

    if(ans1.first < ans2.first && ans1.first < ans3.first)
        return ans1.second;
    if(ans2.first < ans1.first && ans2.first < ans3.first)
        return ans2.second;

    else
        return ans3.second;
   }

void search_route(char *topo[5000], int edge_num, char *demand)
{
    vector<int> ans = solv(topo,edge_num,demand);
    for(int i = 0; i < ans.size(); i++)
        record_result(ans[i]);
}
