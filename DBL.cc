#include "DBL.h"
#include <omp.h>
#include <climits>
#include <xmmintrin.h>
#include <algorithm>
#include <cstdint>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <bitset>
#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <random>

using namespace std;
#define rep(i, n) for (int i = 0; i < (int)(n); i++)
//#define test


double get_wall_time() {
	struct timeval time;
	if (gettimeofday(&time, NULL)) {
		return 0;
	}
	return (double) time.tv_sec + (double) time.tv_usec * .000001;
}

namespace {

int get_max_threads() {
#ifdef _OPENMP
	return omp_get_max_threads();
#else
	return 1;
#endif
}

int get_thread_id() {
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}

template<typename T>
struct parallel_vector {
	parallel_vector(size_t size_limit) :
			v(get_max_threads(), vector<T>(size_limit)), n(get_max_threads(), 0) {
	}

	void push_back(const T &x) {
		int id = get_thread_id();
		v[id][n[id]++] = x;
	}

	void clear() {
		rep (i, get_max_threads())
			n[i] = 0;
	}

	vector<vector<T> > v;
	vector<size_t> n;
};
}


void DBL::insert_edges(const char *filename,
		int num) {
	std::ifstream ifs(filename);
	_CHECK(!ifs.bad());
	insert_edges(ifs, num);
}

void DBL::insert_edges(
		std::istream &ifs,
		int num)
{
	printf("Reading Edges......\n");
	std::vector<std::tuple<int, int> > es;
	for (int v, w, i = 0; i < num; i++) {
		ifs >> v >> w;
		es.emplace_back(v, w);
	}
	insert_edges(es);
}



void DBL::insert_edges(
		vector<tuple<int, int>> &es) {

	printf("Insertion Start\n");
	int size = (int) es.size();
	double time1 = get_wall_time();
	for (int i = 0; i < size; i++)
       	{
		//ignore error input
		if(get<0>(es[i]) >= (int)adj.size() || get<1>(es[i]) >= (int)adj.size())
		{
			continue;
		}
		
		adj[get < 0 > (es[i])].insert(get < 1 > (es[i]));
		rev_adj[get < 1 > (es[i])].insert(get < 0 > (es[i]));
                
		if ((dl_out[get < 0 > (es[i])] & dl_in[get < 1 > (es[i])]).any()) {
			continue;
		}

		vector<int> curr_que;
		vector<int> next_que;

		curr_que.push_back(get < 1 > (es[i]));
		dl_in[get < 1 > (es[i])] |= dl_in[get < 0 > (es[i])];
		bl_in[get < 1 > (es[i])] |= bl_in[get < 0 > (es[i])];

		while (!curr_que.empty()) 
		{
			for (int j = 0; j < (int) curr_que.size(); j++) 
			{
				for (auto x = adj[curr_que[j]].begin(); *x != INT_MAX; x++) {
					if ((dl_in[get < 0 > (es[i])] | dl_in[*x]) != dl_in[*x]
					|| (bl_in[get < 0 > (es[i])] | bl_in[*x]) != bl_in[*x])
					{
						dl_in[*x] |= dl_in[get < 0 > (es[i])];
						bl_in[*x] |= bl_in[get < 0 > (es[i])];
						next_que.push_back(*x);
					}
				}
			}

			curr_que.assign(next_que.begin(), next_que.end());
			next_que.clear();
		}

		curr_que.clear();
		next_que.clear();

		curr_que.push_back(get < 0 > (es[i]));
		dl_out[get < 0 > (es[i])] |= dl_out[get < 1 > (es[i])];
		bl_out[get < 0 > (es[i])] |= bl_out[get < 1 > (es[i])];
		
		while (!curr_que.empty()) {
			for (int j = 0; j < (int) curr_que.size(); j++) {
				for (auto x = rev_adj[curr_que[j]].begin(); *x != INT_MAX; x++) {
					if ((dl_out[get < 1 > (es[i])] | dl_out[*x]) != dl_out[*x]
					|| (bl_out[get < 1 > (es[i])] | bl_out[*x]) != bl_out[*x])
					{
						dl_out[*x] |= dl_out[get < 1 > (es[i])];
						bl_out[*x] |= bl_out[get < 1 > (es[i])];
						next_que.push_back(*x);
					}
				}
			}

			curr_que.assign(next_que.begin(), next_que.end());
			next_que.clear();
	 	}
	}

	double time3 = get_wall_time();
	printf("The Insertion Time = %f\n\n\n", time3 - time1);

}



void DBL::test_query(vector<int> &u,   //Using BFS to validate the correctness of the query
		vector<int> &v, int sz, bool *test_result, bool *result) {

	double time1 = get_wall_time();
	const int threadnum = 1;

#pragma omp parallel for num_threads(threadnum),schedule(static, 1)
	for (int i = 0; i < sz; i++) {

		vector<vector<bool>> visited(2, vector<bool>(V,false));
		queue<int> q[2]; 
		q[0].push(u[i]);
		q[1].push(v[i]);
		visited[0][u[i]] = true;
	       	visited[1][v[i]] = true;

		for (int t = 0; q[t].size() > 0; t = ((int)q[0].size() <= (int)q[1].size() ? 0 : 1)) {
			for (int qsz = q[t].size(), qi = 0; qi < qsz; ++qi) {
				int f = q[t].front(); 
				q[t].pop();
				int pfsm = t == 0 ? prefixsum[f] : rev_prefixsum[f];
				if (visited[!t][f]) 
				{
					test_result[i] = true;
					break;		
				}
				
				int w = (t == 0? csr[pfsm] : rev_csr[pfsm]); 
				
				while(w != INT_MAX)
				{
					if(t==0)
					{
						if (!visited[t][w]) 
						{
							visited[t][w] = true;
							q[t].push(w);
						}
					}
					
					if(t==1)
					{
						if (!visited[t][w]) 
                  		                {
					            visited[t][w] = true;
						    q[t].push(w);
					        }
					
					}
					pfsm++;
					w = (t == 0? csr[pfsm] : rev_csr[pfsm]);
				}
			}
			if(test_result[i]) break;

		}
		
		if (result[i] != test_result[i])
			printf("error %d:%d,%d\n", i, u[i], v[i]);

	}

	double time2 = get_wall_time();

	printf("Test query Time:%f\n", time2 - time1);

}

void DBL::query_label(
		vector<int> &u,
		vector<int> &v,
		int sz,
		bool *sure,
		bool *result)
{

#pragma omp parallel for num_threads(query_threadNum),schedule(guided, 1)
	for (int i = 0; i < sz; i++) {
		if ((dl_out[u[i]] & dl_in[v[i]]).any()) {
			sure[i] = true;
			result[i] = true;
			continue;
		} 
		if ( (bl_in[v[i]] | bl_in[u[i]]) != bl_in[v[i]]
		  || (bl_out[v[i]] | bl_out[u[i]]) != bl_out[u[i]]) {
			sure[i] = true;
			continue;
		} 
		if ((dl_out[v[i]] & dl_in[u[i]]).any()) {
			sure[i] = true;
			continue;
		}
		if ((dl_out[u[i]] & dl_in[u[i]]).any()
		|| (dl_out[v[i]] & dl_in[v[i]]).any()) {
			sure[i] = true;
			continue;
		}
		if (u[i] == v[i]) {
			sure[i] = true;
			result[i] = true;
			continue;
		}

	}

}

void DBL::query(vector<int> &u,
		vector<int> &v, int sz, bool *result) {
	bool *sure = new bool[sz];
	fill_n(sure, sz, false);

	double time1 = get_wall_time();
	query_label(u, v, sz, sure, result);
	double time2 = get_wall_time();

	printf("Label Query Time:%f\n", time2 - time1);
	double time3 = get_wall_time();

#pragma omp parallel for num_threads(query_threadNum),schedule(guided, 1)
	for (int i = 0; i < sz; i++) {
		if (sure[i] == false && result[i] == false) {
	
			vector<int> curr_que;
			vector<int> next_que;
			vector<bool> visited(V, false);
			curr_que.push_back(u[i]);
			while (!curr_que.empty()) {
				for (int j = 0; j < (int) curr_que.size(); j++) {
					int d = prefixsum[curr_que[j]];
					while (csr[d] != INT_MAX) {
						if (csr[d] == v[i]) {
							result[i] = true;
							break;
						}

						if (visited[csr[d]] == false
						&& (dl_out[u[i]] & dl_in[csr[d]]).none()
		 				&& ((bl_in[v[i]] | bl_in[csr[d]]) == bl_in[v[i]]
					    && (bl_out[v[i]] | bl_out[csr[d]]) == bl_out[csr[d]])
					    )
						{
							next_que.push_back(csr[d]);
							visited[csr[d]] = true;
					//		countnode++;
						}
						d++;
					}
					if (result[i] == true) {
						break;
					}

				}
				if (result[i] == true) {
					break;
				}
				curr_que.assign(next_que.begin(), next_que.end());
				next_que.clear();
			}
		}
	}

	double time4 = get_wall_time();
	printf("BFS Query Time:%f\n", time4 - time3);
	printf("Total Query Time:%f\n", time4 - time1);
	int trueresult = count(result,result+sz,true);
  	printf("Number of Reachable Query = %d\n", trueresult);
	printf("\n\n\n");

#ifdef test
	printf("Test Correctness\n");

	bool *test_result = new bool[sz];
	fill_n(test_result, sz, false);
	test_query(u, v, sz, test_result, result);
	
	int trueresult2 = count(test_result, test_result + sz, true);
	printf("Number of test reachable answer = %d\n", trueresult2);
	printf("Test finished\n");
#endif

}

void split(const string& s,vector<int>& sv,const char flag) {
    sv.clear();
    istringstream iss(s);
    string temp;

    while (getline(iss, temp, flag)) {
            sv.push_back(stoi(temp));
        }
    return;
}

void DBL::construct_index(char *filename1)
{
	std::ifstream ifs1(filename1);
	construct_index(ifs1);

}

void DBL::construct_index(istream &ifs1)
{
	std::vector<std::tuple<int, int> > es;

	printf("Reading Graph\n");
	for (int v, w; ifs1 >> v >> w;) {
		es.emplace_back(v, w);
	}

	V = 0;
	E = 0;
	for (const auto &e : es) {
		V = max( { V, get < 0 > (e) + 1, get < 1 > (e) + 1 });
		E++;
	}

	printf("Number of Nodes:%d\n", V);
	printf("Number of Edges:%d\n", E);
	printf("Constructing Adjacent List......\n");

	adj.assign(V, set<int>());
	rev_adj.assign(V, set<int>());

	for (const auto &e : es) {
		adj[get < 0 > (e)].insert(get < 1 > (e));
		rev_adj[get < 1 > (e)].insert(get < 0 > (e));
	}

	int adj_size = 0;
	int *degree = new int[V];
	int *rev_degree = new int[V];
	prefixsum.resize(V);
	rev_prefixsum.resize(V);

	for (int i = 0; i < V; i++) {
		adj[i].insert(INT_MAX);
		rev_adj[i].insert(INT_MAX);
		degree[i] = (int) adj[i].size();
		rev_degree[i] = (int) rev_adj[i].size();
		adj_size += degree[i];
	}
	csr_size = adj_size;
	prefixsum[0] = 0;
	rev_prefixsum[0] = 0;
	for (int i = 1; i < V; i++) {
		prefixsum[i] = prefixsum[i - 1] + degree[i - 1];
		rev_prefixsum[i] = rev_prefixsum[i - 1] + rev_degree[i - 1];
	}

	csr.resize(adj_size);
	rev_csr.resize(adj_size);

	for (int i = 0; i < V; i++) {
		copy(adj[i].begin(), adj[i].end(), csr.begin() + prefixsum[i]);
		copy(rev_adj[i].begin(), rev_adj[i].end(),
				rev_csr.begin() + rev_prefixsum[i]);
	}

	int capacity=0;
	for(int i=0; i<(int)adj.size(); i++)
	{
		capacity+=(int)adj[i].size();
	}

	printf("Construct labels\n");
	double start_time1 = get_wall_time();
	construct_index_dl(csr, prefixsum, rev_csr, rev_prefixsum);
	construct_index_bl(csr, prefixsum, rev_csr, rev_prefixsum);
	double end_time1 = get_wall_time();
	printf("Total Construction Time:%f\n", end_time1 - start_time1);

}




void DBL::construct_index_dl(
		vector<int> &csr, vector<int> &prefixsum, vector<int> &rev_csr,
		vector<int> &rev_prefixsum) {
	double start_time1 = get_wall_time();

	printf("Constructing DL Label......\n");
	dl_in.clear();
	dl_in.resize(V);
	dl_out.clear();
	dl_out.resize(V);
	dl_hash.resize(V);
	fill_n(&dl_hash[0], V,-1);

	get_root_order(ord, rev_ord);
//	get_root_order_max(ord, rev_ord);
//	get_root_order_min(ord, rev_ord);
//	get_root_order_plus(ord, rev_ord);
//	get_root_order_betweeness(ord, rev_ord)



	parallel_vector<int> p_nxt_que_in(V), p_nxt_que_out(V);

	rep (source_i, (int) dl_in[0].size())
	{
		int s = ord[source_i];
		dl_hash[s]=source_i;
		vector<int> crr_que_in, crr_que_out;
		crr_que_in.push_back(s);
		crr_que_out.push_back(s);
		set<int> crr_que_in_set, crr_que_out_set;


		bool *visited_in = new bool[V], *visited_out = new bool[V];
		fill_n(visited_in, V, false);
		fill_n(visited_out, V, false);

		visited_in[s] = true;
		visited_out[s] = true;

		for (int d = 0; !crr_que_out.empty(); ++d) {
			int crr_que_size_out = crr_que_out.size();

#pragma omp parallel for num_threads(construct_threadNum),schedule(static, 1)
			rep (que_i, crr_que_size_out)
			{
				int v = crr_que_out[que_i];

				auto s1 = dl_out[s];
				auto s2 = dl_in[v];

				if ((s1 & s2).none())
				{
					dl_in[v].set(source_i, 1);
					for (int i = prefixsum[v]; csr[i] != INT_MAX; i++) {
						if (visited_out[csr[i]] == false) {
							visited_out[csr[i]] = true;
							p_nxt_que_out.push_back(csr[i]);
						}
					}

				}
			}

#pragma omp flush
			crr_que_out.clear();
			rep (i, get_max_threads())
			{
				rep (j, p_nxt_que_out.n[i])
				{
					int v = p_nxt_que_out.v[i][j];
					crr_que_out_set.insert(v);
				}
			}
			p_nxt_que_out.clear();
			crr_que_out.assign(crr_que_out_set.begin(), crr_que_out_set.end());
			crr_que_out_set.clear();
		}

		for (int d = 0; !crr_que_in.empty(); ++d) {
			int crr_que_size_in = crr_que_in.size();

#pragma omp parallel for num_threads(construct_threadNum),schedule(static, 1)
			rep (que_i, crr_que_size_in)
			{
				int v = crr_que_in[que_i];

				auto s1 = dl_in[s];
				auto s2 = dl_out[v];

				if ((s1 & s2).none())
				{
					dl_out[v].set(source_i, 1);
					for (int i = rev_prefixsum[v]; rev_csr[i] != INT_MAX; i++) {
						if (visited_in[rev_csr[i]] == false) {
							visited_in[rev_csr[i]] = true;
							p_nxt_que_in.push_back(rev_csr[i]);
						}
					}
				}

			}

#pragma omp flush

			crr_que_in.clear();

			rep (i, get_max_threads())
			{
				rep (j, p_nxt_que_in.n[i])
				{
					int v = p_nxt_que_in.v[i][j];
					crr_que_in_set.insert(v);

				}
			}
			p_nxt_que_in.clear();
			crr_que_in.assign(crr_que_in_set.begin(), crr_que_in_set.end());
			crr_que_in_set.clear();

		}

	}

	double end_time1 = get_wall_time();
	printf("DL Label Construction Time:%f\n", end_time1 - start_time1);

}

void DBL::construct_index_bl(vector<int> &csr,
		vector<int> &prefixsum, vector<int> &rev_csr,
		vector<int> &rev_prefixsum) {

	printf("Constructing Bl Label......\n");
	double start_time1 = get_wall_time();

	parallel_vector<int> p_nxt_que_in(V), p_nxt_que_out(V);
	bl_in.resize(V);
	bl_out.resize(V);

	vector<int> startnode, endnode;
	for (int i = 0; i < V; i++) {
		int adj_size = (int) adj[i].size();
		int rev_adj_size = (int) rev_adj[i].size();
		if (adj_size != 1 && rev_adj_size == 1) {
			startnode.push_back(i);
		}
		if (adj_size == 1 && rev_adj_size != 1) {
			endnode.push_back(i);
		}
	}
	
	shuffle(startnode.begin(), startnode.end(), default_random_engine());
	shuffle(endnode.begin(), endnode.end(), default_random_engine());

	int startsize = (int) startnode.size() / (int) bl_in[0].size() + 1;
	int endsize = (int) endnode.size() / (int) bl_out[0].size() + 1;

	vector<int> crr_que_in, crr_que_out;
	bool *visited_in = new bool[V], *visited_out = new bool[V];

	rep (source_i, (int)bl_in[0].size())
	{
		int hash = source_i;

		fill_n(visited_out, V, false);
		for (int i = source_i * startsize;
				i < (source_i + 1) * startsize && i < (int) startnode.size();
				i++) {
			bl_in[startnode[i]].set(hash, 1);
			crr_que_out.push_back(startnode[i]);
			visited_out[startnode[i]] = true;
		}

		for (int d = 0; !crr_que_out.empty(); ++d) {
			int crr_que_size_out = crr_que_out.size();

#pragma omp parallel for num_threads(construct_threadNum),schedule(guided, 1)
			rep (que_i, crr_que_size_out)
			{
				int v = crr_que_out[que_i];

				bl_in[v].set(hash, 1);
				for (int i = prefixsum[v]; csr[i] != INT_MAX; i++) {
					if (visited_out[csr[i]] == false) {
						visited_out[csr[i]] = true;
						p_nxt_que_out.push_back(csr[i]);
					}
				}

			}
#pragma omp flush

			crr_que_out.clear();
			rep (i, get_max_threads())
			{
				rep (j, p_nxt_que_out.n[i])
				{
					int v = p_nxt_que_out.v[i][j];
					crr_que_out.push_back(v);
				}
			}
			p_nxt_que_out.clear();

		}

	}

	rep (source_i, (int)bl_out[0].size())
	{

		int hash = source_i;

		fill_n(visited_in, V, false);

		for (int i = source_i * endsize;
				i < (source_i + 1) * endsize && i < (int) endnode.size(); i++) {
			bl_out[endnode[i]].set(hash, 1);
			crr_que_in.push_back(endnode[i]);
			visited_in[endnode[i]] = true;
		}

		for (int d = 0; !crr_que_in.empty(); ++d) {
			int crr_que_size_in = crr_que_in.size();

#pragma omp parallel for num_threads(construct_threadNum),schedule(guided, 1)
			rep (que_i, crr_que_size_in)
			{
				int v = crr_que_in[que_i];

				bl_out[v].set(hash, 1);
				for (int i = rev_prefixsum[v]; rev_csr[i] != INT_MAX; i++) {
					if (visited_in[rev_csr[i]] == false) {
						visited_in[rev_csr[i]] = true;
						p_nxt_que_in.push_back(rev_csr[i]);
					}
				}

			}
#pragma omp flush

			crr_que_in.clear();
			rep (i, get_max_threads())
			{
				rep (j, p_nxt_que_in.n[i])
				{
					int v = p_nxt_que_in.v[i][j];
					crr_que_in.push_back(v);
				}
			}
			p_nxt_que_in.clear();

		}

	}

	double end_time1 = get_wall_time();
	printf("BL Label construction time:%f\n", end_time1 - start_time1);

}




void DBL::get_root_order(vector<int> &ord,
		vector<int> &rev_ord) {
	vector<pair<int, int> > deg(V), deg_out(V);
	vector<pair<double, int>> deg_in(V);
	rep (v, V)
		deg_out[v] = make_pair(adj[v].size(), v);
	rep (v, V)
		deg_in[v] = make_pair(rev_adj[v].size(), v);
	rep (v, V)
		deg[v] = make_pair(deg_in[v].first/1000*deg_out[v].first, v);

	sort(deg.begin(), deg.end());
	reverse(deg.begin(), deg.end());

	ord.resize(V);
	rev_ord.resize(V);
	rep (i, V)
		ord[i] = deg[i].second;
	rep (i, V)
		rev_ord[ord[i]] = i;
}

void DBL::get_root_order_max(vector<int> &ord, vector<int> &rev_ord)
{
	vector<pair<int,int>> deg(V), deg_out(V), deg_in(V);
	rep(v, V)
		deg_out[v] = make_pair(adj[v].size(), v);
	rep(v, V)
		deg_in[v] = make_pair(rev_adj[v].size(), v);
	rep(v, V)
		deg[v] = make_pair(deg_in[v].first, v);
	sort(deg.begin(), deg.end());
	reverse(deg.begin(), deg.end());
	
	ord.resize(V);
	rev_ord.resize(V);
	rep (i, V)
		ord[i] = deg[i].second;
	rep (i, V)
		rev_ord[ord[i]] =i;

}


void DBL::get_root_order_min(vector<int> &ord, vector<int> &rev_ord)
{
	vector<pair<int,int>> deg(V), deg_out(V), deg_in(V);
	rep(v, V)
		deg_out[v] = make_pair(adj[v].size(), v);
	rep(v, V)
		deg_in[v] = make_pair(rev_adj[v].size(), v);
	rep(v, V)
		deg[v] = make_pair(deg_out[v].first, v);
	sort(deg.begin(), deg.end());
	reverse(deg.begin(), deg.end());
	
	ord.resize(V);
	rev_ord.resize(V);
	rep (i, V)
		ord[i] = deg[i].second;
	rep (i, V)
		rev_ord[ord[i]] =i;

}




void DBL::get_root_order_plus(vector<int> &ord, vector<int> &rev_ord)
{
	vector<pair<int,int>> deg(V), deg_out(V), deg_in(V);
	rep(v, V)
		deg_out[v] = make_pair(adj[v].size(), v);
	rep(v, V)
		deg_in[v] = make_pair(rev_adj[v].size(), v);
	rep(v, V)
		deg[v] = make_pair(deg_in[v].first+deg_out[v].first, v);
	sort(deg.begin(), deg.end());
	reverse(deg.begin(), deg.end());
	
	ord.resize(V);
	rev_ord.resize(V);
	rep (i, V)
		ord[i] = deg[i].second;
	rep (i, V)
		rev_ord[ord[i]] =i;

}



