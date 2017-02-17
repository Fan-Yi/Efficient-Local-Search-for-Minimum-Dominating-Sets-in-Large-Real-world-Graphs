#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

//#define windows_mode
#define linux_mode

#ifdef linux_mode
#include <sys/times.h>
#endif

#ifdef windows_mode
#include <windows.h>
#endif

#include <cstdlib>
#include  <stdio.h>

#include "graph.h"
#include "myBijection.h"
//#include "myBuckets.h"
#include "ansUpdate.h"
#include "hugeInt.h"
#include "operandSets.h"
#include "dsBuckets.h"
#include "dominationHash.h"
#include "constants.h"

//#define strategy_analysis_mode

//#define only_test_init_mode

#define ans_update_op_mode
#define age_init_i_mode

//#define edge_time_stamp_mode
#define individual_analysis_on_init_sls_mode

#define greedy_construction_mode
//#define domination_hash_mode

//#define nvcc_mode
#define score_changed_mode

//#define probablistic_remove_mode
//#define double_random_remove_mode

//#define consecutive_collision_mode

//#define debug_mode

class StochasticLocalSearch : private Graph{
private:
	int* in_dominating_set;
	int dominating_set_size;
	int* connect_dominating_set_degree;
	//ScoreBasedPartition *ptr_to_partitioned_dominating_set;
	DS_Partition *ptr_to_partitioned_dominating_set;
	//Bijection* ptr_to_uncov_vertices;
	ConfchangedScoreAgePickSet* ptr_to_uncov_vertices;
#ifdef ans_update_op_mode
	AnsChangeSet* ptr_to_moved_v;
#endif
	int *score;
	HugeInt *time_stamp;

#ifdef edge_time_stamp_mode
	HugeInt* edge_time_stamp;
#endif
	int *confChange;
#ifdef probablistic_remove_mode
	int *multi_vertex_set;
#endif

	HugeInt step;
	int time_limit;
#ifdef linux_mode
	tms start, finish;
	int start_time;
#endif
#ifdef windows_mode
	unsigned long start_time;
	unsigned long best_cmp_time;
#endif
	int *best_in_dominating_set;
	int best_dominating_set_size;
#ifdef linux_mode
	double best_cmp_time;
#endif
	HugeInt best_solve_step;

#ifdef individual_analysis_on_init_sls_mode
	#ifdef linux_mode
		double init_time;
		double sls_time;
	#endif
	#ifdef windows_mode
		unsigned long init_time;
		unsigned long sls_time;
	#endif
#endif

#ifdef domination_hash_mode
	DominationHash *ptr_to_hashed_domination;
	int consecutive_collision_num;
#endif

	int check_solution()
	{
		int v;
		int verified_dominating_set_size = 0;
		for(v = 1; v <= v_num; v++)
		{
			if(best_in_dominating_set[v])
			{
				verified_dominating_set_size++;
			}
		}
		if(verified_dominating_set_size != best_dominating_set_size)
		{
			cout << "verified_dominating_set_size: " << verified_dominating_set_size << endl;
			cout << "best_dominating_set_size: " << best_dominating_set_size << endl;
			cout << "the dominating set size is computed incorrectly" << endl;
			return 0;
		}
		for(v = 1; v <= v_num; v++)
		{
			if(best_in_dominating_set[v]) continue;
			int flag = 0;
			for(int i = 0; i < vertices[v].get_degree(); i++)
			{
				int *nbs = vertices[v].get_neighbors();
				if(best_in_dominating_set[nbs[i]])
				{
					flag = 1;
					break;
				}
			}
			if(flag == 0)
			{
				cout << "vertex " << v << " is uncovered" << endl;
				return 0;
			}
		}
		return 1;
	}

#ifndef greedy_construction_mode
void init_solution() //random construction mode
{
#ifdef debug_mode
cout << "rand construction" << endl;
#endif
	for(int v = 1; v <= v_num; v++)
	{
		ptr_to_uncov_vertices->insert_element(v);
	}
	while(ptr_to_uncov_vertices->size())
	{
		int v = ptr_to_uncov_vertices->rand_element();
		add(v);
//show_state();
	}
	ptr_to_partitioned_dominating_set->filter_out_sole_vertices(vertices, v_num);
#ifdef debug_mode
cout << "after filtering, " << ptr_to_partitioned_dominating_set->get_remaining_vertex_num() << " elements exist" << endl;
getchar();
#endif
	update_best_dominating_set();
	//ptr_to_moved_v->efficiently_update_best_dominating_set_vertices();
	best_dominating_set_size = dominating_set_size;
//cout << "after construction, o " << best_dominating_set_size << endl;
#ifdef linux_mode
	times(&finish);
	best_cmp_time = double(finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
#endif
#ifdef windows_mode
	best_cmp_time = timeGetTime() - start_time;
#endif
	best_solve_step = step;

/*
	if(optimal_solution_found)
	{
		if(check_solution() == 1)
			{
				cout << "o " << best_dominating_set_size << endl;
				cout << "c searchSteps " << best_solve_step <<endl;
				cout << "c solveTime " << best_cmp_time << endl;
			}
			else
			{
				cout << "the solution is wrong." << endl;
			}
			exit(0);
	}
*/
#ifdef only_test_init_mode
	//else
	{
		if(check_solution() == 1)
			{
				cout << "o " << best_dominating_set_size << endl;
				cout << "c searchSteps " << best_solve_step <<endl;
				cout << "c solveTime " << best_cmp_time << endl;
			}
			else
			{
				cout << "the solution is wrong." << endl;
			}
			exit(0);
	}
#endif

#ifdef individual_analysis_on_init_sls_mode
	#ifdef linux_mode
		init_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime) / sysconf(_SC_CLK_TCK);
		init_time = round(init_time * 100)  / 100.0;
	#endif
	#ifdef windows_mode
		init_time = timeGetTime() - start_time;
	#endif
#endif
}

#else

void init_solution()
{
	int* v_ascending_degree_array = ptr_to_partitioned_dominating_set->get_v_ascending_degree_array();
	int* ptr_to_v_gain = ptr_to_partitioned_dominating_set->get_ptr_to_v_gain();
	int* ptr_to_v_loss = ptr_to_partitioned_dominating_set->get_ptr_to_v_loss();
////////////////////////////////
	int k = get_max_degree() + 1;
#ifdef debug_mode
cout << "greedy construction" << endl;
cout << "k: " << k << endl;
#endif
	for(int v = 1; v <= v_num; v++)
	{
		ptr_to_uncov_vertices->insert_element(v);
	}
	while(ptr_to_uncov_vertices->size())
	{
		int v;
		while(ptr_to_v_loss[k+1] == ptr_to_v_gain[k]) k--;
//cout << "k: " << k << endl;
		v = v_ascending_degree_array[ptr_to_v_gain[k] + rand() % (ptr_to_v_loss[k+1] - ptr_to_v_gain[k])];
		add(v);
//if(!ptr_to_partitioned_dominating_set->checkPartition(v_num, score)) exit(1);
		//uncov_vertex_num -= score[v];

//show_state();
	}
	ptr_to_partitioned_dominating_set->filter_out_sole_vertices(vertices, v_num);
	update_best_dominating_set();
	//ptr_to_moved_v->efficiently_update_best_dominating_set_vertices();
	best_dominating_set_size = dominating_set_size;
//cout << "after construction, o " << best_dominating_set_size << endl;
#ifdef linux_mode
	times(&finish);
	best_cmp_time = double(finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
#endif
#ifdef windows_mode
	best_cmp_time = timeGetTime() - start_time;
#endif
	best_solve_step = step;

#ifdef only_test_init_mode
	{
		if(check_solution() == 1)
			{
				cout << "o " << best_dominating_set_size << endl;
				cout << "c searchSteps " << best_solve_step << endl;
				cout << "c solveTime " << best_cmp_time << endl;
				//cout << best_dominating_set_size << '\t' << best_cmp_time << endl;
			}
			else
			{
				cout << "the solution is wrong." << endl;
			}
			exit(0);
	}
#endif

#ifdef individual_analysis_on_init_sls_mode
	#ifdef linux_mode
		init_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime) / sysconf(_SC_CLK_TCK);
		init_time = round(init_time * 100)  / 100.0;
	#endif
	#ifdef windows_mode
		init_time = timeGetTime() - start_time;
	#endif
#endif
}

#endif


	void update_best_dominating_set()
	{
		for(int v = 1; v <= v_num; v++)
		{
			best_in_dominating_set[v] = in_dominating_set[v];
		}
	}

	void add(int v)
	{
#ifdef debug_mode
cout << "==================add " << v << " gain: " << score[v] << endl;
#endif
		in_dominating_set[v] = 1;
		if(!connect_dominating_set_degree[v])//may add a covered vertex
			ptr_to_uncov_vertices->delete_element(v);
		dominating_set_size++;
		ptr_to_partitioned_dominating_set->placeInVertexToCover(v, score);
		int* nbs = vertices[v].get_neighbors();
		int* adj_edgs = vertices[v].get_adj_edges();
		int dgr = vertices[v].get_degree();
		for(int i = 0; i < dgr; i++)
		{
			int n = nbs[i];
			connect_dominating_set_degree[n]++;
			int e = adj_edgs[i];
			if(in_dominating_set[n])
			{
				if(connect_dominating_set_degree[v] == 1)// v becomes strongly covered
					ptr_to_partitioned_dominating_set->lossMinusMinus(n, score);
				if(connect_dominating_set_degree[n] == 1)// n becomes strongly covered
					ptr_to_partitioned_dominating_set->lossMinusMinus(n, score);// one vertex may cover two vertices
			}
			else
			{
				if(connect_dominating_set_degree[v] == 0)// v becomes covered
				{
					ptr_to_partitioned_dominating_set->gainMinusMinus(n, score);
#ifdef score_changed_mode
					confChange[n] = 1;
#endif
				}
				if(connect_dominating_set_degree[n] == 2)// n becomes strongly covered
				{
					int* nbs_ = vertices[n].get_neighbors();
					int* adj_edgs_ = vertices[n].get_adj_edges();
					int dgr_ = vertices[n].get_degree();
					for(int j = 0; j < dgr_; j++)
					{
						int n_ = nbs_[j];
						int e_ = adj_edgs_[j];
						if(n_ != v && in_dominating_set[n_])
						{
							ptr_to_partitioned_dominating_set->lossMinusMinus(n_, score);// no longer causes loss of n
							break;
						}
					}							
				}
				else if(connect_dominating_set_degree[n] == 1)// n becomes critically covered
				{
					ptr_to_partitioned_dominating_set->gainMinusMinus(n, score);//n itself no longer causes gain of n
#ifdef score_changed_mode
					confChange[n] = 1;
#endif
					ptr_to_uncov_vertices->delete_element(n);
					int* nbs_ = vertices[n].get_neighbors();
					int* adj_edgs_ = vertices[n].get_adj_edges();
					int dgr_ = vertices[n].get_degree();
					for(int j = 0; j < dgr_; j++)
					{
						int n_ = nbs_[j];
						int e_ = adj_edgs_[j];
						if(!in_dominating_set[n_])//n_ != v
						{
							ptr_to_partitioned_dominating_set->gainMinusMinus(n_, score);//n_ no longer causes gain of n
#ifdef score_changed_mode
							confChange[n_] = 1;
#endif
						}
					}
				}
#ifdef nvcc_mode
				confChange[n] = 1;
#endif
			}
		}
		time_stamp[v] = step;
#ifdef domination_hash_mode
		ptr_to_hashed_domination->update_hash_wrt_add(v);
#endif
#ifdef ans_update_op_mode
		ptr_to_moved_v->ans_update(v);
#endif
//show_state();
	}

	void remove(int v)
	{
#ifdef debug_mode
cout << "==================remove " << v << " loss: " << -score[v] << endl;
#endif
#ifdef debug_mode
//cout << "1. in remove func, in localSearch.h, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
cout << 1 << endl;
if(!ptr_to_partitioned_dominating_set->checkMap(v_num, vertices)) exit(1);
if(!ptr_to_partitioned_dominating_set->checkPartition(v_num, score, vertices)) exit(1);
/*
if(ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] > v_num)
{
	cout << "in localSearch.h, containing abnormal element, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
	exit(1);
}
*/
#endif
		confChange[v] = 0;
		in_dominating_set[v] = 0;
		if(!connect_dominating_set_degree[v])
			ptr_to_uncov_vertices->insert_element(v);
		dominating_set_size--;
		ptr_to_partitioned_dominating_set->placeOutVertexFromCover(v, score);
#ifdef debug_mode
//cout << "1.5 in remove func, in localSearch.h, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
cout << 1.5 << endl;
if(!ptr_to_partitioned_dominating_set->checkMap(v_num, vertices)) exit(1);
if(!ptr_to_partitioned_dominating_set->checkPartition(v_num, score, vertices)) exit(1);
/*
if(ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] > v_num)
{
	cout << "in localSearch.h, containing abnormal element, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
	exit(1);
}
*/
#endif
		int* nbs = vertices[v].get_neighbors();
		int* adj_edgs = vertices[v].get_adj_edges();
		int dgr = vertices[v].get_degree();
		for(int i = 0; i < dgr; i++)
		{
			int n = nbs[i];
#ifdef debug_mode
cout << "***considering " << n << endl;
if(n == 413) {cout << "found 413" << endl; getchar();}
#endif
			connect_dominating_set_degree[n]--;
			int e = adj_edgs[i];
			if(in_dominating_set[n])
			{
#ifdef debug_mode
//cout << "1.5 in remove func, in localSearch.h, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
cout << 1.51 << endl;
cout << "n = " << n << endl;
if(!ptr_to_partitioned_dominating_set->checkMap(v_num, vertices)) exit(1);
if(!ptr_to_partitioned_dominating_set->checkPartition(v_num, score, vertices)) exit(1);
/*
if(ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] > v_num)
{
	cout << "in localSearch.h, containing abnormal element, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
	exit(1);
}
*/
#endif
				if(connect_dominating_set_degree[v] == 1)// v becomes critically covered
					ptr_to_partitioned_dominating_set->lossPlusPlus(n, score);
#ifdef debug_mode
//cout << "1.5 in remove func, in localSearch.h, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
cout << 1.515 << endl;
if(!ptr_to_partitioned_dominating_set->checkMap(v_num, vertices)) exit(1);
if(!ptr_to_partitioned_dominating_set->checkPartition(v_num, score, vertices)) exit(1);
/*
if(ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] > v_num)
{
	cout << "in localSearch.h, containing abnormal element, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
	exit(1);
}
*/
#endif
				if(connect_dominating_set_degree[n] == 0)// n becomes critically covered
					ptr_to_partitioned_dominating_set->lossPlusPlus(n, score);
#ifdef debug_mode
//cout << "1.5 in remove func, in localSearch.h, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
cout << 1.52 << endl;
if(!ptr_to_partitioned_dominating_set->checkMap(v_num, vertices)) exit(1);
if(!ptr_to_partitioned_dominating_set->checkPartition(v_num, score, vertices)) exit(1);
/*
if(ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] > v_num)
{
	cout << "in localSearch.h, containing abnormal element, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
	exit(1);
}
*/
#endif
			}
			else
			{
#ifdef debug_mode
//cout << "1.5 in remove func, in localSearch.h, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
cout << 1.53 << endl;
if(!ptr_to_partitioned_dominating_set->checkMap(v_num, vertices)) exit(1);
if(!ptr_to_partitioned_dominating_set->checkPartition(v_num, score, vertices)) exit(1);
/*
if(ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] > v_num)
{
	cout << "in localSearch.h, containing abnormal element, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
	exit(1);
}
*/
#endif
				if(connect_dominating_set_degree[v] == 0)// v becomes uncovered
				{
					ptr_to_partitioned_dominating_set->gainPlusPlus(n, score);
#ifdef score_changed_mode
					confChange[n] = 1;
#endif
				}
				if(connect_dominating_set_degree[n] == 0)// n becomes uncovered
				{
					ptr_to_partitioned_dominating_set->gainPlusPlus(n, score);
#ifdef score_changed_mode
					confChange[n] = 1;
#endif
					ptr_to_uncov_vertices->insert_element(n);
					int* nbs_ = vertices[n].get_neighbors();
					int* adj_edgs_ = vertices[n].get_adj_edges();
					int dgr_ = vertices[n].get_degree();
					for(int j = 0; j < dgr_; j++)
					{
						int n_ = nbs_[j];
#ifdef debug_mode
cout << "\tconsidering " << n_ << endl;
if(n_ == 413) {cout << "found 413" << endl; getchar();}
#endif
						int e_ = adj_edgs_[j];
						if(!in_dominating_set[n_] && n_ != v)
						{
							ptr_to_partitioned_dominating_set->gainPlusPlus(n_, score);// undirect effects
#ifdef score_changed_mode
							confChange[n_] = 1;
#endif
						}
					}
				}
				else if(connect_dominating_set_degree[n] == 1)
				{
					int* nbs_ = vertices[n].get_neighbors();
					int* adj_edgs_ = vertices[n].get_adj_edges();
					int dgr_ = vertices[n].get_degree();
					for(int j = 0; j < dgr_; j++)
					{
						int n_ = nbs_[j];
#ifdef debug_mode
cout << "\tconsidering " << n_ << endl;
if(n_ == 413) {cout << "found 413" << endl; getchar();}
#endif
						int e_ = adj_edgs_[j];
						if(in_dominating_set[n_])
						{
							ptr_to_partitioned_dominating_set->lossPlusPlus(n_, score);
							break;
						}
					}					
				}
#ifdef nvcc_mode
				confChange[n] = 1;
#endif
#ifdef debug_mode
//cout << "1.5 in remove func, in localSearch.h, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
cout << 1.54 << endl;
if(!ptr_to_partitioned_dominating_set->checkMap(v_num, vertices)) exit(1);
if(!ptr_to_partitioned_dominating_set->checkPartition(v_num, score, vertices)) exit(1);
/*
if(ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] > v_num)
{
	cout << "in localSearch.h, containing abnormal element, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
	exit(1);
}
*/
#endif
			}
		}
		time_stamp[v] = step;
#ifdef domination_hash_mode
		ptr_to_hashed_domination->update_hash_wrt_remove(v);
#endif
#ifdef ans_update_op_mode
		ptr_to_moved_v->ans_update(v);
#endif
//show_state();
#ifdef debug_mode
//cout << "2. in remove func, in localSearch.h, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
cout << "2." << endl;
if(!ptr_to_partitioned_dominating_set->checkMap(v_num, vertices)) exit(1);
if(!ptr_to_partitioned_dominating_set->checkPartition(v_num, score, vertices)) exit(1);
/*
if(ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] > v_num)
{
	cout << "in localSearch.h, containing abnormal element, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
	exit(1);
}
*/
#endif
	}

#ifdef edge_time_stamp_mode
	HugeInt aged_gain(int v)
	{
		HugeInt a_gain = 0;
		int* nbs = vertices[v].get_neighbors();
		int* adj_edgs = vertices[v].get_adj_edges();
		int dgr = vertices[v].get_degree();
		for(int i = 0; i < dgr; i++)
		{
			int n = nbs[i];
			if(in_dominating_set[n] == 1) continue;
			int e = adj_edgs[i];
			a_gain += (step - edge_time_stamp[e]);		
		}
		return a_gain;
	}
#endif

	int rand_vertex_in_dominating_set()
	{
		int v;
		do{
			v = rand() % v_num + 1;
		}while(!in_dominating_set[v]);
		return v;
	}

	int rand_vertex_outside_dominating_set()
	{
		int v;
		do{
			v = rand() % v_num + 1;
		}while(in_dominating_set[v]);
		return v;
	}

	int rand_cov_vertex()
	{
		int v;
		do{
			v = rand() % v_num + 1;
		}while(ptr_to_uncov_vertices->element_in(v));
		return v;
	}

	int rand_dominator(int v)// v must be covered
	{
		int *nbs = vertices[v].get_neighbors();
		int dgr = vertices[v].get_degree();
		int *cand_vertices = new int[dgr + 1];
		int cand_num = 0;

		for(int i = 0; i < dgr; i++)
		{
			int n = nbs[i];
			if(in_dominating_set[n]) cand_vertices[cand_num++] = n;
		}
		if(in_dominating_set[v]) cand_vertices[cand_num++] = v;
		
		int u = cand_vertices[rand() % cand_num];
		delete[] cand_vertices;
		return u;
	}

#ifdef probablistic_remove_mode
	int probabilistic_vertex_in_dominating_set_wrt_to_degree()
	{
		int v;
		do{
			v = multi_vertex_set[rand() % (e_num << 1)];
		}
		while(!in_diminating_set[v])
		return v;
	}
#endif

	void local_move()
	{
		int i;
		int u, v;
		int* nbs;
		int dgr;
		if(!ptr_to_uncov_vertices->size())
		{
#ifdef ans_update_op_mode
				ptr_to_moved_v->efficiently_update_best_dominating_set_vertices();
				best_dominating_set_size = dominating_set_size;
#else
				update_best_dominating_set();
				best_dominating_set_size = dominating_set_size;
#endif
#ifdef linux_mode
				times(&finish);
				best_cmp_time = double((finish.tms_utime + finish.tms_stime - start_time)) / sysconf(_SC_CLK_TCK);
#endif
#ifdef windows_mode
				best_cmp_time = timeGetTime() - start_time;
#endif
				best_solve_step = step;
				u = ptr_to_partitioned_dominating_set->randMinLossVertexInDominatingSet();
				remove(u);
//show_state();
				return;
		}

#ifdef domination_hash_mode
		ptr_to_hashed_domination->mark_hash_entry();
//cout << "curr_hash_entry: " << ptr_to_hashed_domination->get_hash_entry() << endl;
		if(ptr_to_hashed_domination->curr_hash_entry_marked_too_frequently())
		{
//cout << "remarked" << endl;
				int x, y, z;
#ifdef probabilistic_remove_mode
				y = probabilistic_vertex_in_dominating_set_wrt_to_degree();
#endif
#ifdef double_random_remove_mode
				x = rand_cov_vertex();
				y = rand_dominator(x);
#endif
				remove(y);
				z = rand_vertex_outside_dominating_set();
				add(z);
//show_state();
				return;
		}
#endif
#ifdef debug_mode
cout << "step: " << step << endl;
cout << "in localSearch.h, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
if(!ptr_to_partitioned_dominating_set->checkMap(v_num, vertices)) exit(1);
if(!ptr_to_partitioned_dominating_set->checkPartition(v_num, score, vertices)) exit(1);
if(ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] > v_num)
{
	cout << "in localSearch.h, containing abnormal element, at location 0, the element is " << ptr_to_partitioned_dominating_set->get_v_ascending_degree_array()[0] << endl;
	exit(1);
}
#endif
		u = ptr_to_partitioned_dominating_set->randMinLossVertexInDominatingSet();
#ifdef debug_mode
//cout << "to remove " << u << endl;
cout << 2.001 << endl;
#endif
		remove(u);	
#ifdef debug_mode
cout << 2.002 << endl;
#endif
		int best_add_v;

		best_add_v = ptr_to_uncov_vertices->rand_element();
#ifdef debug_mode
cout << "potential best_add_v: " << best_add_v << endl;
#endif
		nbs = vertices[best_add_v].get_neighbors();
		dgr = vertices[best_add_v].get_degree();
#ifdef debug_mode
cout << 2.01 << endl;
#endif
		if(confChange[best_add_v])
		{
#ifdef debug_mode
cout << 2.1 << endl;
#endif
			for(i = 0; i < dgr; i++)
			{
				int n = nbs[i];
				if(!confChange[n]) continue;
	//cout << "considering " << n << endl;
				if(score[n] > score[best_add_v] || (score[n] == score[best_add_v] && time_stamp[n] < score[best_add_v]))
					best_add_v = n;
			}
		}
		else
		{
#ifdef debug_mode
cout << 2.2 << endl;
#endif
			for(i = 0; i < dgr; i++)
			{
				if(confChange[nbs[i]]) break;
			}
//cout << 2.21 << endl;
//cout << "dgr: " << dgr << endl;
//cout << "i: " << i << endl;
//cout << "i should be smaller than dgr" << endl;
			best_add_v = nbs[i];
			i++;
//cout << 2.22 << endl;
			for(; i < dgr; i++)
			{
//cout << 2.23 << endl;
				int n = nbs[i];
				if(!confChange[n]) continue;
//cout << 2.24 << endl;
	//cout << "considering " << n << endl;
				if(score[n] > score[best_add_v] || (score[n] == score[best_add_v] && time_stamp[n] < score[best_add_v]))
					best_add_v = n;
			}
		}
		
		add(best_add_v);
#ifdef debug_mode
cout << 3 << endl;
#endif
		step++;
#ifdef debug_mode
show_state();
#endif
	}



public:

	StochasticLocalSearch(char *filename, int cut_off) : Graph(filename)
	{
//cout << "constructing stochastic local search" << endl;
#ifdef debug_mode
cout << "begin to construct" << endl;
#endif
		step = 0;
		time_limit = cut_off;
#ifdef linux_mode
		times(&start);
		start_time = start.tms_utime + start.tms_stime;
#endif
#ifdef windows_mode
		start_time = timeGetTime();
#endif
#ifdef individual_analysis_on_init_sls_mode
		init_time = 0;
#endif

//cout << 1 << endl;
		in_dominating_set = new int[v_num + 2];
		memset(in_dominating_set, 0, sizeof(int) * (v_num + 2));
		dominating_set_size = 0;
		connect_dominating_set_degree = new int[v_num + 2];
		memset(connect_dominating_set_degree, 0, sizeof(int) * (v_num + 2));
//cout << 1.1 << endl;
//cout << "get_max_degree(): " << get_max_degree() << endl;
		//ptr_to_partitioned_dominating_set = new ScoreBasedPartition(v_num, get_max_degree(), vertices);
		ptr_to_partitioned_dominating_set = new DS_Partition(v_num, get_max_degree(), vertices);
//cout << 1.2 << endl;
		//ptr_to_uncov_edges = new Bijection(e_num);
		ptr_to_uncov_vertices = new ConfchangedScoreAgePickSet(v_num);
//cout << 2 << endl;
		score = new int[v_num + 2];
		for(int v = 1; v <= v_num; v++)
		{
				score[v] = vertices[v].get_degree() + 1;//
		}

		time_stamp = new HugeInt[v_num + 2];
#ifdef edge_time_stamp_mode
		edge_time_stamp = new HugeInt[e_num + 2];
		for(int e = 0; e < e_num; e++)
			edge_time_stamp[e] = 0;
#endif
		confChange = new int[v_num + 2];
		for(int v = 1; v <= v_num; v++)
		{
#ifdef age_init_i_mode
			time_stamp[v] = -v_num - 1 + v;
#endif
			confChange[v] = 1;
		}

#ifdef domination_hash_mode
		ptr_to_hashed_domination = new DominationHash(v_num);
		consecutive_collision_num = 0;
#endif

#ifdef probablistic_remove_mode
		multi_vertex_set = new int[e_num * 2];
		{
			int i = 0;// cursor in the multi_vertex_set
			for(int v = 1; v <= v_num; v++)
			{
				int j = 0;// for repeating vertices
				while(j < vertices[v].get_degree())
				{
					multi_vertex_set[i++] = v;
					j++;
				}
			}
		}
#endif

#ifdef ans_update_op_mode
		ptr_to_moved_v = new AnsChangeSet(v_num);
#endif
//cout << 3 << endl;
		best_in_dominating_set = new int[v_num + 2];
		memset(best_in_dominating_set, 0, sizeof(int) * (v_num + 2));
		best_dominating_set_size = 0;
#ifdef debug_mode
//show_state();
cout << "still in contruction function" << endl;
#endif
		init_solution();
#ifdef debug_mode
cout << "construction completed" << endl;
#endif
	}


	~StochasticLocalSearch()
	{
		delete[] in_dominating_set;
		delete[] connect_dominating_set_degree;
		delete ptr_to_partitioned_dominating_set;
		delete ptr_to_uncov_vertices;
		delete[] time_stamp;
#ifdef edge_time_stamp_mode
		delete[] edge_time_stamp;
#endif
		delete[] score;
		delete[] confChange;
#ifdef domination_hash_mode
		delete ptr_to_hashed_domination;
#endif
#ifdef probablistic_remove_mode
		delete[] multi_vertex_set;
#endif
		delete[] best_in_dominating_set;
#ifdef ans_update_op_mode
		delete ptr_to_moved_v;
#endif
	}

	void cover_sls()
	{
		while(1)
		{
			if(step % TRY_STEP == 0)
			{
#ifdef linux_mode
					times(&finish);
					double elap_time = (finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
#endif
#ifdef windows_mode
					unsigned long elap_time = timeGetTime() - start_time;
#endif
					if(elap_time >= time_limit)
					{
#ifdef ans_update_op_mode
						ptr_to_moved_v->compute_final_answer(v_num, in_dominating_set, best_in_dominating_set);
#endif
						return;
					}
			}
			local_move();
		}
	}

	void show_results()
	{
#ifdef linux_mode
		times(&finish);
		double elap_time = double(finish.tms_utime + finish.tms_stime - start_time) / sysconf(_SC_CLK_TCK);
#endif
#ifdef windows_mode
		double elap_time = timeGetTime() - start_time;
#endif
		sls_time = elap_time - init_time;
		if(check_solution())
		{
			cout << "o " << best_dominating_set_size << endl;
			cout << "c solveTime " << best_cmp_time << endl;
			cout << "c searchSteps " << best_solve_step << endl;
			cout << "c stepSpeed(/ms) " << step / 1000.0 / sls_time << endl;
			cout << "c dominating set: " << endl;
			for(int v = 1; v <= v_num; v++) if(best_in_dominating_set[v]) cout << v << '\t'; 
			cout << endl;
#ifdef strategy_analysis_mode
cout << "the total number of steps: " << step << endl;
cout << "answer update times: " << ans_update_times << endl;
#endif
		}
		else
		{
			cout << "sorry, something is wrong" << endl;
		}
	}	

	void show_state()
	{
		int v, e;
/*
		for(e = 0; e < e_num; e++)
		{
			cout << edges[e].get_v1() << " " << edges[e].get_v2() << endl;
		}

		for(v = 1; v <= v_num; v++)
		{
			cout << "vertex " << v << ": " << endl;
			vertices[v].show_neighbors();
		}
*/
		int i, j;
		//cout << "step: " << step << endl;
#if 0
		cout << "the dominating set: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			if(in_dominating_set[v])
				cout << v << '\t';
		}
		cout << endl;
#endif
#if 0
		cout << "loss: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			if(in_dominating_set[v])
				cout << -score[v] << '\t';
		}
		cout << endl;
		cout << "connect dominating set degrees: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << v << '\t';
		}
		cout << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << connect_dominating_set_degree[v] << '\t';
		}
		cout << endl;
		cout << "vertices: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << v << '\t';
		}
		cout << endl;
		cout << "score: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << score[v] << '\t';
		}
		cout << endl;
		cout << "confChange: " << endl;
		for(v = 1; v <= v_num; v++)
		{
			cout << confChange[v] << '\t';
		}
		cout << endl;
		cout << "uncov vertices: " << endl;
		for(i = ptr_to_uncov_vertices->begin(); i < ptr_to_uncov_vertices->end(); i++)
		{
			int v = ptr_to_uncov_vertices->at(i);
			cout << v << '\t';
		}
		cout << endl;
		cout << "max_gain: " << ptr_to_partitioned_dominating_set->get_max_gain() << endl;
/*
		cout << "confChange: " << endl;
		for(i = ptr_to_uncov_vertices->begin(); i < ptr_to_uncov_vertices->end(); i++)
		{
			int v = ptr_to_uncov_vertices->at(i);
			cout << confChange[v] << '\t';
		}
		cout << endl;
*/
/*
		cout << "uncov edges: " << endl;
		for(i = ptr_to_uncov_edges->begin(); i < ptr_to_uncov_edges->end(); i++)
		{
			int e = ptr_to_uncov_edges->at(i);
			cout << e;
			int v1, v2;
			edges[e].get_vertices(v1, v2);
			cout << ", " << "endpoints: " << v1 << " and " << v2 << endl;
		}
		cout << endl;
*/
	
cout << "**************************************" << endl;
#endif
//getchar();
	}
};
#endif
