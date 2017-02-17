#ifndef DS_BUCKETS_H
#define DS_BUCKETS_H

#include "myBuckets.h"

class DS_Partition : public ScoreBasedPartition
{
private:
	int remaining_vertex_num;	
public:
	DS_Partition(int v_num, int max_degree, Vertex* vertices):ScoreBasedPartition(v_num, max_degree, vertices)
	{
		remaining_vertex_num = v_num;
	}

	~DS_Partition()
	{
	}

	void filter_out_sole_vertices(Vertex* vertices, int v_num)
	{
		int i, j, l;
		l = 1;
		i = j = ptr_to_v_loss[1];
		while(i < v_num)
		{
			while(i < ptr_to_v_loss[l+1])
			{
				int v = v_ascending_degree_array[i++];
				if(l == 1 && vertices[v].get_degree() == 0)
				{
#ifdef debug_mode
//cout << "filter out " << v << endl;
#endif
					continue;
				}
				v_ascending_degree_array[j] = v;
				index_in_v_ascending_degree_array[v] = j;
				j++;
			}
			ptr_to_v_loss[l+1] = j;
			ptr_to_v_gain[l] = j;
			l++;
		}
		while(l < max_degree + 1 + 1)
		{
			ptr_to_v_loss[l+1] = j;
			ptr_to_v_gain[l] = j;
			l++;
		}
		remaining_vertex_num = j - 1;
	}

	int get_remaining_vertex_num()
	{
		return remaining_vertex_num;
	}

int checkMap(int v_num, Vertex* vertices)
//very time-consuming if it exists in frequently called functions or loops
{
#ifdef debug_mode
cout << "checking map in dsBuckets.h" << endl;
#endif
	int v;
	for(v = 1; v <= v_num; ++v)
	{
		if(vertices[v].get_degree() == 0) continue;
		if(v_ascending_degree_array[index_in_v_ascending_degree_array[v]] != v)
		{
			cout << "vertex " << v << " error!" << endl;
			cout << "location " << index_in_v_ascending_degree_array[v] << " error!" << endl;
			return 0;
		}
	}
#ifdef debug_mode
cout << "having checked map in dsBuckets.h" << endl;
#endif
	return 1;
}

int checkPartition(int v_num, int* score, Vertex* vertices)
//very time-consuming if it exists in frequently called functions or loops
{
	int v;

	int i;
#ifdef debug_mode
cout << "checking partition in dsBuckets.h" << endl;
#endif
/*
	cout << "vertices: " << endl;
	for(i = 0; i < remaining_vertex_num; i++)
	{
		cout << v_ascending_degree_array[i] << '\t';
	}
	cout << endl;
	cout << "score: " << endl;
	for(i = 0; i < v_num; i++)
	{
		cout << score[v_ascending_degree_array[i]] << '\t';
	}
	cout << endl;
*/
	for(v = 1; v <= v_num; ++v)
	{
		if(vertices[v].get_degree() == 0) continue;
		if(score[v] > 0)
		{
			if(index_in_v_ascending_degree_array[v] < ptr_to_v_gain[score[v]] || index_in_v_ascending_degree_array[v] >= ptr_to_v_loss[score[v]+1])
			{
				cout << "1. vertex  " << v << " has errors" << endl;
				cout << "index_in_v_ascending_degree_array of " << v << " is " << index_in_v_ascending_degree_array[v] << endl;
				cout << "ptr_to_v_gain: " << ptr_to_v_gain[score[v]] << endl;
				cout << "ptr_to_v_loss: " << ptr_to_v_loss[score[v]+1] << endl;
				return 0;
			}
		}
		else if(score[v] < 0)
		{
			if(index_in_v_ascending_degree_array[v] < ptr_to_v_loss[-score[v]] || index_in_v_ascending_degree_array[v] >= ptr_to_v_gain[-score[v]])
			{
				cout << "2. vertex  " << v << " has errors" << endl;
				cout << "index_in_v_ascending_degree_array of " << v << " is " << index_in_v_ascending_degree_array[v] << endl;
				cout << "the score of " << v << " is " << score[v] << endl;
				cout << "ptr_to_v_loss: " << ptr_to_v_loss[-score[v]] << endl;
				cout << "ptr_to_v_gain: " << ptr_to_v_gain[-score[v]] << endl;
				return 0;
			}
		}
		else
		{
			if(index_in_v_ascending_degree_array[v] >= ptr_to_v_loss[1])
			{
				cout << "3. vertex  " << v << " has errors" << endl;
				return 0;
			}
		}
	}
#ifdef debug_mode
cout << "having checked partition in dsBuckets.h" << endl;
#endif
	return 1;
}
};
#endif
