#ifndef _MY_BUCKETS_H
#define _MY_BUCKETS_H

#include <iostream>

#include "graph.h"
#include "hugeInt.h"

//#define debug_mode

using namespace std;

class ScoreBasedPartition{
protected:
	int		*v_ascending_degree_array;//[MAXV]
	int		*index_in_v_ascending_degree_array;//[MAXV]
	int		*ptr_to_v_gain;//[MAXV]
	int		*ptr_to_v_loss;//[MAXV]
	int		max_degree;
	int		max_gain;
public:
	ScoreBasedPartition(int v_num, int max_d, Vertex* vertices)
	{
//cout << 1.11 << endl;
		max_degree = max_d;
		v_ascending_degree_array = new int[v_num + 2];
		index_in_v_ascending_degree_array = new int[v_num + 2];
		ptr_to_v_gain = new int[max_degree + 1 + 2];
		ptr_to_v_loss = new int[max_degree + 1 + 2];

		int* v_gain_num = new int[max_degree + 1 + 2];
		int *qtr_to_v_gain = new int[max_degree + 1 + 2];
//cout << 1.12 << endl;
		int i;
		int v;
		int s;

		//statistics
		memset(v_gain_num, 0, sizeof(int) * (max_degree + 1 + 2));
		for(v = 1; v <= v_num; ++v)
		{
			s = vertices[v].get_degree() + 1;
			v_gain_num[s]++;
		}
	// partition
		ptr_to_v_gain[0] = 0;
		ptr_to_v_loss[0] = 0;
		qtr_to_v_gain[0] = 0;
//cout << 1.13 << endl;
		for(s = 1; s < max_degree + 1 + 2; s++)//
		{
			ptr_to_v_gain[s] = ptr_to_v_gain[s-1] + v_gain_num[s-1];
			ptr_to_v_loss[s] = ptr_to_v_gain[s];
			qtr_to_v_gain[s] = ptr_to_v_gain[s];
		}
//cout << 1.14 << endl;
	// place in
		for(v = 1; v <= v_num; ++v)
		{
			s = vertices[v].get_degree() + 1;
			// place v in the partition with gain g
			v_ascending_degree_array[qtr_to_v_gain[s]] = v;
			index_in_v_ascending_degree_array[v] = qtr_to_v_gain[s];
			qtr_to_v_gain[s]++;		
		}
//cout << 1.15 << endl;	
		delete[] qtr_to_v_gain;
		delete[] v_gain_num;
//cout << 1.16 << endl;
		max_gain = max_degree + 1;
	}
	~ScoreBasedPartition()
	{
		delete[] v_ascending_degree_array;
		delete[] index_in_v_ascending_degree_array;
		delete[] ptr_to_v_gain;
		delete[] ptr_to_v_loss;		
	}

int* get_v_ascending_degree_array()
{
	return v_ascending_degree_array;
}

int* get_index_in_v_ascending_degree_array()
{
	return index_in_v_ascending_degree_array;
}

int* get_ptr_to_v_gain()
{
	return ptr_to_v_gain;
}

int* get_ptr_to_v_loss()
{
	return ptr_to_v_loss;
}

int get_max_gain()
{
	return max_gain;
}

void placeOutVertexFromCover(int v, int* score)
{
#ifdef debug_mode
//cout << "2. in myBuckets.h, at location 0, the element is " << v_ascending_degree_array[0] << endl;
if(v_ascending_degree_array[0] > 6100)
{
	cout << "0.01 in myBuckets.h, containing abnormal element, at location 0, the element is " << v_ascending_degree_array[0] << endl;
	//exit(1);
}
#endif
	int vPtr = index_in_v_ascending_degree_array[v];
	int boundary_loc = ptr_to_v_gain[-score[v]] - 1;
	int boundary_v = v_ascending_degree_array[boundary_loc];
#ifdef debug_mode
if(boundary_loc < 0)
{
	cout << "1. boundary_loc < 0" << endl;
	exit(1);
}
#endif
	//
	v_ascending_degree_array[vPtr] = boundary_v;// already had a copy, don't need a tmp vertex
	index_in_v_ascending_degree_array[boundary_v] = vPtr;
	//
	v_ascending_degree_array[boundary_loc] = v;
	index_in_v_ascending_degree_array[v] = boundary_loc;
#ifdef debug_mode
if(v_ascending_degree_array[0] > 6100)
{
	cout << "0.015 in myBuckets.h, containing abnormal element, at location 0, the element is " << v_ascending_degree_array[0] << endl;
	cout << "vPtr: " << vPtr << endl;
	cout << "boundary_loc: " << boundary_loc << endl;
	//exit(1);
}
#endif
	//
	ptr_to_v_gain[-score[v]]--;
	//
	score[v] = -score[v];

	if(score[v] > max_gain)
		max_gain = score[v];
#ifdef debug_mode
//cout << "2. in myBuckets.h, at location 0, the element is " << v_ascending_degree_array[0] << endl;
if(v_ascending_degree_array[0] > 6100)
{
	cout << "0.02 in myBuckets.h, containing abnormal element, at location 0, the element is " << v_ascending_degree_array[0] << endl;
	//exit(1);
}
#endif
}

void placeInVertexToCover(int v, int* score)
{
	int vPtr = index_in_v_ascending_degree_array[v];
	int boundary_loc = ptr_to_v_gain[score[v]];
	int boundary_v = v_ascending_degree_array[boundary_loc];
//cout << "1. index of " << v << " in ascending degree array is " << index_in_v_ascending_degree_array[v] << endl;
#ifdef debug_mode
if(boundary_loc < 0)
{
	cout << "2. boundary_loc < 0" << endl;
	exit(1);
}
#endif
	//
	v_ascending_degree_array[vPtr] = boundary_v;// already had a copy, don't need a tmp vertex
	index_in_v_ascending_degree_array[boundary_v] = vPtr;
	//
//cout << "2. index of " << v << " in ascending degree array is " << index_in_v_ascending_degree_array[v] << endl;
	v_ascending_degree_array[boundary_loc] = v;
	index_in_v_ascending_degree_array[v] = boundary_loc;
	//
//cout << "3. index of " << v << " in ascending degree array is " << index_in_v_ascending_degree_array[v] << endl;
	ptr_to_v_gain[score[v]]++;
	//
	score[v] = -score[v];

	while(ptr_to_v_gain[max_gain] == ptr_to_v_loss[max_gain + 1])
		max_gain--;
#ifdef debug_mode
//cout << "1. in myBuckets.h, at location 0, the element is " << v_ascending_degree_array[0] << endl;
if(v_ascending_degree_array[0] > 6100)
{
	cout << "1. in myBuckets.h, containing abnormal element, at location 0, the element is " << v_ascending_degree_array[0] << endl;
	//exit(1);
}
#endif
}

void gainPlusPlus(int v, int* score)
{
#ifdef debug_mode
//cout << "2. in myBuckets.h, at location 0, the element is " << v_ascending_degree_array[0] << endl;
if(v_ascending_degree_array[0] > 6100)
{
	cout << "2.01 in myBuckets.h, containing abnormal element, at location 0, the element is " << v_ascending_degree_array[0] << endl;
	//exit(1);
}
#endif
	score[v]++;
	int vPtr = index_in_v_ascending_degree_array[v];
	int boundary_loc_1 = ptr_to_v_loss[score[v]] - 1;// loc of v41
	int boundary_v_1 = v_ascending_degree_array[boundary_loc_1];
	int boundary_loc_2 = ptr_to_v_gain[score[v]] - 1;// loc of v68
	int boundary_v_2 = v_ascending_degree_array[boundary_loc_2];
#ifdef debug_mode
if(boundary_loc_1 < 0)
{
	cout << "1. boundary_loc_1 < 0" << endl;
	exit(1);
}
if(boundary_loc_2 < 0)
{
	cout << "1. boundary_loc_2 < 0" << endl;
	exit(1);
}
#endif
	//
	v_ascending_degree_array[vPtr] = boundary_v_1;
	v_ascending_degree_array[boundary_loc_1] = boundary_v_2;
	v_ascending_degree_array[boundary_loc_2] = v;
	//
	index_in_v_ascending_degree_array[boundary_v_2] = boundary_loc_1;
	index_in_v_ascending_degree_array[boundary_v_1] = vPtr;
	index_in_v_ascending_degree_array[v] = boundary_loc_2;
	//
	ptr_to_v_loss[score[v]]--;
	ptr_to_v_gain[score[v]]--;
	
	if(score[v] > max_gain)
		max_gain = score[v];
#ifdef debug_mode
//cout << "2. in myBuckets.h, at location 0, the element is " << v_ascending_degree_array[0] << endl;
if(v_ascending_degree_array[0] > 6100)
{
	cout << "2.02 in myBuckets.h, containing abnormal element, at location 0, the element is " << v_ascending_degree_array[0] << endl;
	//exit(1);
}
#endif
}

void lossMinusMinus(int v, int* score)
{
	int vPtr = index_in_v_ascending_degree_array[v];
	int boundary_loc_1 = ptr_to_v_loss[-score[v]];// loc of v68
	int boundary_v_1 = v_ascending_degree_array[boundary_loc_1];
	int boundary_loc_2 = ptr_to_v_gain[-score[v]-1];// loc of v99
	int boundary_v_2 = v_ascending_degree_array[boundary_loc_2];
#ifdef debug_mode
if(boundary_loc_1 < 0)
{
	cout << "2. boundary_loc_1 < 0" << endl;
	exit(1);
}
if(boundary_loc_2 < 0)
{
	cout << "2. boundary_loc_2 < 0" << endl;
	exit(1);
}
#endif
	//
	v_ascending_degree_array[vPtr] = boundary_v_1;
	v_ascending_degree_array[boundary_loc_1] = boundary_v_2;
	v_ascending_degree_array[boundary_loc_2] = v;
	//
	index_in_v_ascending_degree_array[boundary_v_2] = boundary_loc_1;
	index_in_v_ascending_degree_array[boundary_v_1] = vPtr;
	index_in_v_ascending_degree_array[v] = boundary_loc_2;
	//
	ptr_to_v_loss[-score[v]]++;
	ptr_to_v_gain[-score[v]-1]++;
	//
	score[v]++;
#ifdef debug_mode
//cout << "3. in myBuckets.h, at location 0, the element is " << v_ascending_degree_array[0] << endl;
if(v_ascending_degree_array[0] > 6100)
{
	cout << "3. in myBuckets.h, containing abnormal element, at location 0, the element is " << v_ascending_degree_array[0] << endl;
	//exit(1);
}
#endif	
}

void gainMinusMinus(int v, int* score)
{
	int vPtr = index_in_v_ascending_degree_array[v];
	int boundary_loc_1 = ptr_to_v_gain[score[v]];
	int boundary_v_1 = v_ascending_degree_array[boundary_loc_1];
	int boundary_loc_2 = ptr_to_v_loss[score[v]];
	int boundary_v_2 = v_ascending_degree_array[boundary_loc_2];
#ifdef debug_mode
if(boundary_loc_1 < 0)
{
	cout << "3. boundary_loc_1 < 0" << endl;
	exit(1);
}
if(boundary_loc_2 < 0)
{
	cout << "3. boundary_loc_2 < 0" << endl;
	exit(1);
}
#endif
	//
	v_ascending_degree_array[vPtr] = boundary_v_1;
	v_ascending_degree_array[boundary_loc_1] = boundary_v_2;
	v_ascending_degree_array[boundary_loc_2] = v;
	//
	index_in_v_ascending_degree_array[boundary_v_2] = boundary_loc_1;
	index_in_v_ascending_degree_array[boundary_v_1] = vPtr;
	index_in_v_ascending_degree_array[v] = boundary_loc_2;
	//
	ptr_to_v_gain[score[v]]++;
	ptr_to_v_loss[score[v]]++;
	//
	score[v]--;	
//cout << "gain--. index of " << v << " in ascending degree array is " << index_in_v_ascending_degree_array[v] << endl;
//cout << "score of " << v << " is " << score[v] << endl;
	while(ptr_to_v_gain[max_gain] == ptr_to_v_loss[max_gain + 1])
		max_gain--;
#ifdef debug_mode
//cout << "4. in myBuckets.h, at location 0, the element is " << v_ascending_degree_array[0] << endl;
if(v_ascending_degree_array[0] > 6100)
{
	cout << "4. in myBuckets.h, containing abnormal element, at location 0, the element is " << v_ascending_degree_array[0] << endl;
	//exit(1);
}
#endif
}

void lossPlusPlus(int v, int* score)
{
	int vPtr = index_in_v_ascending_degree_array[v];
	int boundary_loc_1 = ptr_to_v_gain[-score[v]] - 1;
	int boundary_v_1 = v_ascending_degree_array[boundary_loc_1];
	int boundary_loc_2 = ptr_to_v_loss[-score[v]+1] - 1;
	int boundary_v_2 = v_ascending_degree_array[boundary_loc_2];
#ifdef debug_mode
if(boundary_loc_1 < 0)
{
	cout << "4. boundary_loc_1 < 0" << endl;
	exit(1);
}
if(boundary_loc_2 < 0)
{
	cout << "4. boundary_loc_2 < 0" << endl;
	exit(1);
}
cout << "v: " << v << endl;
cout << "vPtr: " << vPtr << endl;
cout << "boundary_v_1: " << boundary_v_1 << endl;
cout << "boundary_loc_1: " << boundary_loc_1 << endl;
cout << "boundary_v_2: " << boundary_v_2 << endl;
cout << "boundary_loc_2: " << boundary_loc_2 << endl;
#endif

	//
	v_ascending_degree_array[vPtr] = boundary_v_1;
	v_ascending_degree_array[boundary_loc_1] = boundary_v_2;
	v_ascending_degree_array[boundary_loc_2] = v;
	//
	index_in_v_ascending_degree_array[boundary_v_2] = boundary_loc_1;
	index_in_v_ascending_degree_array[boundary_v_1] = vPtr;
	index_in_v_ascending_degree_array[v] = boundary_loc_2;
	//
	ptr_to_v_gain[-score[v]]--;
	ptr_to_v_loss[-score[v]+1]--;
	//
	score[v]--;
#ifdef debug_mode
cout << "v: " << v << endl;
cout << "vPtr: " << vPtr << endl;
cout << "boundary_v_1: " << boundary_v_1 << endl;
cout << "boundary_loc_1: " << boundary_loc_1 << endl;
cout << "boundary_v_2: " << boundary_v_2 << endl;
cout << "boundary_loc_2: " << boundary_loc_2 << endl;
#endif
}

int randMinLossVertexInDominatingSet()
{
	int i;
	for(i = 0; ptr_to_v_loss[i] == ptr_to_v_gain[i]; i++){}
	int randPtr = ptr_to_v_loss[i] + rand() % (ptr_to_v_gain[i] - ptr_to_v_loss[i]);
#ifdef debug_mode
//cout << "in myBuckets.h, the rand_min_loc: " << randPtr << endl;
//cout << "in myBuckets.h, rand_min_loss_v: " << v_ascending_degree_array[randPtr] << endl;
//cout << "in myBuckets.h, the loss is " << i << endl;
#endif
	return v_ascending_degree_array[randPtr];
}

int ageMinLossVertexInCover(HugeInt *time_stamp)
{
	int i, j;
	for(i = 0; ptr_to_v_loss[i] == ptr_to_v_gain[i]; i++){}

	int best_v = v_ascending_degree_array[ptr_to_v_loss[i]];
	for(j = ptr_to_v_loss[i] + 1; j < ptr_to_v_gain[i]; j++)
	{
		int v = v_ascending_degree_array[j];
		if(time_stamp[best_v] > time_stamp[v]) best_v = v;
	}

	return best_v;
}

int randMinPositiveGainVertex()
{
	int i;
	for(i = 1; ptr_to_v_gain[i] == ptr_to_v_loss[i+1]; i++){}
	int randPtr = ptr_to_v_gain[i] + rand() % (ptr_to_v_loss[i+1] - ptr_to_v_gain[i]);
	return v_ascending_degree_array[randPtr];
}

int checkMap(int v_num)
//very time-consuming if it exists in frequently called functions or loops
{
#ifdef debug_mode
cout << "checking map in myBuckets.h" << endl;
#endif
	int v;
	for(v = 1; v <= v_num; ++v)
	{
		if(v_ascending_degree_array[index_in_v_ascending_degree_array[v]] != v)
		{
			cout << "vertex " << v << " error!" << endl;
			cout << "location " << index_in_v_ascending_degree_array[v] << " error!" << endl;
			return 0;
		}
	}
#ifdef debug_mode
cout << "having checked map in myBuckets.h" << endl;
#endif
	return 1;
}

int checkPartition(int v_num, int* score)
//very time-consuming if it exists in frequently called functions or loops
{
	int v;

	int i;
	cout << "vertices: " << endl;
	for(i = 0; i < v_num; i++)
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

	for(v = 1; v <= v_num; ++v)
	{
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
	return 1;
}

};
#endif
