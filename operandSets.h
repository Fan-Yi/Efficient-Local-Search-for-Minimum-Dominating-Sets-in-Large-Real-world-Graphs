#ifndef OPERAND_SETS_H
#define OPERAND_SETS_H

#include "hugeInt.h"

class ScoreAgePickSet : public Bijection
{
public:
	ScoreAgePickSet(int sz) : Bijection(sz)
	{
	}
	~ScoreAgePickSet()
	{
	}
/*
	int best_element(Vertex* vertices, HugeInt* time_stamp)
	{
		if(!array_size) return 0;
		int best_v = array[begin()];
		for(int i = begin() + 1; i != end(); i++)
		{
			int v = array[i];
			if(-vertices[v].get_weight() > -vertices[best_v].get_weight() || (vertices[v].get_weight() == vertices[best_v].get_weight() && time_stamp[v] < time_stamp[best_v]))
				best_v = v;
		}
		return best_v;
	}
*/
	int best_element(long double* score, HugeInt* time_stamp)
	{
		if(!array_size) return 0;
		int best_e = at(begin());
		int i;
		for(i = begin() + 1; i != end(); i++)
		{
			int e = array[i];
			if(score[e] > score[best_e] || (score[e] == score[best_e] && time_stamp[e] < time_stamp[best_e]))
				best_e = e;
		}
		return best_e;
	}
};
class ConfchangedScoreAgePickSet : public Bijection
{
public:	
	ConfchangedScoreAgePickSet(int sz) : Bijection(sz)
	{
	}
	~ConfchangedScoreAgePickSet()
	{
	}

	int best_element(int* confChange, int* score, HugeInt* time_stamp)
	{
		if(!array_size) return 0;
		int best_e = 0;
		int i;
		for(i = begin(); i != end(); i++)
		{
			int e = array[i]; 
			if(confChange[e])
			{
				best_e = e;
				i++;
				break;
			}
		}
		for(; i != end(); i++)
		{
			int e = array[i];
			if(!confChange[e]) continue;
			if(score[e] > score[best_e] || (score[e] == score[best_e] && time_stamp[e] < time_stamp[best_e]))
				best_e = e;
		}
		return best_e;
	}
/*
	int best_element(int* confChange, Vertex* vertices, HugeInt* time_stamp)
	{
		if(!array_size) return 0;
		int best_v = 0;
		int i;
		for(i = begin(); i != end(); i++)
		{
			int v = array[i]; 
			if(confChange[v])
			{
				best_v = v;
				i++;
				break;
			}
		}
		for(; i != end(); i++)
		{
			int v = array[i];
			if(!confChange[v]) break;
			if(vertices[v].get_weight() > vertices[best_v].get_weight() || (vertices[v].get_weight() == vertices[best_v].get_weight() && time_stamp[v] < time_stamp[best_v]))
				best_v = v;
		}
		return best_v;
	}
*/
};
#endif
