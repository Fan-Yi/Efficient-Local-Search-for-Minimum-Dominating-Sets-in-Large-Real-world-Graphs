#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <string.h>
#include <unordered_set>
#include <math.h>

//#define test_graph_construction_mode

using namespace std;

typedef long long LL;

class Vertex
{
private:
	int *neighbors;
	int *adj_edges;
	int degree;
	int weight;
public:
	Vertex()
	{
		neighbors = NULL;
		adj_edges = NULL;
		degree = 0;
		weight = 0;
	}
	~Vertex()
	{
		free_neighborhood_space();
	}
	void allocate_neighborhood_space(int dgr)
	{
		neighbors = new int[dgr];
		adj_edges = new int[dgr];
	}
	void free_neighborhood_space()
	{
		delete[] neighbors;
		delete[] adj_edges;
	}
	void add_neighbor(int name, int index)
	{
		neighbors[index] = name;
	}
	void add_adj_edge(int name, int index)
	{
		adj_edges[index] = name;
	}
	int *get_neighbors()
	{
		return neighbors;
	}
	int *get_adj_edges()
	{
		return adj_edges;
	}
	void set_degree(int d)
	{
		degree = d;
	}
	int get_degree()
	{
		return degree;
	}
	int get_weight()
	{
		return weight;
	}
	void set_weight(int w)
	{
		weight = w;
	}
	void show_neighbors()
	{
		cout << "neighbors: ";
		for(int i = 0; i < degree; i++)
		{
			cout << neighbors[i] << '\t';
		}
		cout << endl;
	}
};

class Edge
{
private:
	int v1, v2;
public:
	Edge(){}
	void set_vertices(int u1, int u2)
	{
		v1 = u1;
		v2 = u2;
	}
	void get_vertices(int &u1, int &u2)
	{
		u1 = v1;
		u2 = v2;
	}
	~Edge(){}
};

class Graph
{
private:
	int max_degree;
protected:
	Vertex *vertices;
	Edge	*edges;
	int v_num;
	int e_num;
	unordered_set<LL> edge_hash_id_set;
private:
	void encode_pairID(LL &pairID, LL n1, LL n2)
	{
		pairID = ((n1 + n2 + 1) * (n1 + n2) >> 1) + n2;
	}
	void decode_pairID(LL pairID, LL &n1, LL &n2)
	{
		LL w = LL((sqrt(double((pairID << 3) + 1)) - 1) / 2);
		LL t = (w * w + w) >> 1;
		n2 = pairID - t;
		n1 = w - n2; 
	}
	void encode_unordered_pairID(LL &pairID, LL n1, LL n2)
	{
		LL u, v;
		if(n1 < n2)
		{
			u = n1; v = n2;
		}
		else
		{
			u = n2; v = n1;
		}
		encode_pairID(pairID, u, v);		
	}
	void insertEdgeHashIDToSet(int n1, int n2)
	{
		LL edge_hash_id;
		encode_unordered_pairID(edge_hash_id, LL(n1), LL(n2));
		edge_hash_id_set.insert(edge_hash_id);
	}
#ifdef test_graph_construction_mode
void output_info_on_vertex(int u)
{
cout << "the degree of " << u << " is " << vertices[u].get_degree() << endl;
cout << "its neighbor list: \n";
for(int i = 0; i < vertices[u].get_degree(); i++)
{
	cout << vertices[u].get_neighbors()[i] << '\n';
}
cout << "its incident edge list: \n";
for(int i = 0; i < vertices[u].get_degree(); i++)
{
	cout << vertices[u].get_adj_edges()[i] << '\n';
}
}
#endif
public:
	Graph(char *filename)
	{
		ifstream infile(filename);
		if(infile == NULL)
		{
			cout << "File " << filename << " cannot be opened" << endl;
			exit(1);
		}

		char line[1024];
		infile.getline(line, 1024);

		while(line[0] != 'p')
			infile.getline(line, 1024);

		char tempstr1[1024], tempstr2[1024];
		sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &v_num, &e_num);

		if(strcmp(tempstr1, "p") != 0 || strcmp(tempstr2, "edge") != 0)
		{
			cout << "format error occurs in reading p lines" << endl;
			exit(1);
		}

		vertices = new Vertex[v_num + 2];
		edges = new Edge[e_num + 2];

		char ch_tmp;
		int v;

		int v1, v2;
		int *v_degree_tmp = new int[v_num + 2];
		memset(v_degree_tmp, 0, sizeof(int) * (v_num + 2));

		int e;
		for(e = 0; e < e_num; e++)
		{
			infile >> ch_tmp >> v1 >> v2;
			edges[e].set_vertices(v1, v2);
			v_degree_tmp[v1]++;
			v_degree_tmp[v2]++;
			insertEdgeHashIDToSet(v1, v2);
		}

		for(v = 1; v <= v_num; v++)
		{
			vertices[v].allocate_neighborhood_space(v_degree_tmp[v]);
			vertices[v].set_degree(v_degree_tmp[v]);
		}
		max_degree = v_degree_tmp[1];
		for(int i = 2; i <= v_num; i++)
		{
			if(v_degree_tmp[i] > max_degree)
				max_degree = v_degree_tmp[i];
		}

		memset(v_degree_tmp, 0, sizeof(int) * (v_num + 2));
		for(e = 0; e < e_num; e++)
		{
			edges[e].get_vertices(v1, v2);
			vertices[v1].add_neighbor(v2, v_degree_tmp[v1]);
			vertices[v2].add_neighbor(v1, v_degree_tmp[v2]);
			vertices[v1].add_adj_edge(e, v_degree_tmp[v1]);
			vertices[v2].add_adj_edge(e, v_degree_tmp[v2]);
			v_degree_tmp[v1]++;
			v_degree_tmp[v2]++;
		}

		delete[] v_degree_tmp;
		infile.close();

#ifdef test_graph_construction_mode
int u;
u = 1;
output_info_on_vertex(u);
u = 2;
output_info_on_vertex(u);
u = v_num - 1;
output_info_on_vertex(u);
u = v_num;
output_info_on_vertex(u);
cout << "**************" << endl;
cout << "max degree: " << max_degree << endl;
float d_sum = 0;
float d_avg = 0;
for(u = 1; u <= v_num; u++)
{
	d_sum += float(vertices[u].get_degree());
}
d_avg = d_sum / float(v_num);
cout << "avg degree: " << d_avg << endl;
#endif

	}
	Vertex* get_vertices()
	{
		return vertices;
	}
	int get_max_degree()
	{
		return max_degree;
	}
	int isConnected(int n1, int n2)
	{
		LL edge_hash_id;
		encode_unordered_pairID(edge_hash_id, LL(n1), LL(n2));
		return edge_hash_id_set.count(edge_hash_id);
	}
	~Graph()
	{
		delete[] vertices;
		delete[] edges;
	}
};
#endif
