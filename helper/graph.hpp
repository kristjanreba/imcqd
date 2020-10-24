#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <map>

class Graph {
public:
    Graph();
    Graph(const int num_nodes);
    int num_nodes;
    int num_edges;
    std::vector<std::vector<int>> adj_list;
    std::vector<std::pair<int, int>> edge_list;

    void add_edge(int a, int b);
    void print();
};

#endif
