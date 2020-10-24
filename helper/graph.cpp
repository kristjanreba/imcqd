#include "graph.hpp"

#include <iostream>
#include <cassert>

Graph::Graph() {}

Graph::Graph(const int num_nodes) : num_nodes(num_nodes), num_edges(0), adj_list(num_nodes) {}

void Graph::add_edge(int a, int b) {
    assert(0 <= a && a < num_nodes);
    assert(0 <= b && b < num_nodes);
    if (a > b) std::swap(a, b);
    adj_list[a].push_back(b);
    adj_list[b].push_back(a);
    edge_list.emplace_back(a, b);
    num_edges++;
}

void Graph::print() {
    std::cout << num_nodes << " " << num_edges << std::endl;
    for (int i = 0; i < num_edges; i++) {
        std::cout << edge_list[i].first << " " << edge_list[i].second << std::endl;
    }
}
