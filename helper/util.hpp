#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>

#include "helper/array2d.hpp"
#include "graph.hpp"

void array2d_to_adj_list(Array2d<bool> g, Graph *g_new);
void adj_list_to_array2d(Graph g, Array2d<bool> *g_new);
void save_data(std::string path, std::vector<Graph> graphs, std::vector<float> Tlimits);
void load_data(std::string path, std::vector<Graph> *graphs, std::vector<float> *Tlimits);
void load_data(std::string path, std::vector<Graph> *graphs, std::vector<float> *Tlimits, bool shuffle);
void random_graph(Array2d<bool> *g, int n, float p);
void print_array2d(Array2d<bool> a);
void load_predictions(std::string path, std::vector<float> *predictions);
void save_results(std::string path, std::vector<int> steps, std::vector<float> times);
void save_dist(std::string path, std::vector<float> Tlimits, std::vector<float> times);
void save_vector(std::string path, std::vector<float> v);
std::vector<float> make_vector(float beg, float end, float step);
std::vector<float> make_vector(int size, float val);

void print_vector(std::vector<int> v);
std::vector<int> make_vector(int beg, int end, int step);
std::map<int,int> get_shuffled_vertices(int n);

void load_data_weighted(std::string path, std::vector<Graph> *graphs, std::vector< std::vector<double> > *W, std::vector<float> *Tlimits, bool shuffle);

#endif