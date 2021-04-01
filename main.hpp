#ifndef MAIN_H
#define MAIN_H

#include <vector>
#include <string>
#include "helper/array2d.hpp"
#include "helper/graph.hpp"
#include "helper/util.hpp"

void get_time_MCQD(Array2d<bool>, float, int *, float *, int, float);
float get_best_Tlimit(Array2d<bool>);

Graph array2d_to_adj_list(Array2d<bool>);
Array2d<bool> adj_list_to_array2d(Graph);
void save_data(std::string, std::vector<Graph>, std::vector<float>);
void load_data(std::string, std::vector<Graph>, std::vector<float>);
void generate_train_data(std::string, int, int, int);
void generate_test_data(std::string);
void generate_train_data_from_csv(std::string, std::string);
void get_times(std::vector<Graph>, std::vector<int> *, std::vector<float> *);
void get_times(std::vector<Graph>, std::vector<float>, std::vector<int> *, std::vector<float> *);

// weighted



int main();


#endif // MAIN_H