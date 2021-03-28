#include "main.hpp"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>

#include "graph/mcqd.hpp"
#include "helper/benchmark.hpp"
#include "helper/array2d.hpp"
#include "helper/graph.hpp"
#include "helper/util.hpp"

#define TLIMIT_DEFAULT 0.025

void get_time_MCQD(Array2d<bool> g, float Tlimit, int *steps, float *time, int max_steps, float max_time) {
    Benchmark benchmark;
    MaxClique::SearchType search = MaxClique::MCQD;
    MaxClique::HowToSort sort = MaxClique::DESC_DEGREE;
    int n_cliques = 100;
    MaxClique mc(g, sort, n_cliques, max_steps, Tlimit, max_time);
    benchmark.reset();
    MaxClique::Clique c = mc.maximum_clique(search);
    *time = benchmark.seconds_from_start();
    *steps = mc.steps();
    //std::cout << *time << std::endl;
    if (fabs(*time - max_time) < 0.1) {
        *time = -1;
        //std::cout << "set to -1" << std::endl;
    }
}

float get_best_Tlimit(Array2d<bool> g) {
    float Tlimits[] = {0.00001, 0.0001, 0.001, 0.002, 0.0025, 0.005, 0.01, 0.02, 0.025, 0.03, 0.035, 0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1.0};
    float best_Tlimit = -1;
    float best_time = -1;
    for (int j = 0; j < sizeof(Tlimits)/sizeof(Tlimits[0]); j++) { // test multiple Tlimit values
        float time;
        int steps;
        int max_steps = -1;
        float max_time = 5.0;
        get_time_MCQD(g, Tlimits[j], &steps, &time, max_steps, max_time);
        if (fabs(time - max_time) < 0.1) continue;
        if (time < best_time || best_time == -1) {
            best_time = time;
            best_Tlimit = Tlimits[j];
        }
    }
    return best_Tlimit;
}

void get_time_dist(std::string path_dist, Array2d<bool> g){
    std::vector<float> Tlimits{0.0, 0.00001, 0.0001, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.007, 0.01, 0.02, 0.025, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0};
    std::vector<float> times;
    for (int j = 0; j < Tlimits.size(); j++) { // test multiple Tlimit values
        float time;
        int steps;
        int max_steps = -1;
        float max_time = 60.0;
        get_time_MCQD(g, Tlimits[j], &steps, &time, max_steps, max_time);
        if (fabs(time - max_time) < 0.1) {
            time = -1.;
        }
        times.push_back(time);
    }
    save_dist(path_dist, Tlimits, times);
}

void generate_train_data(std::string path, int dataset_size, int max_vertices, int max_edges) {
    std::vector<Graph> graphs;
    std::vector<float> Tlimits;

    std::random_device rd;  // will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (int i = 0; i < dataset_size; i++) {
        int n = (int)(dis(gen) * max_vertices + 1); // number of vertices
        float p = dis(gen); // edge probability
        if ((n*n*p) > max_edges) {
            i--;
            continue;
        }
        Array2d<bool> g(n);
        random_graph(&g, n, p);
        float best_Tlimit = get_best_Tlimit(g);
        Graph g_adj(n);
        array2d_to_adj_list(g, &g_adj);
        graphs.push_back(g_adj);
        Tlimits.push_back(best_Tlimit);
        
        std::cout << i << std::endl;
    }
    save_data(path, graphs, Tlimits);
}

void generate_test_data(std::string path) {
    int sizes[] = {100, 100, 100, 100, 100, 
                    150, 150, 150, 150, 150, 150,
                    200, 200, 200, 200, 200, 200, 200,
                    300, 300, 300, 300,
                    500, 500, 500, 500,
                    1000, 1000, 1000, 1000};
    float edge_prob[] = {0.6, 0.7, 0.8, 0.9, 0.95,
                            0.5, 0.6, 0.7, 0.8, 0.9, 0.95,
                            0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,
                            0.4, 0.5, 0.6, 0.7,
                            0.3, 0.4, 0.5, 0.6,
                            0.2, 0.3, 0.4, 0.5};
    std::vector<Graph> graphs;
    std::vector<float> Tlimits;
    for (int i = 0; i < sizeof(sizes)/sizeof(sizes[0]); i++) {
        int n = sizes[i]; // number of vertices
        float p = edge_prob[i]; // edge probability
        Array2d<bool> g(n);
        random_graph(&g, n, p);
        float best_Tlimit = 0.0;
        Graph g_adj(n);
        array2d_to_adj_list(g, &g_adj);
        graphs.push_back(g_adj);
        Tlimits.push_back(best_Tlimit);
    }
    save_data(path, graphs, Tlimits);
}

void generate_dense_graphs(std::string path, int dataset_size, int max_vertices, bool calculate_best_Tlimit) {
    std::vector<Graph> graphs;
    std::vector<float> Tlimits;

    std::random_device rd;  // will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (int i = 0; i < dataset_size; i++) {
        int n = (int)(dis(gen) * max_vertices + 1); // number of vertices
        float p = 0.99 + (dis(gen)/100.0); // edge probability
        std::cout << p << std::endl;
        Array2d<bool> g(n);
        random_graph(&g, n, p);
        float best_Tlimit = 0;
        if (calculate_best_Tlimit) best_Tlimit = get_best_Tlimit(g);
        Graph g_adj(n);
        array2d_to_adj_list(g, &g_adj);
        graphs.push_back(g_adj);
        Tlimits.push_back(best_Tlimit);
    }
    save_data(path, graphs, Tlimits);
}

void generate_sparse_graphs(std::string path, int dataset_size, int max_vertices) {
    std::vector<Graph> graphs;
    std::vector<float> Tlimits;

    std::random_device rd;  // will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 0.3);
    
    for (int i = 0; i < dataset_size; i++) {
        int n = (int)(dis(gen) * max_vertices + 1); // number of vertices
        float p = dis(gen); // edge probability
        std::cout << p << std::endl;
        Array2d<bool> g(n);
        random_graph(&g, n, p);
        float best_Tlimit = 0.0;
        Graph g_adj(n);
        array2d_to_adj_list(g, &g_adj);
        graphs.push_back(g_adj);
        Tlimits.push_back(best_Tlimit);
    }
    save_data(path, graphs, Tlimits);
}

void generate_train_data_from_csv(std::string path_in, std::string path_out) {
    std::vector<Graph> graphs;
    std::vector<float> Tlimits;
    std::vector<float> best_Tlimits;
    std::cout << "Loading data from " << path_in << std::endl;
    load_data(path_in, &graphs, &Tlimits);
    int n = graphs.size();
    for (int i = 0; i < n; i++) {
        std::cout << i << std::endl;
        Array2d<bool> g_new(graphs[i].num_nodes);
        adj_list_to_array2d(graphs[i], &g_new);
        float tlimit = get_best_Tlimit(g_new);
        best_Tlimits.push_back(tlimit);
    }
    save_vector(path_out, best_Tlimits);
}

void get_times(std::vector<Graph> graphs, std::vector<float> predictions, std::vector<int> *steps, std::vector<float> *times) {
    if (graphs.size() != predictions.size()) {
        std::cout << "graphs and predictions are not the same size" << std::endl;
        return;
    }
    int n = graphs.size();
    //int n = 5;
    for (int i = 0; i < n; i++) {
        std::cout << i << std::endl;
        Array2d<bool> g_new(graphs[i].num_nodes);
        adj_list_to_array2d(graphs[i], &g_new);
        float time;
        int s;
        int max_steps = -1;
        float max_time = 2000.0;
        get_time_MCQD(g_new, predictions[i], &s, &time, max_steps, max_time);
        times->push_back(time);
        steps->push_back(s);
    }
}

void get_times(std::vector<Graph> graphs, std::vector<int> *steps, std::vector<float> *times) {
    std::vector<float> predictions(graphs.size(), TLIMIT_DEFAULT);
    get_times(graphs, predictions, steps, times);
}

void test_model_on_dataset(std::string dataset, std::string model_name) {
    std::string path_predictions = "pred/" + dataset + "_" + model_name + ".csv";
    std::string path_results_pred = "results/" + dataset + "_" + model_name + ".csv";
    std::string path_data = "datasets/" + dataset + "_test.csv";

    std::vector<float> predictions;
    std::vector<float> times_pred;
    std::vector<int> steps_pred;
    std::vector<Graph> graphs;
    std::vector<float> _Tlimits; // we don't need this actually but the function requires this argumet

    std::cout << "Testing " << model_name << " on " << dataset << " data" << std::endl;
    load_data(path_data, &graphs, &_Tlimits); // load graphs
    if (model_name.compare("default") == 0) predictions = make_vector(graphs.size(), TLIMIT_DEFAULT);
    else load_predictions(path_predictions, &predictions); // load predictions
    get_times(graphs, predictions, &steps_pred, &times_pred);
    save_results(path_results_pred, steps_pred, times_pred);
}

int main()
{
    // Novi eksperimenti 26.3.2021
    /*
    int n = 5;
    std::map<int, int> d = get_shuffled_vertices(n);
    std::cout << "Test 3" << std::endl;
    for (auto itr = d.find(3); itr != d.end(); itr++) {
        cout << itr->first << '\t' << itr->second << '\n';
    }
    */

    



    // -----------------------------



    // Stari eksperimenti
    //generate_train_data_from_csv("datasets/docking_train.csv", "datasets/docking_train_tlimits.csv");
    //generate_train_data_from_csv("datasets/product_train.csv", "datasets/product_train_tlimits.csv");

    //std::vector<std::string> test_datasets = {"protein", "product", "docking", "rand", "dense"};
    //std::vector<std::string> models = {"default", "gcn", "xgb", "svr", "gbr"};
    /*
    std::vector<std::string> models = {"gcn", "gin", "gat"};
    for (int i = 0; i < models.size(); i++) {
        for (int j = 0; j < test_datasets.size(); j++) {
            test_model_on_dataset(test_datasets[j], models[i]);
        }
    }
    */

    

    //test_model_on_dataset("rand", "default");
    //test_model_on_dataset("dimacs", "default");


    /*
    int i = 2;
    std::vector<Graph> graphs;
    std::vector<float> _Tlimits;
    load_data("datasets/sparse_data.csv", &graphs, &_Tlimits); // load graphs
    std::cout << graphs[i].num_nodes << ' ' << graphs[i].num_edges << std::endl;
    std::cout << (float)graphs[i].num_edges / (float)(graphs[i].num_nodes*(graphs[i].num_nodes-1)/2.0) << std::endl;
    Array2d<bool> g_new(graphs[i].num_nodes);
    adj_list_to_array2d(graphs[i], &g_new);
    get_time_dist("data/dist.csv", g_new);
    */

    // models: gnn, default, svm, adaboost, gin, ...
    // datasets: rand, dimacs, train

    //test_model_on_dataset("rand", "default");
    //test_model_on_dataset("rand", "gnn");
    //test_model_on_dataset("rand", "gin");
    //test_model_on_dataset("rand", "adaboost");
    //test_model_on_dataset("rand", "gbr");
    //test_model_on_dataset("rand", "svr");

    //test_model_on_dataset("dimacs", "default");
    //test_model_on_dataset("dimacs", "gnn");
    //test_model_on_dataset("dimacs", "adaboost");
    //test_model_on_dataset("dimacs", "gbr");
    //test_model_on_dataset("dimacs", "svr");


    
    // ------------- LEGACY CODE --------------
    
    /*
    // Generate dense data
    string path = "data/dense_data.csv";
    int dataset_size = 20;
    int max_vertices = 1000;
    generate_dense_graphs(path, dataset_size, max_vertices);
    std::cout << "Done creating dataset." << std::endl;
    */
    
    /*
    // Generate train data
    string path = "data/train_data_4.csv";
    int dataset_size = 100'000;
    int max_vertices = 2000;
    int max_edges = 20'000;
    generate_train_data(path, dataset_size, max_vertices, max_edges);
    std::cout << "Done creating dataset." << std::endl;
    */

    
    // Generate random test data
    //string path = "data/rand_data.csv";
    //generate_test_data(path);
    //std::cout << "Done creating test dataset." << std::endl;
    

    /*
    // Load data
    std::vector<Graph> graphs;
    std::vector<float> Tlimits;
    load_data(path, &graphs, &Tlimits);
    Array2d<bool> g_load(graphs[0].num_nodes);
    adj_list_to_array2d(graphs[0], &g_load);
    cout << "Done loading dataset." << endl;
    */
    
    /*
    string graph_name = "johnson8-2-4.clq";
    string path = folder + graph_name;

    // parameters for MCQD
    Array2d<bool> g;
    MaxClique::SearchType search = MaxClique::MCQD;
    MaxClique::HowToSort sort = MaxClique::DESC_DEGREE;
    int n_cliques = 100;
    int max_steps = -1;
    float Tlimit = 0.025;

    // read graph and find max clique
    MaxClique::read_graph(g, path);
    const Array2d<bool> conn = g;
    MaxClique mc(conn, sort, n_cliques, max_steps, Tlimit);
    benchmark.reset();
    MaxClique::Clique c = mc.maximum_clique(search);
    double time = benchmark.seconds_from_start();
    c.print();
    cout << "Time to find maximum clique: " << time << " s" << endl;
    */

    cout << "Got to the end." << endl;
    return 0;
} 
