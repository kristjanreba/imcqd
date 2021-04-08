#include "util.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <map>

#include "helper/array2d.hpp"
#include "graph.hpp"


void array2d_to_adj_list(Array2d<bool> g, Graph *g_new) {
    for (int i = 0; i < g.get_szi(); i++) {
        for (int j = 0; j < i; j++) {
            if (g.get(i, j) == 1) {
                g_new->add_edge(i, j);
            }
        }
    }
}

void adj_list_to_array2d(Graph g, Array2d<bool> *g_new) {
    for (auto& p : g.edge_list) {
        int u = p.first, v = p.second;
        g_new->set(v, u);
        g_new->set(u, v);
    }
}

// -------------- Loading data ------------------

void load_data(std::string path, std::vector<Graph> *graphs, std::vector<float> *Tlimits) {
    std::ifstream f(path);
    int n; // number of graphs
    f >> n;

    // shuffle vertices
    std::map<int, int> d = get_shuffled_vertices(n);

    for (int i = 0; i < n; i++) {
        int n_vertices, m; // number of vertices in graph
        f >> n_vertices >> m;
        Graph g(n_vertices);
        for (int j = 0; j < m; j++) {
            int u, v;
            f >> u >> v;
            g.add_edge(d.find(u)->second, d.find(v)->second);
        }
        graphs->push_back(g);
        float Tlimit;
        f >> Tlimit;
        Tlimits->push_back(Tlimit);
    }
    f.close();
}

void load_predictions(std::string path, std::vector<float> *predictions) {
    std::ifstream f(path);
    int n; // number of predictions
    f >> n;
    for (int i = 0; i < n; i++) {
        float p;
        f >> p;
        predictions->push_back(p);
    }
    f.close();
}

void random_graph(Array2d<bool> *g, int n, float p) {
    std::random_device rd;  // will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    // init edges so that the percantage of all possible edges is equal to p.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            if (dis(gen) <= p) {
                g->set(i,j);
                g->set(j,i);
            }
        }
    }
}

void print_array2d(Array2d<bool> a) {
    for (int i = 0; i < a.get_szi(); i++) {
        for (int j = 0; j < a.get_szj(); j++) {
            std::cout << a.get(i, j);
        }
        std:cout << std::endl;
    }
}

void print_vector(std::vector<int> v) {
    for (int i = 0; i < v.size(); i++) {
        std::cout << v[i] << " ";
    }
    std:cout << std::endl;
}

// -------------- Saving data ------------------

void save_data(std::string path, std::vector<Graph> graphs, std::vector<float> Tlimits) {
    if (graphs.size() != Tlimits.size()) return;
    int n = graphs.size();
    std::ofstream f(path);
    f << n << '\n';
    for (int i = 0; i < n; i++) {
        f << graphs[i].num_nodes << ' ' << graphs[i].num_edges << '\n';
        for (auto& p : graphs[i].edge_list) {
            int u = p.first, v = p.second;
            f << u << ' ' << v << '\n';
        }
        f << Tlimits[i] << '\n';
    }
    f.close();
}

void save_results(std::string path, std::vector<int> steps, std::vector<float> times) {
    int n = times.size();
    std::ofstream f(path);
    f << n << '\n';
    for (int i = 0; i < n; i++) {
        f << steps[i] << ' ' << times[i] << '\n';
    }
    f.close();
}

void save_dist(std::string path, std::vector<float> Tlimits, std::vector<float> times) {
    int n = times.size();
    std::ofstream f(path);
    f << n << '\n';
    for (int i = 0; i < n; i++) {
        f << Tlimits[i] << ' ' << times[i] << '\n';
    }
    f.close();
}

void save_vector(std::string path, std::vector<float> v) {
    int n = v.size();
    std::ofstream f(path);
    f << n << '\n';
    for (int i = 0; i < n; i++) {
        f << v[i] << '\n';
    }
    f.close();
}

// vectors

std::vector<int> make_vector(int beg, int end, int step) {
    std::vector<int> vec;
    vec.reserve((end - beg) / step + 1);
    while (beg <= end) {
        vec.push_back(beg);
        beg += step;
    }
    return vec;
}

std::vector<float> make_vector(int size, float val) {
    std::vector<float> vec;
    for (int i = 0; i < size; i++)
        vec.push_back(val);
    return vec;
}

std::vector<double> make_vector(int size, double val) {
    std::vector<double> vec;
    for (int i = 0; i < size; i++)
        vec.push_back(val);
    return vec;
}

// create shuffled dict
std::map<int, int> get_shuffled_vertices(int n) {
    std::map<int, int> d;
    std::vector<int> vertices = make_vector(0, n-1, 1);
    //print_vector(vertices);
    auto rd = std::random_device {}; 
    auto rng = std::default_random_engine { rd() };
    std::shuffle(std::begin(vertices), std::end(vertices), rng);
    //print_vector(vertices);
    for (int i = 0; i < n; i++) {
        d.insert({i, (int) vertices[i]});
    }
    return d;
}


// weighted graphs

void load_data_weighted(std::string path, std::vector<Graph> *graphs, std::vector< std::vector<double> > *W, std::vector<float> *Tlimits, bool shuffle) {
    std::ifstream f(path);
    int n; // number of graphs
    f >> n;
    
    for (int i = 0; i < n; i++) {
        int n_vertices, m; // number of vertices in graph
        f >> n_vertices >> m;

        // shuffle vertices
        std::map<int, int> d;
        if (shuffle) { d = get_shuffled_vertices(n_vertices); }

        std::vector<double> ws = make_vector(n_vertices, 0.0);
        for (int k = 0; k < n_vertices; k++) {
            double w;
            f >> w;
            ws[d.find(k)->second] = w;
        }
        W->push_back(ws);

        Graph g(n_vertices);
        for (int j = 0; j < m; j++) {
            int u, v;
            f >> u >> v;
            //std::cout << u << " " << v << " " << n_vertices << std::endl;
            if (shuffle) {
                g.add_edge(d.find(u)->second, d.find(v)->second);
            } else {
                g.add_edge(u, v);
            }
        }
        graphs->push_back(g);
        float Tlimit;
        f >> Tlimit;
        Tlimits->push_back(Tlimit);
    }
    f.close();
}


