/*
    Copyright 2007-2012 Janez Konc

    If you use this program, please cite:
    Janez Konc and Dusanka Janezic. An improved branch and bound algorithm for
   the maximum clique problem. MATCH Commun. Math. Comput. Chem., 2007, 58,
   569-590.

    More information at: http://www.sicmm.org/~konc

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MCQD_H
#define MCQD_H

#include "helper/array2d.hpp"
#include "helper/error.hpp"
#include <algorithm>
#include <assert.h>
#include <exception>
#include <functional>
#include <iostream>
#include <queue>
#include <vector>
#include <string>

using namespace std;

class MaxClique {
public:
  enum SearchType { MCQ, MCQD, MCQW, MCQDW };
  enum HowToSort {
    DESC_DEGREE,
    DESC_ENERGY,
    INCR_ENERGY,
    INCR_DEGREE_INCR_ENERGY,
    DESC_DEGREE_DESC_ENERGY,
    DESC_DEGREE_INCR_ENERGY
  };

  class Clique {
  public:
    typedef int VertexType;
    vector<VertexType> vertices;
    double energy = 0.0;

  public:
    double get_energy() const { return energy; }
    void set_energy(const double ene) { energy = ene; }
    const vector<VertexType> &get_vertices() const { return vertices; }
    void increase_energy(const double e) { energy += e; }
    void decrease_energy(const double e) { energy -= e; }
    void push(const int vertex) { vertices.push_back(vertex); }
    void pop() { vertices.pop_back(); }
    int size() const { return vertices.size(); }
    void print(const vector<double> & = vector<double>());
  };

private:
  static constexpr double BIG_NUM = 1000000.0;
  const Array2d<bool> &__conn;
  const HowToSort __how_to_sort;
  vector<double> __energy;
  const int __max_steps;
  const float __max_time;
  const float Tlimit;
  int pk = 0, level = 1;
  float t; // current time of execution

  class Vertices {
    class Vertex {
      int i, d;
      double e;

    public:
      void set_i(const int ii) { i = ii; }
      int get_i() const { return i; }
      void set_degree(int dd) { d = dd; }
      int get_degree() const { return d; }
      void set_energy(double ee) { e = ee; }
      double get_energy() const { return e; }
    };
    Vertex *v;
    int sz;

    static bool desc_degree(const Vertex &vi, const Vertex &vj) {
      return (vi.get_degree() > vj.get_degree());
    }
    static bool desc_energy(const Vertex &vi, const Vertex &vj,
                            const vector<double> &energy) {
      return energy[vi.get_i()] > energy[vj.get_i()];
    }
    static bool incr_energy(const Vertex &vi, const Vertex &vj,
                            const vector<double> &energy) {
      return energy[vi.get_i()] < energy[vj.get_i()];
    }
    static bool incr_degree_incr_energy(const Vertex &vi, const Vertex &vj,
                                        const vector<double> &energy) {
      return vi.get_degree() < vj.get_degree() ||
             (vi.get_degree() == vj.get_degree() &&
                 energy[vi.get_i()] < energy[vj.get_i()]);
    }
    static bool desc_degree_desc_energy(const Vertex &vi, const Vertex &vj,
                                        const vector<double> &energy) {
      return vi.get_degree() > vj.get_degree() ||
             (vi.get_degree() == vj.get_degree() &&
                 energy[vi.get_i()] > energy[vj.get_i()]);
    }
    static bool desc_degree_incr_energy(const Vertex &vi, const Vertex &vj,
                                        const vector<double> &energy) {
      return vi.get_degree() > vj.get_degree() ||
             (vi.get_degree() == vj.get_degree() &&
                 energy[vi.get_i()] < energy[vj.get_i()]);
    }

  public:
    Vertices(int size) : sz(0) { v = new Vertex[size]; }
    ~Vertices() {}
    void dispose() {
      if (v)
        delete[] v;
    }
    void sort(MaxClique::HowToSort, const vector<double> & = vector<double>());
    void init_colors();
    void set_degrees(MaxClique &);
    int size() const { return sz; }
    void push(const int ii) { v[sz++].set_i(ii); }
    void pop() { sz--; }
    Vertex &at(const int ii) const { return v[ii]; }
    Vertex &end() const { return v[sz - 1]; }
    void print(const vector<double> & = vector<double>());
  };

  class ColorClass {
    int *i;
    int sz;
    double min_energy;

  public:
    ColorClass() : sz(0), i(0), min_energy(BIG_NUM) {}
    ColorClass(const int sz) = delete;
    ~ColorClass() {
      if (i)
        delete[] i;
    }
    ColorClass(const ColorClass &) = delete;
    void init(const int sz) {
      i = new int[sz];
      rewind();
    }
    void push(const int ii) { i[sz++] = ii; };
    void pop() { sz--; };
    void rewind() {
      sz = 0;
      min_energy = BIG_NUM;
    };
    int size() const { return sz; }
    int &at(const int ii) const { return i[ii]; }
    void set_min_energy(const double ee) {
      if (ee < min_energy)
        min_energy = ee;
    }
    double get_min_energy() const { return min_energy; }
  };

  Vertices V;
  ColorClass *C;

  class BestCliques {
  private:
    struct incr_energy {
      bool operator()(const Clique &lhs, const Clique &rhs) const {
        return lhs.get_energy() < rhs.get_energy();
      }
    };
    typedef priority_queue<Clique, vector<Clique>, incr_energy> Queue;
    const int maxsz;
    Queue QMAXES;

  public:
    BestCliques(const int maxsz) : maxsz(maxsz) {}
    void insert_if_better(const Clique &clique);
    void insert(const Clique &clique) { QMAXES.push(clique); }
    vector<Clique> get_k_cliques() const;
    double get_highest_energy() const {
      return QMAXES.empty() ? BIG_NUM : QMAXES.top().get_energy();
    }
  };

  BestCliques __best_cliques;
  Clique Q, QMAX;

  class StepCount {
    int i1, i2;

  public:
    StepCount() : i1(0), i2(0) {}
    void set_i1(const int ii) { i1 = ii; }
    int get_i1() const { return i1; }
    void set_i2(const int ii) { i2 = ii; }
    int get_i2() const { return i2; }
    void inc_i1() { i1++; }
  };

  StepCount *S;

  const bool connection(const int i, const int j) const {
    return __conn.get(i, j);
  }
  void init(const MaxClique::SearchType);
  bool cut1(const int, const ColorClass &);
  void cut2(const Vertices &, Vertices &);

  void color_sort_mcq(Vertices &);
  void color_sort_mcqw(Vertices &);
  void color_sort_kcqw(Vertices &);
  double sum_kcqw(const vector<double> &, int);
  void expand_mcq(Vertices);
  void expand_mcqd(Vertices);
  void expand_mcqw(Vertices);
  void expand_mcqdw(Vertices);
  void expand_kcq(Vertices);
  void expand_kcqd(Vertices);
  void expand_kcqw(Vertices);
  void expand_kcqdw(Vertices);
  void degree_sort(Vertices &R) {
    R.set_degrees(*this);
    R.sort(__how_to_sort, __energy);
  }
  bool ub_is_breached(const vector<double> &, const double);
  void check_ub(const vector<double> &);
  double tighten_ub(const vector<double> &, const double);

public:
  MaxClique(const Array2d<bool> &, const HowToSort = DESC_DEGREE,
            const int = 100, const int = -1, const float = 0.025, const float = 5.0);
  MaxClique(const Array2d<bool> &, const vector<double> &,
            const HowToSort = DESC_DEGREE, const int = 100, const int = -1,
            const float = 0.025, const float = 5.0);
  ~MaxClique() {
    if (C)
      delete[] C;
    if (S)
      delete[] S;
    V.dispose();
  };

  int steps() const { return pk; }

  Clique maximum_clique(SearchType); // find one maximum clique (MCQ or MCQD)
  vector<Clique>
  k_cliques(const int,
            const SearchType); // find all k-cliques (larger than minimum size)
                               // with best weights (modified MCQ or MCQW)

  static void output_graph(const Array2d<bool> &, const std::vector<double> &);
  static void read_graph(Array2d<bool> &, std::vector<double> &);
  static void read_graph(Array2d<bool> &, std::string);
  static void print_color_classes(const int);
  static void print_vertices(Vertices R);
};
#endif
