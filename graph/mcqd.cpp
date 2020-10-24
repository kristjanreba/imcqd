#include "mcqd.hpp"
#include "helper/help.hpp"
#include "helper/inout.hpp"
#include <numeric>
#include <utility>
#include <random>
#include <time.h>
#include <iostream>


/*
 * Initialize maximum clique algorithm.
 */
MaxClique::MaxClique(const Array2d<bool> &conn, const HowToSort h,
                     const int n_cliques, const int max_steps, const float tt, const float max_time)
    : __conn(conn), __how_to_sort(h), __best_cliques(n_cliques),
      __max_steps(max_steps), Tlimit(tt), V(conn.get_szi()), __max_time(max_time) {
  if (conn.get_szi() == 0)
    throw Error("NOTE : Graph is empty (meaning that the molecule is built of "
                "non-seed fragments only)."); // fixes issue #116
  for (int i = 0; i < conn.get_szi(); i++)
    V.push(i);
  C = new ColorClass[conn.get_szi() + 1];
  for (int i = 0; i < conn.get_szi() + 1; i++)
    C[i].init(conn.get_szi() + 1);
  S = new StepCount[conn.get_szi() + 1];
  t = clock();
}

/*
 * Initialize maximum weight clique algorithm, requires energy for each vertex.
 */
MaxClique::MaxClique(const Array2d<bool> &conn, const vector<double> &energy,
                     const HowToSort h, const int n_cliques,
                     const int max_steps, const float tt)
    : MaxClique(conn, h, n_cliques, max_steps, tt) {
  __energy = energy;
}

void MaxClique::Vertices::sort(MaxClique::HowToSort h,
                               const vector<double> &energy) {
  using namespace std::placeholders;
  switch (h) {
  case DESC_DEGREE:
    std::sort(v, v + sz, desc_degree);
    break;
  case DESC_ENERGY:
    std::sort(v, v + sz, std::bind(&Vertices::desc_energy, _1, _2, energy));
    break;
  case INCR_ENERGY:
    std::sort(v, v + sz, std::bind(&Vertices::incr_energy, _1, _2, energy));
    break;
  case INCR_DEGREE_INCR_ENERGY:
    std::sort(v, v + sz,
              std::bind(&Vertices::incr_degree_incr_energy, _1, _2, energy));
    break;
  case DESC_DEGREE_DESC_ENERGY:
    std::sort(v, v + sz,
              std::bind(&Vertices::desc_degree_desc_energy, _1, _2, energy));
    break;
  case DESC_DEGREE_INCR_ENERGY:
    std::sort(v, v + sz,
              std::bind(&Vertices::desc_degree_incr_energy, _1, _2, energy));
    break;
  };
}

void MaxClique::Vertices::init_colors() {
  const int max_degree = v[0].get_degree();
  for (int i = 0; i < max_degree; i++)
    v[i].set_degree(i + 1);
  for (int i = max_degree; i < sz; i++)
    v[i].set_degree(max_degree + 1);
}

void MaxClique::Vertices::set_degrees(MaxClique &m) {
  for (int i = 0; i < sz; i++) {
    int d = 0;
    for (int j = 0; j < sz; j++)
      if (m.connection(v[i].get_i(), v[j].get_i()))
        d++;
    v[i].set_degree(d);
  }
}

void MaxClique::init(const MaxClique::SearchType t) {
  V.set_degrees(*this);
  V.sort(__how_to_sort, __energy);
  V.init_colors();
  if (t == MCQD || t == MCQDW) {
    for (int i = 0; i < V.size() + 1; i++) {
      S[i].set_i1(0);
      S[i].set_i2(0);
    }
  }
}

bool MaxClique::cut1(const int pi, const ColorClass &A) {
  for (int i = 0; i < A.size(); i++)
    if (connection(pi, A.at(i)))
      return true;
  return false;
}

void MaxClique::cut2(const Vertices &A, Vertices &B) {
  for (int i = 0; i < A.size() - 1; i++) {
    if (connection(A.end().get_i(), A.at(i).get_i()))
      B.push(A.at(i).get_i());
  }
}

void MaxClique::color_sort_mcq(Vertices &R) {
  int j = 0;
  int maxno = 1;
  int min_k = QMAX.size() - Q.size() + 1;
  C[1].rewind();
  C[2].rewind();
  int k = 1;
  for (int i = 0; i < R.size(); i++) {
    int pi = R.at(i).get_i();
    k = 1;
    while (cut1(pi, C[k]))
      k++;
    if (k > maxno) {
      maxno = k;
      C[maxno + 1].rewind();
    }
    C[k].push(pi);
    if (k < min_k) {
      R.at(j++).set_i(pi);
    }
  }
  if (j > 0)
    R.at(j - 1).set_degree(0);
  if (min_k <= 0)
    min_k = 1;

  for (k = min_k; k <= maxno; k++) {
    for (int i = 0; i < C[k].size(); i++) {
      R.at(j).set_i(C[k].at(i));
      R.at(j++).set_degree(k);
    }
  }
}

void MaxClique::color_sort_mcqw(Vertices &R) {
  int j = 0;
  int maxno = 1;
  int min_k = QMAX.size() - Q.size() + 1;
  C[1].rewind();
  C[2].rewind();
  int k = 1;
  for (int i = 0; i < R.size(); i++) {
    int pi = R.at(i).get_i();
    k = 1;
    while (cut1(pi, C[k]))
      k++;
    if (k > maxno) {
      maxno = k;
      C[maxno + 1].rewind();
    }
    C[k].set_min_energy(__energy[pi]);
    C[k].push(pi);
    if (k < min_k) {
      R.at(j++).set_i(pi);
    }
  }
  if (j > 0)
    R.at(j - 1).set_degree(0);
  if (min_k <= 0)
    min_k = 1;
  vector<double> min_energy_k;
  for (k = 1; k < min_k; ++k) {
    min_energy_k.push_back(C[k].get_min_energy());
  }
  for (k = min_k; k <= maxno; k++) {
    min_energy_k.push_back(C[k].get_min_energy());
    const double ene =
        std::accumulate(min_energy_k.begin(), min_energy_k.end(), 0.0);
    for (int i = 0; i < C[k].size(); i++) {
      R.at(j).set_energy(ene);
      R.at(j).set_i(C[k].at(i));
      R.at(j++).set_degree(k);
    }
  }
}

void MaxClique::color_sort_kcqw(Vertices &R) {
  int j = 0;
  int maxno = 1;
  int min_k = QMAX.size() - Q.size() + 1;
  C[1].rewind();
  C[2].rewind();
  int k = 1;
  for (int i = 0; i < R.size(); i++) {
    int pi = R.at(i).get_i();
    k = 1;
    while (cut1(pi, C[k]))
      k++;
    if (k > maxno) {
      maxno = k;
      C[maxno + 1].rewind();
    }
    C[k].set_min_energy(__energy[pi]);
    C[k].push(pi);
    if (k < min_k) {
      R.at(j++).set_i(pi);
    }
  }
  if (j > 0)
    R.at(j - 1).set_degree(0);
  if (min_k <= 0)
    min_k = 1;
  vector<double> min_energy_k;
  for (k = 1; k < min_k; ++k) {
    min_energy_k.push_back(C[k].get_min_energy());
  }
  for (k = min_k; k <= maxno; k++) {
    min_energy_k.push_back(C[k].get_min_energy());
    const double ene = sum_kcqw(min_energy_k, Q.size() + k - (QMAX.size() + 1));
    for (int i = 0; i < C[k].size(); i++) {
      R.at(j).set_energy(ene);
      R.at(j).set_i(C[k].at(i));
      R.at(j++).set_degree(k);
    }
  }
}

double MaxClique::sum_kcqw(const vector<double> &min_energy_k, int discard) {
  vector<double> energy(min_energy_k.begin(), min_energy_k.end());
  if (discard > 0) {
    std::sort(energy.begin(), energy.end());
    // remove vertices with worse energies
    while (discard-- > 0)
      energy.pop_back();
  }
  return std::accumulate(energy.begin(), energy.end(), 0.0);
}

bool MaxClique::ub_is_breached(const vector<double> &E, const double ene) {
  for (int i = 0; i < E.size(); ++i) {
    if (ene < E[i] - 0.001) {
      clog << "UB BREACHED ENERGY(Q)=" << ene << " UB(" << i << ")=" << E[i]
           << endl;
      return true;
    }
  }
  return false;
}

void MaxClique::check_ub(const vector<double> &E) {
  for (int i = 0; i < E.size(); ++i) {
    clog << "LEVEL=" << i << " UB=" << E[i] << endl;
    if (i < E.size() - 1 && E[i] > E[i + 1] + 0.0001) {
      clog << "INCREASING UPPER BOUND: E[" << i << "]=" << E[i] << " E["
           << i + 1 << "]=" << E[i + 1] << endl;
      //			throw Error("Upper bounds should be decreasing
      //with levels");
    }
  }
}

double MaxClique::tighten_ub(const vector<double> &E, const double ub) {
  for (int i = 0; i < E.size(); ++i) {
    if (E[i] > ub) {
      //			clog << "FOUND A TIGHTER UPPER BOUND: E[" << i
      //<< "]=" << E[i] << " than the current UB=" << ub << endl;
      return E[i];
    }
  }
  return ub;
}

void MaxClique::expand_mcq(Vertices R) {
  if (__max_steps != -1 && pk > __max_steps)
    return;
  while (R.size()) {
    if (Q.size() + R.end().get_degree() > QMAX.size()) {
      Q.push(R.end().get_i());
      Vertices Rp(R.size());
      cut2(R, Rp);
      if (Rp.size()) {
        color_sort_mcq(Rp);
        pk++;
        expand_mcq(Rp);
      } else if (Q.size() > QMAX.size()) {
        QMAX = Q;
      }
      Rp.dispose();
      Q.pop();
    } else {
      return;
    }
    R.pop();
  }
}

void MaxClique::expand_mcqd(Vertices R) {
  if ((__max_steps != -1 && pk > __max_steps) || (((clock()-t) / CLOCKS_PER_SEC) > __max_time && __max_time != -1))
    return;
  S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() -
                  S[level].get_i2());
  S[level].set_i2(S[level - 1].get_i1());
  while (R.size() && (((clock()-t) / CLOCKS_PER_SEC) <  __max_time && __max_time != -1)) {
    if (Q.size() + R.end().get_degree() > QMAX.size()) {
      Q.push(R.end().get_i());
      Vertices Rp(R.size());
      cut2(R, Rp);
      if (Rp.size()) {
        if ((float)S[level].get_i1() / ++pk < Tlimit) {
          degree_sort(Rp);
        }
        color_sort_mcq(Rp);
        S[level].inc_i1();
        level++;
        expand_mcqd(Rp);
        level--;
      } else if (Q.size() > QMAX.size()) {
        QMAX = Q;
      }
      Rp.dispose();
      Q.pop();
    } else {
      return;
    }
    R.pop();
  }
}

void MaxClique::expand_mcqw(Vertices R) {
  if (__max_steps != -1 && pk > __max_steps)
    return;
  while (R.size()) {
    if (Q.get_energy() + R.end().get_energy() < QMAX.get_energy()) {
      Q.push(R.end().get_i());
      Q.increase_energy(__energy[R.end().get_i()]);
      Vertices Rp(R.size());
      cut2(R, Rp);
      if (Rp.size()) {
        color_sort_mcqw(Rp);
        pk++;
        expand_mcqw(Rp);
      } else if (Q.get_energy() > QMAX.get_energy()) {
        QMAX = Q;
      }
      Rp.dispose();
      Q.pop();
      Q.decrease_energy(__energy[R.end().get_i()]);
    } else {
      return;
    }
    R.pop();
  }
}

void MaxClique::expand_mcqdw(Vertices R) {
  if (__max_steps != -1 && pk > __max_steps)
    return;
  S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() -
                  S[level].get_i2());
  S[level].set_i2(S[level - 1].get_i1());
  while (R.size()) {
    if (Q.get_energy() + R.end().get_energy() < QMAX.get_energy()) {
      Q.push(R.end().get_i());
      Q.increase_energy(__energy[R.end().get_i()]);
      Vertices Rp(R.size());
      cut2(R, Rp);
      if (Rp.size()) {
        if ((float)S[level].get_i1() / ++pk < Tlimit) {
          degree_sort(Rp);
        }
        color_sort_mcqw(Rp);
        S[level].inc_i1();
        level++;
        expand_mcqdw(Rp);
        level--;
      } else if (Q.get_energy() > QMAX.get_energy()) {
        QMAX = Q;
      }
      Rp.dispose();
      Q.pop();
      Q.decrease_energy(__energy[R.end().get_i()]);
    } else {
      return;
    }
    R.pop();
  }
}

void MaxClique::expand_kcq(Vertices R) {
  if (__max_steps != -1 && pk > __max_steps)
    return;
  while (R.size()) {
    if (Q.size() + R.end().get_degree() > QMAX.size()) {
      Q.push(R.end().get_i());
      Vertices Rp(R.size());
      cut2(R, Rp);
      if (Rp.size()) {
        color_sort_mcq(Rp);
        pk++;
        expand_kcq(Rp);
      } else if (Q.size() == QMAX.size() + 1) {
        __best_cliques.insert(Q);
      }
      Rp.dispose();
      Q.pop();
    } else {
      return;
    }
    R.pop();
  }
}

void MaxClique::expand_kcqd(Vertices R) {
  if (__max_steps != -1 && pk > __max_steps)
    return;
  S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() -
                  S[level].get_i2());
  S[level].set_i2(S[level - 1].get_i1());
  while (R.size()) {
    if (Q.size() + R.end().get_degree() > QMAX.size()) {
      Q.push(R.end().get_i());
      Vertices Rp(R.size());
      cut2(R, Rp);
      if (Rp.size()) {
        if ((float)S[level].get_i1() / ++pk < Tlimit) {
          degree_sort(Rp);
        }
        color_sort_mcq(Rp);
        S[level].inc_i1();
        level++;
        expand_kcqd(Rp);
        level--;
      } else if (Q.size() == QMAX.size() + 1) {
        __best_cliques.insert(Q);
      }
      Rp.dispose();
      Q.pop();
    } else {
      return;
    }
    R.pop();
  }
}

void MaxClique::expand_kcqw(Vertices R) {
  if (__max_steps != -1 && pk > __max_steps)
    return;
  while (R.size()) {
    if (Q.size() + R.end().get_degree() > QMAX.size() &&
        Q.get_energy() + R.end().get_energy() <
            __best_cliques.get_highest_energy()) {
      Q.push(R.end().get_i());
      Q.increase_energy(__energy[R.end().get_i()]);
      Vertices Rp(R.size());
      cut2(R, Rp);
      if (Rp.size()) {
        color_sort_kcqw(Rp);
        pk++;
        expand_kcqw(Rp);
      } else if (Q.size() == QMAX.size() + 1) {
        __best_cliques.insert_if_better(Q);
      }
      Rp.dispose();
      Q.pop();
      Q.decrease_energy(__energy[R.end().get_i()]);
    } else {
      return;
    }
    R.pop();
  }
}

void MaxClique::expand_kcqdw(Vertices R) {
  if (__max_steps != -1 && pk > __max_steps)
    return;
  S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() -
                  S[level].get_i2());
  S[level].set_i2(S[level - 1].get_i1());
  while (R.size()) {
    if (Q.size() + R.end().get_degree() > QMAX.size() &&
        Q.get_energy() + R.end().get_energy() <
            __best_cliques.get_highest_energy()) {
      Q.push(R.end().get_i());
      Q.increase_energy(__energy[R.end().get_i()]);
      Vertices Rp(R.size());
      cut2(R, Rp);
      if (Rp.size()) {
        if ((float)S[level].get_i1() / ++pk < Tlimit) {
          degree_sort(Rp);
        }
        color_sort_kcqw(Rp);
        S[level].inc_i1();
        level++;
        expand_kcqdw(Rp);
        level--;
      } else if (Q.size() == QMAX.size() + 1) {
        __best_cliques.insert_if_better(Q);
      }
      Rp.dispose();
      Q.pop();
      Q.decrease_energy(__energy[R.end().get_i()]);
    } else {
      return;
    }
    R.pop();
  }
}

void MaxClique::BestCliques::insert_if_better(const MaxClique::Clique &clique) {
  if (QMAXES.size() < maxsz) {
    QMAXES.push(clique);
  } else {
    const double highest_energy = QMAXES.top().get_energy();
    if (clique.get_energy() < highest_energy) {
      QMAXES.pop();        // remove worst energy clique
      QMAXES.push(clique); // add a better clique
      dbgmsg("inserting clique of size " + help::to_string(clique.size()) +
             " with energy " + help::to_string(clique.get_energy()));
    }
  }
}

vector<MaxClique::Clique> MaxClique::BestCliques::get_k_cliques() const {
  Queue qmaxes(QMAXES);
  vector<Clique> result;
  while (!qmaxes.empty()) {
    result.push_back(qmaxes.top());
    qmaxes.pop();
  }
  return result;
}

/*
 * Find one maximum (weight) clique.
 * MCQ ... use maximum clique algorithm
 * MCQD ... use maximum clique with dynamic upper bounds algorithm
 * (MaxCliqueDyn) MCQW ... use maximum weight clique algorithm MCQWD ... use
 * maximum weight clique aklgorithm with dynamic upper bounds
 * (MaxCliqueDynWeights)
 */
MaxClique::Clique MaxClique::maximum_clique(MaxClique::SearchType t) {
  init(t);
  switch (t) {
  case MCQ:
    expand_mcq(V);
    break;
  case MCQD:
    expand_mcqd(V);
    break;
  case MCQW:
    color_sort_mcqw(V);
    expand_mcqw(V);
    break;
  case MCQDW:
    color_sort_mcqw(V);
    expand_mcqdw(V);
    break;
  };
  return QMAX;
}

/*
 * Find all cliques with size >= k
 * MCQ ... use maximum clique algorithm
 * MCQD ... use maximum clique with dynamic upper bounds algorithm
 * (MaxCliqueDyn) MCQW ... fuse maximum weight clique algorithm MCQWD ... use
 * maximum weight clique aklgorithm with dynamic upper bounds
 * (MaxCliqueDynWeights)
 */
vector<MaxClique::Clique> MaxClique::k_cliques(const int k,
                                               const MaxClique::SearchType t) {
  // search for all cliques of size >= k
  for (int i = 0; i < k - 1; ++i)
    QMAX.push(0);
  init(t);
  switch (t) {
  case MCQ:
    expand_kcq(V);
    break;
  case MCQD:
    expand_kcqd(V);
    break;
  case MCQW:
    color_sort_kcqw(V);
    expand_kcqw(V);
    break;
  case MCQDW:
    color_sort_kcqw(V);
    expand_kcqdw(V);
    break;
  };
  return __best_cliques.get_k_cliques();
}
/*
void MaxClique::output_graph(const Array2d<bool> &conn,
                             const vector<double> &energy) {
  stringstream ss;
  int nodes = conn.get_szi();
  int edges = 0;
  for (int i = 0; i < conn.get_szi(); i++) {
    ss << "n " << i << " " << energy[i] << endl;
  }
  for (int i = 0; i < conn.get_szi(); i++) {
    for (int j = i + 1; j < conn.get_szj(); j++) {
      if (conn.get(i, j)) {
        ss << "e " << i << " " << j << endl;
        ++edges;
      }
    }
  }
  stringstream result;
  result << "p edge " << nodes << " " << edges << endl << ss.str();

  inout::output_file(result.str(), "graph.txt");
}
*/
/*
void MaxClique::read_graph(Array2d<bool> &conn, vector<double> &energy) {

  vector<string> graph;
  inout::read_file("graph.txt", graph);

  for (auto &line : graph) {
    stringstream ss(line);
    char ch;
    if (line[0] == 'p') {
      string keyword;
      int nodes, edges;
      ss >> ch >> keyword >> nodes >> edges;
      conn.init(nodes, nodes);
    } else if (line[0] == 'n') {
      int node;
      double ene;
      ss >> ch >> node >> ene;
      energy.push_back(ene);
    } else if (line[0] == 'e') {
      int node1, node2;
      ss >> ch >> node1 >> node2;
      conn.set(node1, node2);
      conn.set(node2, node1);
    }
  }
}

void MaxClique::read_graph(Array2d<bool> &conn, std::string path) {

  vector<string> graph;
  inout::read_file(path, graph);

  for (auto &line : graph) {
    stringstream ss(line);
    char ch;
    if (line[0] == 'p') {
      string keyword;
      int nodes, edges;
      ss >> ch >> keyword >> nodes >> edges;
      conn.init(nodes, nodes);
    } else if (line[0] == 'e') {
      int node1, node2;
      ss >> ch >> node1 >> node2;
      conn.set(node1, node2);
      conn.set(node2, node1);
    }
  }
}
*/
void MaxClique::Clique::print(const vector<double> &energy) {
  clog << "CLIQUE SIZE=" << this->vertices.size() << " TOTAL_E=" << this->energy
       << " " << endl;
  for (int i = 0; i < this->vertices.size(); ++i) {
    if (!energy.empty()) {
      clog << this->vertices[i] << "(" << energy.at(this->vertices[i]) << ") ";
    } else {
      clog << this->vertices[i] << " ";
    }
  }
  clog << endl;
}

void MaxClique::Vertices::print(const vector<double> &energy) {
  clog << "PRINTING VERTICES" << endl;
  for (int i = 0; i < this->size(); ++i) {
    const int pi = this->at(i).get_i();
    if (!energy.empty()) {
      clog << "VERTEX=" << pi << " ENERGY=" << energy.at(pi)
           << " UB=" << this->at(i).get_energy()
           << " DEGREE=" << this->at(i).get_degree() << endl;
    } else {
      clog << "VERTEX=" << pi << " DEGREE=" << this->at(i).get_degree() << endl;
    }
  }
}