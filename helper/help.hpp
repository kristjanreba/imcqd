#pragma once

#include "debug.hpp"
#include "error.hpp"
#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

namespace help {
string memusage(const string &);
std::tuple<double, double, double> gnuplot(const double &x1, const double &x2,
                                           const string &datapoints);

struct edge {
  string atom_property1;
  string atom_property2;
  string bond_property;
};

typedef vector<edge> smiles;

ostream &operator<<(ostream &os, const smiles &edges);

struct rename_rule {
  smiles pattern;
  vector<string> rule;
};

typedef vector<rename_rule> rename_rules;

ostream &operator<<(ostream &os, const rename_rule &rule);

enum IdatmGeometry {
  Ion = 0,
  Single = 1,
  Linear = 2,
  Planar = 3,
  Tetrahedral = 4,
  TrigonalBipyramidal = 5,
  TetragonalBipyramidal = 6,
  High = 9
};
struct IdatmEntry {
  IdatmGeometry geometry;
  int substituents;
  string description;
};
typedef map<const string, const IdatmEntry> IdatmInfoMap;

template <class Set1, class Set2>
bool is_disjoint(const Set1 &set1, const Set2 &set2) {
  if (set1.empty() || set2.empty())
    return true;
  typename Set1::const_iterator it1 = set1.begin(), it1End = set1.end();
  typename Set2::const_iterator it2 = set2.begin(), it2End = set2.end();
  if (*it1 > *set2.rbegin() || *it2 > *set1.rbegin())
    return true;
  while (it1 != it1End && it2 != it2End) {
    if (*it1 == *it2)
      return false;
    if (*it1 < *it2) {
      it1++;
    } else {
      it2++;
    }
  }
  return true;
}

template <typename T> string to_string(T num) {
  ostringstream ss;
  ss << num;
  return ss.str();
}

string dtos(double d, int precision = 3); // double to string of precision

string replace_str_char(string str, const string &replace, char ch);

template <typename T = string>
auto ssplit(
    const string &source, const char *delimiter = " ", bool keepEmpty = false,
    T convert_to(const string &) = [](const string &s) { return s; }) {
  vector<T> results;
  size_t prev = 0;
  size_t next = 0;
  while ((next = source.find_first_of(delimiter, prev)) != string::npos) {
    if (keepEmpty || (next - prev != 0)) {
      results.push_back(convert_to(source.substr(prev, next - prev)));
    }
    prev = next + 1;
  }
  if (prev < source.size()) {
    results.push_back(convert_to(source.substr(prev)));
  }
  return results;
}

char get_one_letter(string s);
string get_three_letter(const char c);

const IdatmEntry &get_info_map(const string &name);
const vector<string> &get_gaff_replacement(const string &name);
vector<vector<string>> get_replacement(const vector<string> &initial);

const string EW = "^N|^O|F|Cl|Br"; // electron withdrawing atoms
const string XX = "^C|^N|^O|^S|^P";
const string XA = "^O|^S";
const string XB = "^N|^P";
const string XC = "F|Cl|Br|I";
const string XD = "^S|^P";

extern const map<const string, const string> gaff_flip;
extern const set<string> gaff_group_1;
extern const set<string> gaff_group_2;

extern const vector<string> idatm_unmask;
extern const map<pair<string, string>, double> repulsion_idx;
extern const map<const string, const int> idatm_mask;
extern const map<const string, const double> element_radius;
extern const double vdw_radius[];

extern const set<string> amino_acids;
extern const set<string> non_standard_amino_acids;
extern const set<string> nucleic_acids;
extern const set<string> ions;
extern const map<const string, const char> heavy_atoms;
extern const map<const string, const char> protein_hydrogen_atoms;
extern const map<const string, const char> dna_hydrogen_atoms;
extern const map<const string, const string> buffer_ions;
extern const map<const string, const pair<string, string>> non_specific_binders;
extern const map<const string, const pair<string, string>>
    non_specific_ion_binders;
extern const map<const string, const map<string, string>> standard_residues;
extern const IdatmInfoMap infoMap;

extern const map<const char, const string> three_letter;
extern const map<const string, string> sybyl;
extern const map<string, const char> one_letter;

extern const rename_rules special;
extern const rename_rules refine;
extern const rename_rules bond_gaff_type;

extern const rename_rules rotatable;
extern const map<string, vector<string>> gaff_replacement;
extern const rename_rules gaff;
extern const rename_rules atomic_penalty_scores;

}; // namespace help
