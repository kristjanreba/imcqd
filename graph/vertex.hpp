#ifndef VERTEX_H
#define VERTEX_H
#include "helper/debug.hpp"
#include "helper/error.hpp"
#include "molecule/it.hpp"
#include <algorithm>
#include <boost/regex.hpp>
#include <functional>
#include <memory>
#include <queue>
#include <set>

namespace glib {

template <typename T>
class Vertex : public template_vector_container<Vertex<T>> {
  T &__vertex;

public:
  template <typename T> Vertex(const T &vertex) : __vertex(vertex) {}

  T &get_vertex() const { return __vertex; }
  friend ostream &operator<<(ostream &stream, const Vertex &v);

  bool compatible(const Vertex &other) const {
    return __vertex.compatible(other.get_vertex());
  }
  const string &get_label() const { return __vertex.get_label(); }
  const string &print() const { return __vertex.get_label(); }
  const int weight() const { return __vertex.weight; }
};
}; // namespace glib

#endif