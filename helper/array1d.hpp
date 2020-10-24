#ifndef ARRAY1D_H
#define ARRAY1D_H

#include <assert.h>
#include <stdio.h>
#include <string.h>

template <typename T> struct Array1d {
  T *data;
  int sz;

  Array1d operator-() const {
    Array1d response(*this);
    for (int i = 0; i < this->sz; i++) {
      response.data[i] = -response.data[i];
    }
    return response;
  }

  Array1d operator/(const Array1d &other) {
    dbgmsg("Array1d divison operator");
    assert(this->sz == other.sz);
    Array1d response(*this);
    for (int i = 0; i < this->sz; i++) {
      if (other.data[i] != 0) {
        response.data[i] /= other.data[i];
      }
    }
    return response;
  }

  Array1d operator+(const Array1d &other) {
    dbgmsg("Array1d plus operator");
    assert(this->sz == other.sz);
    Array1d response(other);
    for (int i = 0; i < this->sz; i++) {
      response.data[i] += this->data[i];
    }
    return response;
  }

  Array1d operator*(const double &right) const {
    Array1d response(*this);
    for (int i = 0; i < this->sz; i++) {
      response.data[i] *= right;
    }
    return response;
  }

  void init(int SZ) {
    sz = SZ;
    data = new T[sz];
    memset(data, 0, sz * sizeof(T));
  }

  Array1d(int SZ) : data(nullptr) { init(SZ); }

  Array1d() : data(nullptr), sz(0) {}

  Array1d(const Array1d &other) { // copy
    dbgmsg("Array1d copy constructor");
    sz = other.sz;
    data = new T[sz];
    memcpy(data, other.data, sz * sizeof(T));
  }

  ~Array1d() {
    if (data) {
      delete[] data;
      data = nullptr;
    }
  }

  void reset() { memset(data, 0, sz * sizeof(T)); }
};

template <typename T>
ostream &operator<<(ostream &stream, const Array1d<T> &s) {
  stream << "(";
  for (int i = 0; i < s.sz; ++i) {
    stream << s.data[i] << (i + 1 == s.sz ? "" : ",");
  }
  stream << ")" << endl;
  return stream;
}

#endif
