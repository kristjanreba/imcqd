#ifndef IOMANIP_H
#define IOMANIP_H

#include <iostream>
using namespace std;

namespace molib {

struct IOmanip {
public:
  enum IOtype {
    default_t = 0,
    pdb_t = 1,
    mol2_t = 2,
    extended_t = 4,
    cif_t = 8
  };

  static void or_IOtype(ios_base &s, IOtype o) { flag(s) |= o; }
  static void set_IOtype(ios_base &s, IOtype o) { flag(s) = o; }
  static IOtype get_IOtype(ios_base &s) { return (IOtype)flag(s); }

  static ostream &reset(ostream &os);
  static ostream &pdb(ostream &os);
  static ostream &mol2(ostream &os);
  static ostream &extended(ostream &os);
  static ostream &cif(ostream &os);

private:
  static long &flag(ios_base &s) {
    static int n = ios_base::xalloc();
    return s.iword(n);
  }
};

}; // namespace molib

#endif
