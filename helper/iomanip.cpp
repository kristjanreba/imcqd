#include "iomanip.hpp"
#include <iostream>

namespace molib {

ostream &IOmanip::reset(ostream &os) {
  IOmanip::set_IOtype(os, IOmanip::default_t);
  return os;
}

ostream &IOmanip::pdb(ostream &os) {
  IOmanip::set_IOtype(os, IOmanip::pdb_t);
  return os;
}

ostream &IOmanip::mol2(ostream &os) {
  IOmanip::set_IOtype(os, IOmanip::mol2_t);
  return os;
}

ostream &IOmanip::extended(ostream &os) {
  IOmanip::or_IOtype(os, IOmanip::extended_t);
  return os;
}

ostream &IOmanip::cif(ostream &os) {
  IOmanip::or_IOtype(os, IOmanip::cif_t);
  return os;
}

}; // namespace molib
