#pragma once

#include <exception>
#include <sstream>
#include <string>
using namespace std;

class Error : public exception {
  const string __msg;

public:
  Error(const string &msg) : __msg(msg) {}
  ~Error() throw() {}
  const char *what() const noexcept { return __msg.c_str(); }
};
