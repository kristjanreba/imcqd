#include "path.hpp"
#include <iostream>
#include <stdio.h>
#include <string>

namespace path {

bool file_exists(const std::string &file_path) {
  std::filesystem::path p(file_path);
  return std::filesystem::exists(p);
}

std::string join(const std::string &str1, const std::string &str2) {
  std::filesystem::path p1(str1);
  std::filesystem::path p2(str2);
  p1 /= p2;
  return p1.string();
}

std::string get_stem(const std::string &str) {
  std::filesystem::path p(str);
  return p.stem().string();
}

bool has_suffix(const std::string &filename, const std::string &suffix) {

  std::string filename_upper(filename);

  std::transform(filename_upper.begin(), filename_upper.end(),
                 filename_upper.begin(),
                 [](unsigned char c) { return std::toupper(c); });

  std::string suffix_upper(suffix);

  std::transform(suffix_upper.begin(), suffix_upper.end(), suffix_upper.begin(),
                 [](unsigned char c) { return std::toupper(c); });

  const bool ret =
      filename_upper.size() >= suffix_upper.size() &&
      filename_upper.compare(filename_upper.size() - suffix_upper.size(),
                             suffix_upper.size(), suffix_upper) == 0;

  dbgmsg("file " << filename << " has suffix " << suffix << " = " << boolalpha
                 << ret);

  return ret;
}

}; // namespace path
