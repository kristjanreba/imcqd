#pragma once

#include "helper/help.hpp"
#include <filesystem>
#include <string>

namespace path {

bool file_exists(const std::string &file_path);

std::string join(const std::string &str1, const std::string &str2);

template <typename... Targs>
std::string join(const std::string &str1, const std::string &str2,
                 Targs... Fargs) {
  std::filesystem::path p1(str1);
  std::filesystem::path p2(str2);
  p1 /= p2;
  return path::join(p1.string(), Fargs...);
}

std::string get_stem(const std::string &str);

template <typename T>
std::string insert_to_filename(const std::string &filename, T value) {
  std::string fn(filename);
  return fn.insert(fn.find_first_of("."),
                   std::string("_" + help::to_string(value)));
}

template <typename T, typename... Targs>
std::string insert_to_filename(const std::string &filename, T value,
                               Targs... args) {
  std::string fn(insert_to_filename(filename, value));
  return insert_to_filename(fn, args...);
}

bool has_suffix(const std::string &filename, const std::string &suffix);

}; // namespace path
