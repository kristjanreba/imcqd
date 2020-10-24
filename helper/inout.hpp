#pragma once

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#include <sys/locking.h>
#else
#include <sys/file.h> /* for flock(2) */
#include <unistd.h>   /* for close(2) prototypes */
#endif

#include <algorithm>
#include <errno.h> /* for errno prototype */
#include <fstream>
#include <functional>
#include <iterator>
#include <sstream>
#include <stdio.h>  /* for fprintf(3),printf(3),stderr protype */
#include <string.h> /* for strerror(3) prototype */
#include <string>
#include <sys/stat.h> /* for S_* constants */
#include <vector>

namespace inout {

enum f_not_found { panic = 0, no_panic = 1 };

void read_file(const std::string &name, std::vector<std::string> &, f_not_found = panic,
               const int num_occur = -1,
               const std::string &pattern = ""); // throws Error

void read_file(const std::string &name, std::string &, f_not_found = panic,
               const int num_occur = -1,
               const std::string &pattern = ""); // throws Error

void read_file(const std::string &name, std::string &, std::streampos &pos_in_file,
               f_not_found = panic, const int num_occur = -1,
               const std::string &pattern = ""); // throws Error

void read_file(const std::string &name, std::vector<std::string> &, std::streampos &pos_in_file,
               f_not_found = panic, const int num_occur = -1,
               const std::string &pattern = ""); // throws Error

void file_open_put_contents(const std::string &name, const std::vector<std::string> &v,
                            std::ios_base::openmode = std::ios_base::out);

void file_open_put_stream(const std::string &name, const std::stringstream &ss,
                          std::ios_base::openmode = std::ios_base::out);

std::vector<std::string> files_matching_pattern(const std::string &,
                                      const std::string &pattern = "");

template <typename S> struct out_manipulator : public std::function<S &(S &)> {
  template <typename T>
  out_manipulator(T &&t)
      : std::function<S &(S &)>([=](S &i) -> S & { return i << t; }) {}
  template <typename T>
  out_manipulator(T *t)
      : std::function<S &(S &)>([=](S &i) -> S & { return i << t; }) {
  } // for g++
  template <typename U> friend U &operator<<(U &u, out_manipulator &a) {
    return static_cast<U &>(a(u));
  }
};

template <class T>
void output_file(const T &anything, const std::string &filename,
                 const std::vector<inout::out_manipulator<std::ostream>> &manips =
                     std::vector<inout::out_manipulator<std::ostream>>(),
                 std::ios_base::openmode mode = std::ios_base::out) {
  std::stringstream ss;
  for (auto m : manips)
    ss << m;
  ss << anything;
  inout::file_open_put_stream(filename, ss, mode);
}

}; // namespace inout
