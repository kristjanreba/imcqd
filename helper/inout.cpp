#include "inout.hpp"
#include "debug.hpp"
#include "error.hpp"
#include "helper/gzip.hpp"
#include "helper/help.hpp"
#include "helper/path.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/stream.hpp> // for stream
#include <boost/regex.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace inout {

void __mkdir(const std::string &dir_path) {
  // makes a path "janez/aska/mia from e.g. "janez/aska/mia/test.txt"
  std::filesystem::path dir(dir_path);
  dir.remove_filename();
  if (!dir.string().empty()) {
    if (!std::filesystem::exists(dir)) {
      if (!std::filesystem::create_directories(dir)) {
        throw Error("WHOOPS : cannot create directory " + dir.string());
      }
    }
  }
}

int __lock(const std::string &name) {
#ifdef _WIN32
  //			int fd;
  //			_sopen_s(&fd, name.c_str(), _O_RDWR | _O_CREAT, 0,
  //_S_IREAD | _S_IWRITE);
  return 0;
#else
  int fd = open(name.c_str(), O_RDWR | O_CREAT,
                S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
  if (flock(fd, LOCK_EX) == -1) {
  }
  return fd;
#endif
}

void __unlock(int fd) {
#ifdef _WIN32
//			if (_close(fd) == -1) {
//			}
#else
  if (flock(fd, LOCK_UN) == -1) {
  }
  if (close(fd) == -1) {
  }
#endif
}

void __read_txt_file(const std::string &name, std::vector<std::string> &s,
                     std::streampos &pos_in_file, const int num_occur,
                     const std::string &pattern) {

  dbgmsg(pos_in_file);

  int fd = __lock(name);
  std::ifstream in(name);

  in.seekg(pos_in_file);

  std::string line;
  int i = 0;
  std::streampos pos = in.tellg();

  while (std::getline(in, line)) {
    if (num_occur != -1 &&
        (line.find(pattern) != std::string::npos && i++ == num_occur)) {
      dbgmsg("breaking on line = " << line);
      break;
    } else
      pos = in.tellg();
    s.push_back(line);
  }

  pos_in_file = pos;

  in.close();
  __unlock(fd);
}

void __read_gz_file(const std::string &name, std::vector<std::string> &s,
                    std::streampos &pos_in_file, const int num_occur,
                    const std::string &pattern) {

  dbgmsg("first position (line number) = "
         << pos_in_file); // gzip stream doesn't support seeking

  int fd = __lock(name);
  std::ifstream file(name, std::ios_base::in | std::ios_base::binary);

  boost::iostreams::filtering_istream in;
  in.push(boost::iostreams::gzip_decompressor());
  in.push(file);

  std::string line;
  std::streampos ip = 0;
  while (ip != pos_in_file && std::getline(in, line)) {
    ip += 1;
  }

  int i = 0;
  std::streampos pos = ip;

  while (std::getline(in, line)) {
    dbgmsg("unzipped line = " << line);
    dbgmsg("current position (adding 2) = " << pos);
    if (num_occur != -1 &&
        (line.find(pattern) != std::string::npos && i++ == num_occur)) {
      dbgmsg("breaking on line = " << line);
      break;
    } else {
      pos += 1;
    }
    s.push_back(line);
  }

  pos_in_file = pos;
  dbgmsg("position after reading = " << pos_in_file);

  boost::iostreams::close(in);
  file.close();
  __unlock(fd);
}

void read_file(const std::string &name, std::string &s, f_not_found w,
               const int num_occur, const std::string &pattern) {
  std::streampos pos_in_file = 0;
  read_file(name, s, pos_in_file, w, num_occur, pattern);
}
void read_file(const std::string &name, std::vector<std::string> &s,
               f_not_found w, const int num_occur, const std::string &pattern) {
  std::streampos pos_in_file = 0;
  read_file(name, s, pos_in_file, w, num_occur, pattern);
}
void read_file(const std::string &name, std::string &s,
               std::streampos &pos_in_file, f_not_found w, const int num_occur,
               const std::string &pattern) {
  std::vector<std::string> vs;
  read_file(name, vs, pos_in_file, w, num_occur, pattern);
  s = boost::algorithm::join(vs, "\n");
}

void read_file(const std::string &name, std::vector<std::string> &s,
               std::streampos &pos_in_file, f_not_found w, const int num_occur,
               const std::string &pattern) {

  std::string filename(name);
  transform(filename.begin(), filename.end(), filename.begin(), ::toupper);

  bool gzipped = path::has_suffix(filename, ".GZ");

  // check if file can be opened, since lock creates the file if it does not
  // exist...
  std::ifstream f(name);
  if (!f.is_open() && w == panic) {
    f.close();
    throw Error("WHOOPS : cannot open file " + name + "\n");
  }
  f.close();

  try {

    if (gzipped) {
      __read_gz_file(name, s, pos_in_file, num_occur, pattern);
    } else {
      __read_txt_file(name, s, pos_in_file, num_occur, pattern);
    }

  } catch (std::exception &e) {
    if (w == panic)
      throw Error("WHOOPS : an error occured during reading file " + name +
                  " with errmsg " + e.what());
  }
}

void file_open_put_contents(const std::string &name,
                            const std::vector<std::string> &v,
                            std::ios_base::openmode mode) {
  std::stringstream ss;
  for (auto &s : v)
    ss << s << endl;
  file_open_put_stream(name, ss, mode);
}

void file_open_put_stream(const std::string &name, const std::stringstream &ss,
                          std::ios_base::openmode mode) {
  __mkdir(name);
  int fd = __lock(name);
  ofstream output_file(name, mode);

  if (!output_file.is_open()) {
    __unlock(fd);
    throw Error("WHOOPS : cannot open output file: " +
                name); // how to emulate $! ?
  }

  std::string filename(name);
  transform(filename.begin(), filename.end(), filename.begin(), ::toupper);
  bool gzipped = path::has_suffix(filename, ".GZ");
  if (gzipped) {
    output_file << Gzip::compress(ss.str());
  } else {
    output_file << ss.str();
  }
  output_file.close();
  __unlock(fd);
}

std::vector<std::string> files_matching_pattern(const std::string &file_or_dir,
                                                const std::string &pattern) {

  std::vector<std::string> fn;
  std::filesystem::path p(file_or_dir);

  if (std::filesystem::is_directory(p)) {

    std::vector<std::filesystem::path> files;
    copy(std::filesystem::directory_iterator(p),
         std::filesystem::directory_iterator(),
         back_inserter(files)); // get all files in directoy...

    for (auto &p : files) {
      fn.push_back(p.string());
    }

  } else {
    throw Error("WHOOPS : cannot open directory " + p.string());
  }

  std::vector<std::string> filenames;
  for (std::string &s : fn) {
    if (boost::regex_search(s, boost::regex(pattern))) {
      filenames.push_back(s);
    }
  }
  if (filenames.empty())
    throw Error("WHOOPS : " + file_or_dir + " does not contain any files...");
  return filenames;
}

}; // namespace inout
