include_directories(helper)

set( helper_SRC 
  array1d.hpp
  array2d.cpp
  array2d.hpp
  array3d.hpp
  benchmark.cpp
  benchmark.hpp
  debug.hpp
  error.hpp
  help.cpp
  help.hpp
  inout.cpp
  inout.hpp
  path.cpp
  path.hpp 
  iomanip.cpp
  iomanip.hpp 
  graph.cpp
  graph.hpp
  util.cpp
  util.hpp
  )

add_library(helper STATIC ${helper_SRC})
