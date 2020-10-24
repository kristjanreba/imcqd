include_directories(graph)

set( graph_SRC
  graph_mcqd.hpp
  mcqd.cpp
  mcqd.hpp
  vertex.hpp
  )

add_library(graph STATIC ${graph_SRC})

