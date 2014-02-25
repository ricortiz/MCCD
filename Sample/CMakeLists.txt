
add_library(loader ../Geometry/plyfile.cpp ../Srcs/loader.cpp)

add_executable(pmain pmain.cpp)
target_link_libraries(pmain mccd loader GL)
