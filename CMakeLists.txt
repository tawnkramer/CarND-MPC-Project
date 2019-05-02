project(MPC)

cmake_minimum_required (VERSION 3.5)

# -g allows for gdb debugging
# turn on -03 for best performance
add_definitions(-std=c++11 -O3)

set(CXX_FLAGS "-Wall -Wno-int-in-bool-context")
set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")

#set(CMAKE_BUILD_TYPE Debug)

set(sources src/MPC.cpp src/main.cpp)
set(sources_z src/MPC.cpp src/main_z.cpp)

include_directories(/usr/local/include)
link_directories(/usr/local/lib)
include_directories(src/Eigen-3.3)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1*/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

add_executable(mpc ${sources})
add_executable(mpc_zmq ${sources_z})

target_link_libraries(mpc ipopt z ssl uv uWS)
target_link_libraries(mpc_zmq ipopt z uv zmq)

