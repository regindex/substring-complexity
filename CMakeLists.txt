cmake_minimum_required(VERSION 2.6)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (substring-complexity)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(~/include) #SDSL headers are here
INCLUDE_DIRECTORIES(include/cpp-HyperLogLog/include include/cpp-HyperLogLog/extlib/igloo include/cpp-HyperLogLog/extlib/igloo-TapTestListener)

LINK_DIRECTORIES(~/lib) #SDSL lib are here

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -Wno-deprecated")
set(CMAKE_CXX_FLAGS_RELEASE "-g -ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native -Wno-deprecated")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native -Wno-deprecated")

add_executable(delta delta.cpp)
add_executable(delta-stream delta-stream.cpp)
TARGET_LINK_LIBRARIES(delta sdsl)
TARGET_LINK_LIBRARIES(delta-stream sdsl)
TARGET_LINK_LIBRARIES(delta divsufsort)
TARGET_LINK_LIBRARIES(delta divsufsort64)