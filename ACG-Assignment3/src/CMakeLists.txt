cmake_minimum_required(VERSION 2.8)
project(Assignment3)

if(APPLE)
include_directories("/System/Library/Frameworks")
find_library(OpenGL_LIBRARY OpenGL)
link_libraries("/usr/local/Cellar/freeglut/2.8.1/lib/libglut.dylib")
set(EXTRA_LIBS ${OpenGL_LIBRARY})
else()
link_directories("/usr/lib/x86_64-linux-gnu")
link_libraries("libGl.so")
link_libraries("libGlut.so")
endif(APPLE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(SOURCE_FILES "main.cpp")

add_executable(Assignment3 ${SOURCE_FILES})
target_link_libraries(Assignment3 ${EXTRA_LIBS})