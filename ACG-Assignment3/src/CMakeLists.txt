cmake_minimum_required(VERSION 2.8)
project(Assignment3)

link_directories("/usr/lib/x86_64-linux-gnu")
link_libraries("libGl.so")
link_libraries("libGlut.so")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(SOURCE_FILES "main.cpp")

add_executable(Assignment3 ${SOURCE_FILES})