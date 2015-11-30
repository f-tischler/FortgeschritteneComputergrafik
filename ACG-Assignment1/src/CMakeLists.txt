cmake_minimum_required(VERSION 2.8)
project(Assignment1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
set(SOURCE_FILES "Radiosity.cpp")

add_executable(Assignment1 ${SOURCE_FILES})