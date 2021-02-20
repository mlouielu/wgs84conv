cmake_minimum_required(VERSION 3.15)
project(wgs84conv LANGUAGES CXX)

add_subdirectory(pybind11)
pybind11_add_module(wgs84conv wgs84conv.cpp wrapper.cpp)

target_include_directories(wgs84conv PUBLIC carma/include)
