cmake_minimum_required(VERSION 3.10)
project(mahalonobis)

set(CMAKE_CXX_STANDARD 17)

add_executable(mahalonobis main.cpp)
target_link_libraries(mahalonobis m)
