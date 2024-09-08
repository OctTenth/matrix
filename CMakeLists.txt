cmake_minimum_required(VERSION 3.10)

project(Matrix)

set(CMAKE_CXX_STANDARD 20)

file(GLOB SRC_FILE
          "./src/*.cpp"
          "./src/*.hpp"
          "./src/*.tpp"
          )

add_executable(Matrix ${SRC_FILE})