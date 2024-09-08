cmake_minimum_required(VERSION 3.10)

project(Matrix)

file(GLOB SRC_FILE
          "./src/*.cpp"
          "./src/*.hpp"
          "./src/*.tpp"
          )

add_executable(Matrix ${SRC_FILE})