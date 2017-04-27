cmake_minimum_required(VERSION 3.6)
project(mixed_code)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        cdn/lib/lib_io.h
        cdn/lib/lib_time.h
        cdn/lib/SKW.h.h
        cdn/cdn.cpp
        cdn/deploy.cpp
        cdn/deploy.h
        cdn/io.cpp
        cdn/MCMF.cpp
        cdn/MCMF.h
        cdn/minmax.h)

add_executable(mixed_code ${SOURCE_FILES})