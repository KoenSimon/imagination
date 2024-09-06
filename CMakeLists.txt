cmake_minimum_required(VERSION 3.16)
project(KNN_Server)
#set(NAME "KNN_Server")
#find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
add_compile_options(-march=native)
#find_package(emp-ot REQUIRED)

#INCLUDE_DIRECTORIES(${EMP-OT_INCLUDE_DIRS})

#add_executable(KNN_Server src/main.cpp src/crypto.h src/crypto.cpp src/SetupPhase.cpp src/SetupPhase.h src/readData.cpp src/readData.h src/MyTime.cpp src/MyTime.h src/KNN_single.cpp src/KNN_single.h)
add_executable(KNN_Server src/main.cpp src/crypto.cpp src/SetupPhase.cpp src/readData.cpp src/MyTime.cpp src/KNN_single.cpp src/Poly.cpp src/Poly.h src/Poly_c.cpp src/Poly_c.h src/polynomial.h src/spec_sort.h)

include_directories("/usr/lib/x86_64-linux-gnu/")

target_link_libraries(KNN_Server gmp gmpxx)

target_link_libraries (${PROJECT_NAME} pthread)
#target_link_libraries(KNN_Server emp-ot)