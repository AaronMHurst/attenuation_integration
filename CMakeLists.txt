# Setup the project
cmake_minimum_required(VERSION 3.10)
project(attenuation_integration)

# Define the application executable
set(APP_NAME attenuation)

# C++ compilation options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_COMPILER_ID "g++")
#set(CMAKE_CXX_COMPILER_ID "Clang")
set(CMAKE_CXX_FLAGS "-g -O3 -Wall -Wextra -Wpedantic")
#set(CMAKE_CXX_FLAGS "-g -O3 -Wall -Wpedantic")

# Locate the sources and headers for this project
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)
include_directories(${PROJECT_SOURCE_DIR}/include)

# Copy the data folders to the project-build directory
file(COPY ${PROJECT_SOURCE_DIR}/mass_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/phys_data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Add the executable
add_executable(${APP_NAME} ${SOURCES})

# Copy the example scripts to the build directory
set(EXAMPLE_SCRIPTS
  gamma_attn_tungsten.mac
  gamma_attn_water.mac
  )

foreach(_script ${EXAMPLE_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script} ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
endforeach()
    
