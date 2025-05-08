cmake_minimum_required(VERSION 3.20)
project(spoc)

# Enable testing with CTest
enable_testing()

# Add subdirectories
add_subdirectory(tools)

message(STATUS "TOP CMAKE_BINARY_DIR : ${CMAKE_BINARY_DIR}")
message(STATUS "TOP CMAKE_SOURCE_DIR : ${CMAKE_SOURCE_DIR}")
message(STATUS "TOP CMAKE_CURRENT_SOURCE_DIR : ${CMAKE_CURRENT_SOURCE_DIR}")

# Register the Python style check as a test
add_test(NAME ${PROJECT_NAME}_python_coding_norms 
         COMMAND ${CMAKE_BINARY_DIR}/bin/spoc_py_lint.sh
	         ${CMAKE_SOURCE_DIR}/dump/mapping
                 ${CMAKE_SOURCE_DIR})

