

list(APPEND _scripts
  spoc_py_lint.sh)

message(STATUS "TOOLS CMAKE_BINARY_DIR : ${CMAKE_BINARY_DIR}")
message(STATUS "TOOLS CMAKE_CURRENT_SOURCE_DIR_DIR : ${CMAKE_CURRENT_SOURCE_DIR}")

foreach(_f IN LISTS _scripts)
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/${_f} ${CMAKE_BINARY_DIR}/bin/${_f} )
endforeach()
