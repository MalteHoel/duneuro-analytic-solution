# File for module specific CMake tests.
find_package(PythonLibs REQUIRED)
message("CMAKE_SOURCE_DIR : ${CMAKE_SOURCE_DIR}")
include_directories(${PYTHON_INCLUDE_DIR})
