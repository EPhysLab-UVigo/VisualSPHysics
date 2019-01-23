cmake_minimum_required(VERSION 2.8)

project(VisualSPHysics)

set(CMAKE_SHARED_LIBRARY_PREFIX "")

if (MSVC)
	add_definitions(-D_MSVC)
endif()

add_subdirectory(foamsimulator)
add_subdirectory(vtkimporter)
