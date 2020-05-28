cmake_minimum_required(VERSION 2.8)

project(VisualSPHysics)

set(CMAKE_SHARED_LIBRARY_PREFIX "")

if (MSVC)
	add_definitions(-D_MSVC)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

add_subdirectory(foamsimulator)
add_subdirectory(vtkimporter)

INSTALL(FILES ${PROJECT_SOURCE_DIR}/blendermodule/VisualSPHysics.py DESTINATION ".")

SET(CPACK_GENERATOR "ZIP")
SET(CPACK_ARCHIVE_FILE_NAME "visualsphysics")
SET(CPACK_INCLUDE_TOPLEVEL_DIRECTORY OFF)

INCLUDE(CPack)
