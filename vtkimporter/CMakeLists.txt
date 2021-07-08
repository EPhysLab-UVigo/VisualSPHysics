cmake_minimum_required(VERSION 2.8)
 
find_package(PythonLibs 3.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(VTK COMPONENTS vtkIOLegacy)

if (VTK_VERSION VERSION_GREATER_EQUAL "8.90.0")
  add_definitions(-DVTK9)
endif()

include(${VTK_USE_FILE})

set(SRCS vtkimportermodule.cpp)
 
add_library(vtkimporter SHARED ${SRCS})

if (NOT WIN32)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math -std=c++11 -I${PYTHON_INCLUDE_DIRS}")
endif ()

if(VTK_LIBRARIES)
  target_link_libraries(vtkimporter ${PYTHON_LIBRARIES} ${VTK_LIBRARIES})
else()
  target_link_libraries(vtkimporter ${PYTHON_LIBRARIES} vtkHybrid )
endif()

if (WIN32)
  set_target_properties(vtkimporter PROPERTIES SUFFIX ".pyd")
endif()

INSTALL(TARGETS vtkimporter DESTINATION ".")

INCLUDE(CPack)