cmake_minimum_required(VERSION 2.8)
 
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(SRCS FluidData.cpp PlyWriter.cpp VtkDWriter.cpp Ops.cpp DiffuseCalculator.cpp rply/rply.c diffuseparticlesmodule.cpp)
 
#add_executable(diffuse_parts ${SRCS})

add_library(diffuseparticles SHARED ${SRCS})

add_subdirectory( rply )

#set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if (NOT WIN32) 
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -ffast-math -std=c++11 -I${PYTHON_INCLUDE_DIRS}")
endif ()
 
if(VTK_LIBRARIES)
  target_link_libraries(diffuseparticles ${VTK_LIBRARIES} ${PYTHON_LIBRARIES})
else()
  target_link_libraries(diffuseparticles vtkHybrid ${PYTHON_LIBRARIES})
endif()

