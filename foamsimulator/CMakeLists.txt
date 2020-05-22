cmake_minimum_required(VERSION 2.8)
 
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(VTK COMPONENTS vtkFiltersModeling vtkIOGeometry)
include(${VTK_USE_FILE})

set(SRCS FluidData.cpp VtkDWriter.cpp Ops.cpp DiffuseCalculator.cpp diffuseparticlesmodule.cpp)
 
add_library(diffuseparticles SHARED ${SRCS})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if (NOT WIN32) 
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -fopenmp-simd -I${PYTHON_INCLUDE_DIRS}")
endif ()
 
if(VTK_LIBRARIES)
  target_link_libraries(diffuseparticles ${VTK_LIBRARIES} ${PYTHON_LIBRARIES})
else()
  target_link_libraries(diffuseparticles vtkHybrid ${PYTHON_LIBRARIES})
endif()

if (WIN32)
  set_target_properties(diffuseparticles PROPERTIES SUFFIX ".pyd")
endif()

INSTALL(TARGETS diffuseparticles DESTINATION ".")

INCLUDE(CPack)