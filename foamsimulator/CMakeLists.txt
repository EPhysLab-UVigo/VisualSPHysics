cmake_minimum_required(VERSION 2.8)
 
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(VTK COMPONENTS vtkFiltersModeling vtkIOGeometry)

if (VTK_VERSION VERSION_GREATER_EQUAL "8.90.0")
  add_definitions(-DVTK9)
endif()

include(${VTK_USE_FILE})

set(SRCS FluidData.cpp VtkDWriter.cpp Ops.cpp DiffuseCalculator.cpp diffuseparticlesmodule.cpp)
 
add_library(diffuseparticles SHARED ${SRCS})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -fopenmp-simd -I${PYTHON_INCLUDE_DIRS}")
   set(CXX_LDFLAGS "-lstdc++fs")
endif()
 
if(VTK_LIBRARIES)
  target_link_libraries(diffuseparticles ${VTK_LIBRARIES} ${PYTHON_LIBRARIES} ${CXX_LDFLAGS})
else()
  target_link_libraries(diffuseparticles vtkHybrid ${PYTHON_LIBRARIES} ${CXX_LDFLAGS})
endif()

if (WIN32)
  set_target_properties(diffuseparticles PROPERTIES SUFFIX ".pyd")
endif()

INSTALL(TARGETS diffuseparticles DESTINATION ".")

INCLUDE(CPack)