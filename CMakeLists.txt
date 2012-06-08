cmake_minimum_required(VERSION 2.6)

Project(KMeansClustering)
ENABLE_TESTING()

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

if( "${VTK_MAJOR_VERSION}" LESS 6 )
  MESSAGE(FATAL_ERROR "You must build this code with VTK >= 6.0!")
endif( "${VTK_MAJOR_VERSION}" LESS 6 )

add_library(libKMeansClustering vtkKMeansClustering.cxx)
TARGET_LINK_LIBRARIES(libKMeansClustering ${VTK_LIBRARIES})

ADD_EXECUTABLE(KMeansClustering Example.cxx)
TARGET_LINK_LIBRARIES(KMeansClustering ${VTK_LIBRARIES} libKMeansClustering)

## test
ADD_EXECUTABLE(Test Test.cxx)
TARGET_LINK_LIBRARIES(Test ${VTK_LIBRARIES} libKMeansClustering)

ADD_TEST(Test ${CMAKE_CURRENT_BINARY_DIR}/Test)