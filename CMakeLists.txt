cmake_minimum_required(VERSION 2.6)

Project(KMeansClustering)
ENABLE_TESTING()

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(KMeansClustering Example.cxx vtkKMeansClustering.cxx)
TARGET_LINK_LIBRARIES(KMeansClustering ${VTK_LIBRARIES})

## test
ADD_EXECUTABLE(Test Test.cxx vtkKMeansClustering.cxx )
TARGET_LINK_LIBRARIES(Test ${VTK_LIBRARIES})

ADD_TEST(Test ${CMAKE_CURRENT_BINARY_DIR}/Test)