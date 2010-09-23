cmake_minimum_required(VERSION 2.6)

Project(KMeansClustering)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(KMeansClustering vtkKMeansClustering.cxx Example.cxx)
TARGET_LINK_LIBRARIES(KMeansClustering vtkHybrid)
