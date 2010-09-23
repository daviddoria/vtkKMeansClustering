cmake_minimum_required(VERSION 2.6)

Project(AgglomerativeClustering)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(AgglomerativeClustering vtkAgglomerativeClustering.cxx Example.cxx)
TARGET_LINK_LIBRARIES(AgglomerativeClustering vtkHybrid)
