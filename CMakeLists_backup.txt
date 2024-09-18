cmake_minimum_required(VERSION 2.6)

PROJECT(FT)
SET (CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)



set (VTK_DIR "~/data/vtk-build")
set (OpenCV_DIR "~/data/openCV/build")
set (NETCDF_INCLUDE_DIR "~/usr/include")



IF(NOT VTK_BINARY_DIR)
FIND_PACKAGE(VTK COMPONENTS)
IF(NOT VTK_DIR)
  MESSAGE(FATAL_ERROR "Please set VTK_DIR.")
ELSE(NOT VTK_DIR)
  MESSAGE(STATUS "VTK_DIR is found: ${VTK_DIR}")
ENDIF(NOT VTK_DIR)
INCLUDE(${VTK_USE_FILE})
ENDIF(NOT VTK_BINARY_DIR)
find_package(OpenCV REQUIRED)

FIND_PACKAGE (HDF5)
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
SET (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})



set(CMAKE_LIBRARY_PATH /home/ardofski/netcdf_build/CMakeFiles/Export/lib)
set(CMAKE_INCLUDE_PATH /home/ardofski/netcdf_build/include)
find_path(NETCDF_INCLUDE_DIR netcdf.mod DOC "NetCDF include directory (must contain netcdf.mod)")
if(NETCDF_INCLUDE_DIR)
  MESSAGE(STATUS "NETCDF module found in " ${NETCDF_INCLUDE_DIR})
else()
  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")
endif()


find_library(NETCDF_C_LIB netcdf DOC "NetCDF C library" PATHS $ENV{NETCDF_LIB} /usr/lib /usr/lib64 /usr/include $ENV{SARA_NETCDF_LIB})

if(NETCDF_C_LIB)
  set(NETCDF_LIBS ${NETCDF_C_LIB})
  MESSAGE(STATUS "NETCDF C library: " ${NETCDF_C_LIB})
else()
  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")
endif()

if(NETCDF_INCLUDE_DIR)

else(NETCDF_INCLUDE_DIR)

  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")

endif(NETCDF_INCLUDE_DIR)

if(NETCDF_C_LIB)

  set(NETCDF_LIBS ${NETCDF_C_LIB})

else(NETCDF_C_LIB)

  MESSAGE(STATUS "WARNING: No NETCDF bindings are found.")

endif(NETCDF_C_LIB)




# Use the include path and library for Qt that is used by VTK.
INCLUDE_DIRECTORIES(
  ${HDF5_INCLUDE_DIR}
  ${NETCDF_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)


# Set your files and resources here  
SET(FTSrcs cellinfo.cpp classifier.cpp FeatureTrack.cpp mainFeatureTrack.cpp nodeinfo.cpp objectinfo.cpp)
SET(FTHeaders cellinfo.h classifier.h FeatureTrack.h FeatureTrackUtil.h input.h interfaceutil.h nodeinfo.h objectinfo.h ObjectSegment.h objSegmentUtil.h stobject1.h stRGB.h util.h)


SOURCE_GROUP("Resources" FILES
  ${FTResources}
  ${EXE_ICON}
)

SOURCE_GROUP("Generated" FILES
  ${UISrcs}
  ${MOCSrcs}
  ${ResourceSrcs}
  ${RCS_SOURCES}
)

ADD_DEFINITIONS(-DQT_GUI_LIBS -DQT_CORE_LIB -DQT3_SUPPORT)
SET_SOURCE_FILES_PROPERTIES(${FTSrcs} PROPERTIES)


ADD_EXECUTABLE( FT MACOSX_BUNDLE ${FTSrcs} ${FTResources} ${MOCSrcs} ${ResourceSrcs} )




TARGET_LINK_LIBRARIES( FT
  ${VTK_LIBRARIES}
  #vtkRendering
  #vtkGraphics
  #vtkIO
  #vtkCommon
  #vtkInfovis
  #vtkViews
  #vtkHybrid
  ${HDF5_LIBRARIES}
  ${NETCDF_LIBS}
  ${OpenCV_LIBS}
)



