cmake_minimum_required(VERSION 3.8)


PROJECT(FT)
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.1)
else()
  MESSAGE(FATAL_ERROR "gcc version must bigger than 5.1. Current version is ${CMAKE_CXX_COMPILER_VERSION}")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
            "Default build type: Debug" FORCE)
endif()

# Set default path
# set (VTK_DIR CACHE PATH "~/data/vtk-build")
# set (HDF5_INCLUDE_DIR "~/data/hdf5")
# set (NETCDF_INCLUDE_DIR "~/usr/Netcdf")
# set(OpenCV_DIR "~/usr/OpenCV")

# Set hdf5 path
FIND_PACKAGE (HDF5)
if(HDF5_FOUND)
  include_directories (${HDF5_INCLUDE_DIR})
  set(LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
  message(STATUS "Hdf5 is found: ${HDF5_INCLUDE_DIR}")
else()
  if(EXISTS ${HDF5_PATH})
    set(HDF5_DIR ${HDF5_PATH})
    FIND_PACKAGE(HDF5 REQUIRED)
    if(NOT HDF5_FOUND)
      message(FATAL_ERROR "hdf5 not found. Please provide the path of hdf5")
    else()
      include_directories (${HDF5_INCLUDE_DIR})
      set(LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
      message(STATUS "Hdf5 is found: ${HDF5_INCLUDE_DIR}")
    endif()
  else()
    message(FATAL_ERROR "hdf5 not found. Please provide the path of hdf5")
    set (HDF5_INCLUDE_PATH ""  CACHE PATH "Path to hdf5")
  endif()
endif()



FIND_PACKAGE (VTK)
if(VTK_FOUND)
  MESSAGE(STATUS "VTK_DIR is found: ${VTK_DIR}")
  INCLUDE(${VTK_USE_FILE})
else()
  if(EXISTS ${VTK_PATH})
    set(VTK_DIR ${VTK_PATH})
    find_package(VTK REQUIRED)
    if(NOT VTK_FOUND)
      message(FATAL_ERROR "VTK not found. P")
    else()
      MESSAGE(STATUS "VTK_DIR is found: ${VTK_DIR}")
    ENDIF()
  else()
    MESSAGE(FATAL_ERROR "VTK not found. Please provide the path of VTK_DIR.")
    set (VTK_PATH "" CACHE PATH "Path to vtk library")
  endif() 
endif()
INCLUDE(${VTK_USE_FILE})


find_package(OpenCV)
if(OpenCV_FOUND)
  message(STATUS "OpenCV is found: ${OpenCV_DIR}")
else()
  if(EXISTS ${OpenCV_PATH})
    set(OpenCV_DIR ${OpenCV_PATH})
    find_package(OpenCV REQUIRED)
    if(NOT OpenCV_FOUND)
      message(FATAL_ERROR FATAL_ERROR "OpenCV not found. Please provide the path of OpenCV.")
    else()
      messsage(STATUS "OpenCV is found: ${OpenCV_DIR}")
    endif()
  else()   
    message(FATAL_ERROR "OpenCV not found. Please provide the path of OpenCV.")
    set(OPENCV_PATH "" CACHE PATH "Path to OpenCV")
  endif()
endif()





# set(CMAKE_LIBRARY_PATH /home/ardofski/netcdf_build/CMakeFiles/Export/lib)
# set(CMAKE_INCLUDE_PATH /home/ardofski/netcdf_build/include)


function(FindNetCDF_get_is_parallel_aware include_dir)
  file(STRINGS "${include_dir}/netcdf_meta.h" _netcdf_lines
    REGEX "#define[ \t]+NC_HAS_PARALLEL[ \t]")
  string(REGEX REPLACE ".*NC_HAS_PARALLEL[ \t]*([0-1]+).*" "\\1" _netcdf_has_parallel "${_netcdf_lines}")
  if (_netcdf_has_parallel)
    set(NetCDF_HAS_PARALLEL TRUE PARENT_SCOPE)
  else()
    set(NetCDF_HAS_PARALLEL FALSE PARENT_SCOPE)
  endif()
endfunction()

# Try to find a CMake-built NetCDF.
find_package(netCDF CONFIG QUIET)
if (netCDF_FOUND)
  # Forward the variables in a consistent way.
  set(NetCDF_FOUND "${netCDF_FOUND}")
  set(NetCDF_INCLUDE_DIRS "${netCDF_INCLUDE_DIR}")
  set(NetCDF_LIBRARIES "${netCDF_LIBRARIES}")
  set(NetCDF_VERSION "${NetCDFVersion}")

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(NetCDF
    REQUIRED_VARS NetCDF_INCLUDE_DIRS NetCDF_LIBRARIES
    VERSION_VAR NetCDF_VERSION)

  if (NOT TARGET NetCDF::NetCDF)
    add_library(NetCDF::NetCDF INTERFACE IMPORTED)
    if (TARGET "netCDF::netcdf")
      # 4.7.3
      set_target_properties(NetCDF::NetCDF PROPERTIES
        INTERFACE_LINK_LIBRARIES "netCDF::netcdf")
    elseif (TARGET "netcdf")
      set_target_properties(NetCDF::NetCDF PROPERTIES
        INTERFACE_LINK_LIBRARIES "netcdf")
    else ()
      set_target_properties(NetCDF::NetCDF PROPERTIES
        INTERFACE_LINK_LIBRARIES "${netCDF_LIBRARIES}")
    endif ()
  endif ()

  FindNetCDF_get_is_parallel_aware("${NetCDF_INCLUDE_DIRS}")
  # Skip the rest of the logic in this file.
  return ()
endif ()




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
  ${NETCDF_LIBRARIES}
  ${OpenCV_LIBS}
)



