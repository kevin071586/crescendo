CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

PROJECT( MyApp )
SET(CMAKE_BUILD_TYPE Release )

# Get Trilinos as one entity
FIND_PACKAGE(Trilinos PATHS /home/kevin/code/my_trilinos/install/lib/cmake/ NO_DEFAULT_PATH REQUIRED)

# Echo trilinos build info just for fun
 MESSAGE("\nFound Trilinos!  Here are the details: ")
 MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
 MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
 MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
 MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
 MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
 MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
 MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
# MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
# MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
# MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
 MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
 MESSAGE("End of Trilinos details\n")

#------------------------------------
#  Make sure to use same compilers 
#  and flags as Trilinos
#------------------------------------

# Compiler
SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER} )
SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER} )
SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER} )

# Compiler flags
SET(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
SET(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
SET(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

#------------------------------------
#  Define:
#    1) Include files
#    2) Library files
#    3) Linking instructions
#    4) Executables
#------------------------------------
SET(CRESCENDO_DIR "/home/kevin/code/crescendo")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CRESCENDO_DIR}/bin)

INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS} 
  ${Trilinos_TPL_INCLUDE_DIRS}
  ${CRESCENDO_DIR}/include
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Print out include directories:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
# foreach(dir ${dirs})
#   message(STATUS "dir='${dir}'")
# endforeach()

LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

LINK_LIBRARIES(${Trilinos_LIBRARIES})

ADD_EXECUTABLE( crescendo 
  ${CRESCENDO_DIR}/src/crescendo.cpp
  ${CRESCENDO_DIR}/src/integrateMatrices.cpp
  )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Some boiler plate for future reference
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ADD_DEFINITIONS( "-std=c++11" ) 
#INCLUDE_DIRECTORIES( /home/kevin/fea/my_trilinos_stk/install/include/ )
#SET( TRILINOS_LIB /home/kevin/fea/my_trilinos_stk/install/lib/ )
#find_library( EXODUS_LIB NAMES exodus PATHS ${TRILINOS_LIB} )
#find_library( KOKKOSCORE_LIB NAMES kokkoscore PATHS ${TRILINOS_LIB} )
#find_library( TEUCHOSCOMM_LIB NAMES teuchoscomm PATHS ${TRILINOS_LIB} )
#find_library( TEUCHOSCORE_LIB NAMES teuchoscore PATHS ${TRILINOS_LIB} )
#find_library( TEUCHOSNUMERICS_LIB NAMES teuchosnumerics PATHS ${TRILINOS_LIB} )
#find_library( INTREPID_LIB NAMES intrepid PATHS ${TRILINOS_LIB} )
#find_library( SHARDS_LIB NAMES shards  PATHS ${TRILINOS_LIB} )
#LINK_LIBRARIES( netcdf )
#LINK_LIBRARIES( ${EXODUS_LIB} )
#LINK_LIBRARIES( ${KOKKOSCORE_LIB} )
#LINK_LIBRARIES( ${TEUCHOSCOMM_LIB} )
#LINK_LIBRARIES( ${TEUCHOSCORE_LIB} )
#LINK_LIBRARIES( ${TEUCHOSNUMERICS_LIB} )
#LINK_LIBRARIES( ${INTREPID_LIB} )
#LINK_LIBRARIES( ${SHARDS_LIB} )


