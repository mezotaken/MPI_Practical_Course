project( scalar_multiplication )
message( STATUS "-- " ${PROJECT_NAME} )

add_executable( ${PROJECT_NAME} main.cpp )

if (UNIX)
    SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11")
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

if( MPI_COMPILE_FLAGS )
    set_target_properties( ${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}" )
endif( MPI_COMPILE_FLAGS )

if( MPI_LINK_FLAGS )
    set_target_properties( ${PROJECT_NAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}" )
endif( MPI_LINK_FLAGS )

target_link_libraries( ${PROJECT_NAME} ${MPI_LIBRARIES} )
