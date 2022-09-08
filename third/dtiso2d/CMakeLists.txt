cmake_minimum_required (VERSION 2.6)
project(dtiso2d)
#best opt
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
#
message("building libdtiso2d")
#
add_definitions(-D _LIB)

add_library(libdtiso2d STATIC
genmesh2d.cpp
genmesh2d.h
iso2d.cpp
iso2d.h
iso2d_error.cpp
iso2d_error.h
../geom/geom_func.h
)
set_target_properties(libdtiso2d PROPERTIES OUTPUT_NAME dtiso2d)
# 设置库的输出路径
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# 设置可执行程序的输出路径
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

install(TARGETS libdtiso2d
	RUNTIME DESTINATION bin
	ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    )