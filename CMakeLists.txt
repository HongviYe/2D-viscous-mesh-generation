cmake_minimum_required(VERSION 3.16)
project(expansion)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
add_subdirectory(third/triangle)
add_subdirectory(third/dtiso2d)
add_subdirectory(third/geom)
# Libigl
option(LIBIGL_GLFW "Build target igl::glfw" ON)
include(libigl)

# Add your project files
file(GLOB INC_FILES  ./include/*.h)
file(GLOB SRC_FILES  ./src/*.cpp)
include_directories(third/triangle)
include_directories(third/dtiso2d)
include_directories(third/geom)
include_directories(third/CLI11/include)

set(ALL_FILES 
${INC_FILES}
${SRC_FILES} 
)
add_executable(${PROJECT_NAME} ${ALL_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw triangle libdtiso2d libgeom) 
