# 设置编译的源文件
set(LIBGEOM_SRC 
		intersection_func.cpp
		geom_func.cpp
		geom_func.h
		predicates.c
		ttotr.c  )

#AUX_SOURCE_DIRECTORY(. DIR_SRCS)

# 定义宏
add_definitions(
    #-D_LIB
    -D_OPT_INTERSEC_ALGO
    #-D_USE_STATIC_FILTER
    )

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffloat-store")

# 设置生成库
if(WIN32)  # 暂时windows生成静态库
    add_library(libgeom ${LIBGEOM_SRC})
else()
if(SHARED)
    add_library(libgeom SHARED ${LIBGEOM_SRC})
else()
    add_library(libgeom STATIC ${LIBGEOM_SRC})
endif()
endif()

# 设置库的输出路径
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# 设置库的输出名字
set_target_properties(libgeom PROPERTIES OUTPUT_NAME "geom")

install(TARGETS libgeom
	ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    )
