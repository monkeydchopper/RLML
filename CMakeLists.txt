cmake_minimum_required(VERSION 3.10)
project(RLML)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Ensure C++17 is required
set(CMAKE_CXX_EXTENSIONS OFF)        # Use only standard-compliant C++17


if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build" FORCE)
endif()
if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
        set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)


find_package(catkin REQUIRED COMPONENTS
        roscpp
        tf
        sensor_msgs
        geometry_msgs
        nav_msgs
        visualization_msgs
        tf_conversions
        pcl_conversions
        )
find_package(Cholmod REQUIRED)
find_package(PCL REQUIRED COMPONENTS common io)


catkin_package(
        CATKIN_DEPENDS roscpp tf sensor_msgs geometry_msgs nav_msgs visualization_msgs tf_conversions pcl_conversions
)



include_directories(
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/include/RLML
        /usr/include/eigen3
        ${PROJECT_SOURCE_DIR}/Thirdparty
        ${PROJECT_SOURCE_DIR}/Thirdparty/g2o
        ${PCL_INCLUDE_DIRS}
        ${catkin_INCLUDE_DIRS}
        ${CHOLMOD_INCLUDE_DIR}
)

add_executable(RLML_node
        src/utility.cpp
        src/node.cpp
        src/RLML_ros.cpp
        src/map_builder.cpp
        src/correlative_scan_matcher.cpp
        src/gauss_newten_scan_matcher.cpp
        src/correlative_grid.cpp
        src/occupancy_grid_map.cpp
        src/probability_grid_map.cpp
        src/landmark.cpp
        )



target_link_libraries(RLML_node
        ${catkin_LIBRARIES}
        ${CHOLMOD_LIBRARIES}
        ${PCL_COMMON_LIBRARIES}
        ${PCL_IO_LIBRARIES}
        ${PROJECT_SOURCE_DIR}/Thirdparty/g2o/lib/libg2o_core.so
        ${PROJECT_SOURCE_DIR}/Thirdparty/g2o/lib/libg2o_solver_cholmod.so
        ${PROJECT_SOURCE_DIR}/Thirdparty/g2o/lib/libg2o_types_slam2d.so
        ${PROJECT_SOURCE_DIR}/Thirdparty/g2o/lib/libg2o_stuff.so
        )

