cmake_minimum_required(VERSION 3.21)

# vcpkg toolchain should be set before project() when provided via environment.
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    if(DEFINED ENV{VCPKG_ROOT})
        set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
    endif()
endif()

project(HeatEquationVTK CXX)

# Prefer static triplet, but allow fallback if only dynamic packages are installed locally.
set(VCPKG_TARGET_TRIPLET "x64-windows-static" CACHE STRING "vcpkg target triplet")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# vcpkg package directory - auto-detect or use environment variable
if(DEFINED ENV{VCPKG_INSTALLED_DIR})
    set(VCPKG_INSTALLED_DIR "$ENV{VCPKG_INSTALLED_DIR}")
elseif(EXISTS "${CMAKE_SOURCE_DIR}/vcpkg_installed")
    set(VCPKG_INSTALLED_DIR "${CMAKE_SOURCE_DIR}/vcpkg_installed")
endif()

if(DEFINED VCPKG_INSTALLED_DIR)
    set(_vcpkg_triplet_path "${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}")

    if(EXISTS "${_vcpkg_triplet_path}")
        list(APPEND CMAKE_PREFIX_PATH "${_vcpkg_triplet_path}")
        message(STATUS "vcpkg packages found at: ${_vcpkg_triplet_path}")
    elseif(EXISTS "${VCPKG_INSTALLED_DIR}/x64-windows")
        message(WARNING "Requested triplet '${VCPKG_TARGET_TRIPLET}' not found under local vcpkg_installed. Falling back to x64-windows.")
        set(VCPKG_TARGET_TRIPLET "x64-windows" CACHE STRING "vcpkg target triplet" FORCE)
        list(APPEND CMAKE_PREFIX_PATH "${VCPKG_INSTALLED_DIR}/x64-windows")
    endif()
endif()

find_package(Eigen3 REQUIRED CONFIG)
find_package(pugixml REQUIRED CONFIG)

add_executable(heat_equation heat_solver.h heat_solver.cpp heat_solver_main.cpp vtk_writer.h vtk_writer.cpp)

target_link_libraries(heat_equation PRIVATE Eigen3::Eigen pugixml::pugixml)

if(WIN32)
    add_custom_command(TARGET heat_equation POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_RUNTIME_DLLS:heat_equation>
                $<TARGET_FILE_DIR:heat_equation>
        COMMAND_EXPAND_LISTS
        COMMENT "Copy runtime DLLs next to heat_equation"
    )
endif()
