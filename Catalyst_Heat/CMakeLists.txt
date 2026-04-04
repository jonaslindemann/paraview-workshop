cmake_minimum_required(VERSION 3.20)
project(HeatEquationCatalyst CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------------------
# Catalyst SDK
#
# The Catalyst SDK is NOT bundled with the ParaView binary installer.
# You need to build it from source:
#
#   git clone https://gitlab.kitware.com/paraview/catalyst.git
#   cmake -B build -S catalyst -DCATALYST_BUILD_STUB_IMPLEMENTATION=ON
#   cmake --build build --config Release
#   cmake --install build --prefix C:/catalyst-sdk
#
# Then point cmake here:
#   cmake -B build -DCATALYST_DIR=C:/catalyst-sdk/lib/cmake/catalyst-2.x
#
# At runtime, place catalyst-paraview.dll (from your ParaView install) next
# to the executable, or set CATALYST_IMPLEMENTATION_PATHS to point to it.
# ---------------------------------------------------------------------------
find_package(catalyst REQUIRED
    HINTS
        "${CATALYST_DIR}"
        "C:/catalyst-sdk/lib/cmake/catalyst-2.0"
        "$ENV{CATALYST_DIR}"
)

# ---------------------------------------------------------------------------
# Eigen3 (header-only linear algebra)
#
# Install via vcpkg:  vcpkg install eigen3
# Or set EIGEN3_DIR:  cmake -B build -DEIGEN3_DIR=C:/vcpkg/packages/eigen3_x64-windows/share/eigen3
# ---------------------------------------------------------------------------
find_package(Eigen3 REQUIRED NO_MODULE
    HINTS
        "${EIGEN3_DIR}"
        "C:/vcpkg/packages/eigen3_x64-windows/share/eigen3"
        "$ENV{VCPKG_ROOT}/packages/eigen3_x64-windows/share/eigen3"
)

# ---------------------------------------------------------------------------
# Simulation executable
# ---------------------------------------------------------------------------
add_executable(heat_equation heat_equation.cpp)

target_link_libraries(heat_equation
    PRIVATE
        catalyst::catalyst
        Eigen3::Eigen
)

# Copy catalyst DLLs next to the executable on Windows so the loader finds them
if(WIN32)
    add_custom_command(TARGET heat_equation POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:catalyst::catalyst>
            $<TARGET_FILE_DIR:heat_equation>
        COMMENT "Copying catalyst stub DLL"
    )
    # Print a reminder about the runtime implementation DLL
    message(STATUS
        "REMINDER: copy catalyst-paraview.dll from your ParaView install/bin "
        "into the build output directory before running."
    )
endif()
