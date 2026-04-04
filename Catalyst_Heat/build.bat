@echo off
rem --------------------------------------------------------------------
rem Build script for HeatEquationCatalyst (Windows / MSVC)
rem
rem Prerequisites
rem -------------
rem 1. Catalyst SDK (build from source — NOT included in ParaView installer)
rem
rem    git clone https://gitlab.kitware.com/paraview/catalyst.git catalyst-src
rem    cmake -B catalyst-build -S catalyst-src ^
rem          -DCATALYST_BUILD_STUB_IMPLEMENTATION=ON ^
rem          -DCMAKE_INSTALL_PREFIX=C:/catalyst-sdk ^
rem          -DCMAKE_BUILD_TYPE=Release
rem    cmake --build catalyst-build --config Release
rem    cmake --install catalyst-build --config Release
rem
rem 2. Eigen3  (available in vcpkg)
rem    vcpkg install eigen3
rem
rem 3. Visual Studio 2019/2022 with C++ workload (for cmake --build)
rem
rem After building, copy catalyst-paraview.dll from your ParaView installation
rem into the build\Release directory so the simulation can find it at runtime:
rem
rem    copy "C:\Program Files\ParaView 6.0.1\bin\catalyst-paraview.dll" build\Release\
rem --------------------------------------------------------------------

set CATALYST_SDK=C:\catalyst-sdk
set VCPKG_EIGEN=C:\vcpkg\packages\eigen3_x64-windows\share\eigen3
set BUILD_DIR=build

cmake -B %BUILD_DIR% -S . ^
    -DCATALYST_DIR="%CATALYST_SDK%\lib\cmake\catalyst-2.0" ^
    -DEIGEN3_DIR="%VCPKG_EIGEN%" ^
    -DCMAKE_BUILD_TYPE=Release

if errorlevel 1 (
    echo [ERROR] CMake configuration failed.
    exit /b 1
)

cmake --build %BUILD_DIR% --config Release

if errorlevel 1 (
    echo [ERROR] Build failed.
    exit /b 1
)

echo.
echo Build complete.  Executable: %BUILD_DIR%\Release\heat_equation.exe
echo.
echo Before running, copy the ParaView Catalyst implementation DLL:
echo   copy "C:\Program Files\ParaView 6.0.1\bin\catalyst-paraview.dll" %BUILD_DIR%\Release\
echo.
echo Then run:
echo   %BUILD_DIR%\Release\heat_equation.exe catalyst_pipeline.py
