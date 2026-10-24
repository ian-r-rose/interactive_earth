cmake_minimum_required(VERSION 3.25)
project(interactive_earth LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(IE_PROFILE "Enable per-phase timers" ON)
option(IE_EXPECT_WASM "Fail unless configured with the Emscripten toolchain" OFF)

if(IE_EXPECT_WASM AND NOT EMSCRIPTEN)
  message(FATAL_ERROR "This preset targets WebAssembly; configure it with `emcmake cmake --preset <name>`")
endif()

if(EMSCRIPTEN)
  add_compile_options(-msimd128)
endif()

# --- FFTW -------------------------------------------------------------------
# Natively, use the installed FFTW. For WebAssembly, build a pinned release
# from source with the same toolchain.
if(EMSCRIPTEN)
  include(ExternalProject)
  ExternalProject_Add(fftw_wasm
    URL https://www.fftw.org/fftw-3.3.11.tar.gz
    URL_HASH SHA256=5630c24cdeb33b131612f7eb4b1a9934234754f9f388ff8617458d0be6f239a1
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    CMAKE_ARGS
      -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
      -DCMAKE_BUILD_TYPE=Release
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_INSTALL_LIBDIR=lib
      -DCMAKE_POLICY_VERSION_MINIMUM=3.5  # FFTW asks for CMake 3.0
      -DCMAKE_C_FLAGS=-msimd128
      -DBUILD_SHARED_LIBS=OFF
      -DBUILD_TESTS=OFF
      -DDISABLE_FORTRAN=ON
    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/libfftw3.a
  )
  ExternalProject_Get_Property(fftw_wasm INSTALL_DIR)
  file(MAKE_DIRECTORY ${INSTALL_DIR}/include)
  add_library(fftw3 STATIC IMPORTED)
  set_target_properties(fftw3 PROPERTIES
    IMPORTED_LOCATION ${INSTALL_DIR}/lib/libfftw3.a
    INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_DIR}/include)
  add_dependencies(fftw3 fftw_wasm)
else()
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
  add_library(fftw3 INTERFACE)
  target_link_libraries(fftw3 INTERFACE PkgConfig::FFTW3)
endif()

# --- Solver -----------------------------------------------------------------
add_library(ie_core STATIC src/convection.cpp src/color.cpp)
target_include_directories(ie_core PUBLIC src)
target_link_libraries(ie_core PUBLIC fftw3)
if(NOT IE_PROFILE)
  target_compile_definitions(ie_core PUBLIC IE_NO_PROFILE)
endif()
if(NOT EMSCRIPTEN)
  find_package(OpenMP REQUIRED)
  target_link_libraries(ie_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# --- Interactive app --------------------------------------------------------
# Windowing, input, and the GL context come from the vendored sokol_app.h.
add_executable(interactive_earth
  src/main.cpp src/rendering.cpp src/rendering_plugins.cpp src/sokol_app_impl.cpp)
target_include_directories(interactive_earth SYSTEM PRIVATE third_party/sokol)
target_link_libraries(interactive_earth PRIVATE ie_core)
if(EMSCRIPTEN)
  target_compile_definitions(interactive_earth PRIVATE SOKOL_GLES3)
  target_link_options(interactive_earth PRIVATE
    -sMIN_WEBGL_VERSION=2 -sMAX_WEBGL_VERSION=2
    -sINITIAL_MEMORY=64MB -sALLOW_MEMORY_GROWTH=1
    -lembind)
  configure_file(web/dev.html ${CMAKE_CURRENT_BINARY_DIR}/index.html COPYONLY)
else()
  target_compile_definitions(interactive_earth PRIVATE SOKOL_GLCORE)
  find_package(OpenGL REQUIRED)
  find_package(X11 REQUIRED)
  find_package(Threads REQUIRED)
  if(NOT X11_Xi_FOUND OR NOT X11_Xcursor_FOUND)
    message(FATAL_ERROR "sokol_app needs the Xi and Xcursor X11 extensions")
  endif()
  target_link_libraries(interactive_earth PRIVATE
    OpenGL::GL X11::X11 X11::Xi X11::Xcursor Threads::Threads ${CMAKE_DL_LIBS} m)
endif()

# --- Benchmark and tests ----------------------------------------------------
add_executable(ie_bench bench/ie_bench.cpp)
target_link_libraries(ie_bench PRIVATE ie_core)

add_executable(test_physics tests/test_physics.cpp)
target_link_libraries(test_physics PRIVATE ie_core)

if(EMSCRIPTEN)
  # Run under node with access to the real filesystem (for ie_bench --dump).
  foreach(target ie_bench test_physics)
    target_link_options(${target} PRIVATE -sALLOW_MEMORY_GROWTH=1 -sNODERAWFS=1)
  endforeach()
endif()

enable_testing()
foreach(test conduction stokes waves convection threads)
  add_test(NAME ${test} COMMAND test_physics ${test})
endforeach()
set_tests_properties(convection PROPERTIES TIMEOUT 600)
