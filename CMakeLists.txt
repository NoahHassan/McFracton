cmake_minimum_required(VERSION 3.21)
project(McFracton LANGUAGES CXX)

# Visual Studio (McFracton.sln) stays the primary build on Windows.
# This file exists so the core and the headless runner also build on the Linux cluster.
# Sources are globbed on purpose: adding a file in Visual Studio needs no edit here.

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(MCF_BUILD_GUI "Build the SFML/ImGui viewer (needs sfml + imgui + imgui-sfml)" ${WIN32})

set(MCF_DIR ${CMAKE_CURRENT_SOURCE_DIR}/McFracton)

# ---------------------------------------------------------------- core library
file(GLOB_RECURSE MCF_CORE_SOURCES CONFIGURE_DEPENDS ${MCF_DIR}/core/src/*.cpp)
file(GLOB_RECURSE MCF_CORE_HEADERS CONFIGURE_DEPENDS ${MCF_DIR}/core/include/*.h)

add_library(mcf_core STATIC ${MCF_CORE_SOURCES} ${MCF_CORE_HEADERS})
target_include_directories(mcf_core PUBLIC
    ${MCF_DIR}/core/include
    ${MCF_DIR}/core/include/systems)

if(MSVC)
    target_compile_options(mcf_core PUBLIC /W3)
else()
    target_compile_options(mcf_core PUBLIC -Wall -Wextra -Wno-unused-parameter)
    target_compile_options(mcf_core PUBLIC $<$<CONFIG:Release>:-O3>)
endif()

# ------------------------------------------------------- headless runner (CLI)
# Built once cli/main.cpp exists (Phase 0 step 7 / Phase 4).
file(GLOB MCF_CLI_SOURCES CONFIGURE_DEPENDS ${MCF_DIR}/cli/*.cpp)
if(MCF_CLI_SOURCES)
    add_executable(mcf_run ${MCF_CLI_SOURCES})
    target_link_libraries(mcf_run PRIVATE mcf_core)
endif()

# ------------------------------------------------------------------- GUI (opt)
if(MCF_BUILD_GUI)
    file(GLOB MCF_GUI_SOURCES CONFIGURE_DEPENDS ${MCF_DIR}/gui/*.cpp)
    file(GLOB MCF_GUI_HEADERS CONFIGURE_DEPENDS ${MCF_DIR}/gui/*.h)

    find_package(SFML 2 COMPONENTS graphics window system CONFIG REQUIRED)
    find_package(imgui CONFIG REQUIRED)
    find_package(ImGui-SFML CONFIG REQUIRED)

    add_executable(mcf_gui ${MCF_GUI_SOURCES} ${MCF_GUI_HEADERS})
    target_include_directories(mcf_gui PRIVATE ${MCF_DIR}/gui)
    target_link_libraries(mcf_gui PRIVATE mcf_core sfml-graphics sfml-window sfml-system ImGui-SFML::ImGui-SFML)

    # The app loads its font from the working directory, as the VS build does.
    add_custom_command(TARGET mcf_gui POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${MCF_DIR}/DMSans-VariableFont_opsz,wght.ttf" $<TARGET_FILE_DIR:mcf_gui>)
endif()
