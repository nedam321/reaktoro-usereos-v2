cmake_minimum_required(VERSION 3.21)

project(ReaktoroUsereos CXX)

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

find_package(Reaktoro REQUIRED)
find_package(yaml-cpp REQUIRED)

add_library(USEREOS SHARED src/reaktoro-usereos.cpp src/reaktoro-usereos.hpp)
target_link_libraries(USEREOS PUBLIC yaml-cpp Reaktoro::Reaktoro)
target_include_directories(USEREOS PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
                                          "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")

add_custom_command(TARGET USEREOS POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:USEREOS>
                                   $<TARGET_FILE_DIR:USEREOS>
  COMMAND_EXPAND_LISTS
)

add_executable(test-usereos src/main.cpp)
target_link_libraries(test-usereos USEREOS)

add_custom_command(TARGET test-usereos POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:test-usereos>
                                   $<TARGET_FILE_DIR:test-usereos>
  COMMAND_EXPAND_LISTS
)