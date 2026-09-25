cmake_minimum_required(VERSION 3.16)
project(OversampledTH LANGUAGES CXX)

find_package(ROOT CONFIG REQUIRED COMPONENTS Hist ROOTDataFrame)

add_library(OversampledTH INTERFACE)
target_include_directories(OversampledTH INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(OversampledTH INTERFACE cxx_std_17)
target_link_libraries(OversampledTH INTERFACE ROOT::Hist ROOT::ROOTDataFrame)

option(OVERSAMPLEDTH_BUILD_EXAMPLES "Build C++ examples" ON)
if(OVERSAMPLEDTH_BUILD_EXAMPLES)
  foreach(example IN ITEMS sequential slot_local range_aware)
    add_executable(example_${example} examples/cpp/${example}.cpp)
    target_link_libraries(example_${example} PRIVATE OversampledTH)
  endforeach()
endif()

option(OVERSAMPLEDTH_BUILD_STUDIES "Build benchmark executables" ON)
if(OVERSAMPLEDTH_BUILD_STUDIES)
  add_executable(benchmark_actions studies/benchmark.cpp)
  target_link_libraries(benchmark_actions PRIVATE OversampledTH)
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(oversampled_histogram_tests test/test_histograms.cpp)
  target_link_libraries(oversampled_histogram_tests PRIVATE OversampledTH)
  add_test(NAME oversampled_histogram_tests COMMAND oversampled_histogram_tests)

  add_executable(coherence_tests test/test_coherence.cpp)
  target_link_libraries(coherence_tests PRIVATE OversampledTH)
  add_test(NAME coherence_tests COMMAND coherence_tests)

  add_executable(stress_tests test/test_stress.cpp)
  target_link_libraries(stress_tests PRIVATE OversampledTH)
  add_test(NAME stress_tests COMMAND stress_tests)
endif()
