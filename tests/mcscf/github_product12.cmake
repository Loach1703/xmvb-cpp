# Inject only a benchmark target; neither source tree is patched or backported.
add_executable(benchmark_product12 "${CMAKE_CURRENT_LIST_DIR}/github_product12.cpp")
target_compile_features(benchmark_product12 PRIVATE cxx_std_23)
target_link_libraries(benchmark_product12 PRIVATE xmvb_vbscf)
if(XMVB_BENCHMARK_BASELINE)
  target_compile_definitions(benchmark_product12 PRIVATE XMVB_GITHUB_BASELINE=1)
endif()
