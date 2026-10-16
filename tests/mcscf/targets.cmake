# Common solver tests must link without the VBSCF or CASSCF model libraries.
set(_xmvb_mcscf_tests
  block_inverse_bfgs ciah neo neo_globalization
  positive_ritz_secants shifted_metric_preconditioner state_weights)
foreach(name IN LISTS _xmvb_mcscf_tests)
  add_executable(test_${name} ${CMAKE_CURRENT_LIST_DIR}/unit/test_${name}.cpp)
  target_link_libraries(test_${name} PRIVATE xmvb_mcscf)
  if (BUILD_TESTING)
    add_test(NAME ${name} COMMAND test_${name})
    set_tests_properties(${name} PROPERTIES
      ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  endif()
endforeach()
unset(_xmvb_mcscf_tests)
