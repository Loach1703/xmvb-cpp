# Linear algebra is independently buildable without any electronic model.
foreach(name davidson eigen_response minres orthonormal_hvp_basis dense_span)
  add_executable(test_${name} ${CMAKE_CURRENT_LIST_DIR}/unit/test_${name}.cpp)
  target_link_libraries(test_${name} PRIVATE xmvb_linalg)
  if (BUILD_TESTING)
    add_test(NAME ${name} COMMAND test_${name})
    set_tests_properties(${name} PROPERTIES
      ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
  endif()
endforeach()
