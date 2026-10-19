file(MAKE_DIRECTORY "${RUN_DIR}")
configure_file("${INPUT}" "${RUN_DIR}/input.xmi" COPYONLY)
configure_file("${BASIS}" "${RUN_DIR}/cc-pvdz.gbs" COPYONLY)
execute_process(COMMAND "${PROGRAM}" "${RUN_DIR}/input.xmi"
  WORKING_DIRECTORY "${RUN_DIR}"
  RESULT_VARIABLE status OUTPUT_VARIABLE report ERROR_VARIABLE errors)
if (NOT status EQUAL 0)
  message(FATAL_ERROR "RKS .xmi execution failed: ${errors}")
endif()
foreach(pattern IN ITEMS "Restricted Kohn-Sham DFT" "Converged: yes"
    "Libxc: [0-9]" "SCF iterations: [1-9][0-9]*" "Final energy: -[0-9]")
  if (NOT report MATCHES "${pattern}")
    message(FATAL_ERROR "RKS report is missing ${pattern}")
  endif()
endforeach()
if (NOT EXISTS "${RUN_DIR}/input.orb")
  message(FATAL_ERROR "RKS did not write converged orbitals")
endif()
file(READ "${RUN_DIR}/input.orb" orbitals)
if (NOT orbitals MATCHES "# ORBITAL 1 NAO = [1-9]" OR
    orbitals MATCHES "\\$GUS|\\$END|[0-9][eE][-+][0-9]")
  message(FATAL_ERROR "RKS orbital output does not follow fixed-point GUS-body format")
endif()
file(READ "${RUN_DIR}/input.xmi" input)
string(REPLACE "LDA_X,LDA_C_VWN" "GGA_X_PBE" unsupported "${input}")
file(WRITE "${RUN_DIR}/unsupported.xmi" "${unsupported}")
execute_process(COMMAND "${PROGRAM}" "${RUN_DIR}/unsupported.xmi"
  WORKING_DIRECTORY "${RUN_DIR}"
  RESULT_VARIABLE unsupported_status OUTPUT_QUIET ERROR_VARIABLE unsupported_error)
if (unsupported_status EQUAL 0 OR NOT unsupported_error MATCHES "LDA")
  message(FATAL_ERROR "unsupported GGA was not rejected explicitly")
endif()
