file(MAKE_DIRECTORY "${RUN_DIR}")
configure_file("${INPUT}" "${RUN_DIR}/input.xmi" COPYONLY)
configure_file("${BASIS}" "${RUN_DIR}/cc-pvdz.gbs" COPYONLY)
execute_process(COMMAND "${PROGRAM}" "${RUN_DIR}/input.xmi"
  WORKING_DIRECTORY "${RUN_DIR}"
  RESULT_VARIABLE status OUTPUT_VARIABLE report ERROR_VARIABLE errors)
if (NOT status EQUAL 0)
  message(FATAL_ERROR "RKS .xmi execution failed: ${errors}")
endif()
configure_file("${UKS_INPUT}" "${RUN_DIR}/uks.xmi" COPYONLY)
execute_process(COMMAND "${PROGRAM}" "${RUN_DIR}/uks.xmi"
  WORKING_DIRECTORY "${RUN_DIR}"
  RESULT_VARIABLE uks_status OUTPUT_VARIABLE uks_report ERROR_VARIABLE uks_errors)
if (NOT uks_status EQUAL 0)
  message(FATAL_ERROR "UKS .xmi execution failed: ${uks_errors}")
endif()
foreach(pattern IN ITEMS "Unrestricted Kohn-Sham DFT" "Converged: yes"
    "Alpha electrons: 1" "Beta electrons: 0" "Spin square: 0.75"
    "XC: GGA_X_PBE,GGA_C_PBE")
  if (NOT uks_report MATCHES "${pattern}")
    message(FATAL_ERROR "UKS report is missing ${pattern}")
  endif()
endforeach()
foreach(spin IN ITEMS alpha beta)
  if (NOT EXISTS "${RUN_DIR}/uks.${spin}.orb")
    message(FATAL_ERROR "UKS did not write ${spin} orbitals")
  endif()
endforeach()
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
string(REPLACE "LDA_X,LDA_C_VWN" "HYB_GGA_XC_PBEH" unsupported "${input}")
file(WRITE "${RUN_DIR}/unsupported.xmi" "${unsupported}")
execute_process(COMMAND "${PROGRAM}" "${RUN_DIR}/unsupported.xmi"
  WORKING_DIRECTORY "${RUN_DIR}"
  RESULT_VARIABLE unsupported_status OUTPUT_QUIET ERROR_VARIABLE unsupported_error)
if (unsupported_status EQUAL 0 OR NOT unsupported_error MATCHES "LDA|GGA|hybrid|Hybrid")
  message(FATAL_ERROR "unsupported hybrid was not rejected explicitly")
endif()
