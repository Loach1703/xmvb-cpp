file(MAKE_DIRECTORY "${RUN_DIR}")
configure_file("${INPUT}" "${RUN_DIR}/input.xmi" COPYONLY)
execute_process(COMMAND "${PROGRAM}" input.xmi
  WORKING_DIRECTORY "${RUN_DIR}"
  RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
file(WRITE "${RUN_DIR}/output.out" "${output}")
file(WRITE "${RUN_DIR}/run.err" "${error}")
if(NOT status EQUAL 0 OR NOT output MATCHES "Converged: yes")
  message(FATAL_ERROR "CASSCF .xmi smoke failed: ${status}\n${output}\n${error}")
endif()
if(NOT output MATCHES "CASSCF \\(exact, single-state, CIAH\\)" OR
    NOT EXISTS "${RUN_DIR}/input.orb")
  message(FATAL_ERROR "CASSCF method dispatch or converged orbital output is missing")
endif()
file(READ "${RUN_DIR}/input.orb" orbitals)
if(orbitals MATCHES "\\$GUS|\\$END|[0-9][eE][+-][0-9]")
  message(FATAL_ERROR "CASSCF .orb is not the fixed-point GUS body")
endif()
