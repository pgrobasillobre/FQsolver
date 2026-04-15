macro(add_FQSolver_runtest _name _labels)
    add_test(
        ${_name}
        python3 ${PROJECT_BINARY_DIR}/tests/${_name}/test --binary-dir=${PROJECT_BINARY_DIR}  --work-dir=${PROJECT_BINARY_DIR}/tests/${_name} --verbose)
    if(NOT "${_labels}" STREQUAL "")
        set_tests_properties(${_name} PROPERTIES LABELS "${_labels}")
    endif()
endmacro()

# Tests launched with "FQSolver"
add_FQSolver_runtest(integrate_density                                "Integrate Cube File")
add_FQSolver_runtest(potential_from_density                           "Solute - Solvent Potential")
add_FQSolver_runtest(field_from_density                               "Solute - Solvent Field")
add_FQSolver_runtest(potential_field_from_density                     "Solute - Solvent Potential and Field")
add_FQSolver_runtest(potential_field_from_density_PDB                 "Solute - Solvent Potential and Field PDB")
