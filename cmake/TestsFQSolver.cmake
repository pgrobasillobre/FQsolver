macro(add_FQSolver_runtest _name _labels)
    add_test(
        ${_name}
        python3 ${PROJECT_BINARY_DIR}/tests/${_name}/test --binary-dir=${PROJECT_BINARY_DIR}  --work-dir=${PROJECT_BINARY_DIR}/tests/${_name} --verbose)
    if(NOT "${_labels}" STREQUAL "")
        set_tests_properties(${_name} PROPERTIES LABELS "${_labels}")
    endif()
endmacro()

# Tests launched with "FQSolver"
add_FQSolver_runtest(integrate_density                                "FQSolver;Integrate Cube File")
#add_FQSolver_runtest(acceptor_donor_coulomb                           "FQSolver;Acceptor - Donor Coulomb;")
#add_FQSolver_runtest(acceptor_donor_with_overlap_integral             "FQSolver;Acceptor - Donor Coulomb + Overlap;")
#add_FQSolver_runtest(acceptor_np_charges                              "FQSolver;Acceptor - Nanoparticle (charges) Interaction;")
#add_FQSolver_runtest(acceptor_np_charges_donor_coulomb                "FQSolver;Acceptor - Nanoparticle (charges) - Donor Interaction;")
#add_FQSolver_runtest(acceptor_np_charges_dipoles_donor_coulomb        "FQSolver;Acceptor - Nanoparticle (charges + dipoles) - Donor Interaction;")
#add_FQSolver_runtest(acceptor_np_charges_donor_with_overlap_integral  "FQSolver;Acceptor - Nanoparticle (charges) - Donor Interaction + Overlap;")
#add_FQSolver_runtest(acceptor_np_charges_dipoles_donor_rotate_coulomb "FQSolver;Acceptor - Nanoparticle (charges) - Donor Interaction + Rotation;")
