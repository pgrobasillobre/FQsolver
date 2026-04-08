def configure(options, input_files, extra_arg):
    """
    This function is used by runtest to configure runtest
    at runtime for code-specific launch command and file naming.
    """

    from os import path

    launcher = 'FQSolver'
    launcher_full_path = path.abspath(path.join(options.binary_dir, launcher))

    if not input_files:
        raise RuntimeError("No input files provided to runtest.")

    if len(input_files) > 1:
        raise RuntimeError(f"Expected 1 input file, got {len(input_files)}: {input_files}")

    inp = input_files[0]
    command = [launcher_full_path, inp]

    full_command = ' '.join(command)
    print(f"Running command: {full_command}")

    output_prefix = inp[:-4] if inp.endswith('.inp') else inp  # strip .inp or fallback
    print(f"Output prefix: {output_prefix}")

    relative_reference_path = 'reference'

    return launcher, full_command, output_prefix, relative_reference_path

