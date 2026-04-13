#!/usr/bin/env python3

import argparse
import difflib
import math
import os
import re
import shlex
import subprocess

NUMBER_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
DEFAULT_REL_TOLERANCE = 1.0e-12
DEFAULT_ABS_TOLERANCE = 1.0e-12


def get_filter(from_string=None, to_string=None, rel_tolerance=None):
    return {
        "from_string": from_string,
        "to_string": to_string,
        "rel_tolerance": rel_tolerance,
    }


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary-dir", required=True)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def _remove_filtered_blocks(text, filter_specs):
    lines = text.splitlines()
    keep = [True] * len(lines)

    for spec in filter_specs:
        start = spec.get("from_string")
        end = spec.get("to_string")
        if not start or not end:
            continue

        in_block = False
        for idx, line in enumerate(lines):
            if start in line:
                in_block = True
            if in_block:
                keep[idx] = False
            if in_block and end in line:
                in_block = False

    filtered = [line for idx, line in enumerate(lines) if keep[idx]]
    return "\n".join(filtered) + ("\n" if text.endswith("\n") else "")


def _normalize_text(text, replacements=None):
    text = text.replace("\r\n", "\n")
    text = re.sub(
        r"Normal Termination of FQSolver program in date .*",
        "Normal Termination of FQSolver program in date <normalized>",
        text,
    )
    text = re.sub(r"OMP Threads:\s+\d+", "OMP Threads: <normalized>", text)
    text = re.sub(r"/[^\n]*?/tests/[^/\n]*_from_density", "<WORK_DIR>", text)
    for old, new in replacements or []:
        text = text.replace(old, new)
    return text


def _numeric_texts_close(reference_text, actual_text, rel_tolerance, abs_tolerance):
    reference_numbers = NUMBER_RE.findall(reference_text)
    actual_numbers = NUMBER_RE.findall(actual_text)

    if len(reference_numbers) != len(actual_numbers):
        return False

    reference_skeleton = NUMBER_RE.sub("<NUM>", reference_text)
    actual_skeleton = NUMBER_RE.sub("<NUM>", actual_text)
    if reference_skeleton != actual_skeleton:
        return False

    for reference_number, actual_number in zip(reference_numbers, actual_numbers):
        if not math.isclose(
            float(reference_number),
            float(actual_number),
            rel_tol=rel_tolerance,
            abs_tol=abs_tolerance,
        ):
            return False

    return True


def _compare_file(actual_path, reference_path, filter_specs, replacements=None):
    with open(actual_path, encoding="utf-8") as handle:
        actual_text = handle.read()
    with open(reference_path, encoding="utf-8") as handle:
        reference_text = handle.read()

    actual_text = _normalize_text(_remove_filtered_blocks(actual_text, filter_specs), replacements)
    reference_text = _normalize_text(_remove_filtered_blocks(reference_text, filter_specs), replacements)

    if actual_text == reference_text:
        return True, ""

    if _numeric_texts_close(
        reference_text,
        actual_text,
        DEFAULT_REL_TOLERANCE,
        DEFAULT_ABS_TOLERANCE,
    ):
        return True, ""

    diff = "".join(
        difflib.unified_diff(
            reference_text.splitlines(keepends=True),
            actual_text.splitlines(keepends=True),
            fromfile=reference_path,
            tofile=actual_path,
        )
    )
    return False, diff


def run(options, configure, input_files=None, filters=None, output_files=None):
    input_files = input_files or []
    filters = filters or {}
    output_files = output_files or []

    launcher, full_command, output_prefix, relative_reference_path = configure(
        options, input_files, None
    )

    command = shlex.split(full_command)
    completed = subprocess.run(command, cwd=options.work_dir, check=False)
    if completed.returncode != 0:
        return completed.returncode

    for extension, filter_specs in filters.items():
        actual_path = os.path.join(options.work_dir, f"{output_prefix}.{extension}")
        reference_path = os.path.join(
            options.work_dir, relative_reference_path, f"{output_prefix}.{extension}"
        )

        if not os.path.exists(actual_path):
            print(f"Missing actual output file: {actual_path}")
            return 1
        if not os.path.exists(reference_path):
            print(f"Missing reference output file: {reference_path}")
            return 1

        ok, diff = _compare_file(actual_path, reference_path, filter_specs)
        if not ok:
            print(f"Output mismatch for {extension} file.")
            print(diff)
            return 1

    replacements = [(os.path.abspath(options.work_dir), "<WORK_DIR>")]
    for actual_relpath, reference_relpath in output_files:
        actual_path = os.path.join(options.work_dir, actual_relpath)
        reference_path = os.path.join(options.work_dir, reference_relpath)

        if not os.path.exists(actual_path):
            print(f"Missing actual output file: {actual_path}")
            return 1
        if not os.path.exists(reference_path):
            print(f"Missing reference output file: {reference_path}")
            return 1

        ok, diff = _compare_file(actual_path, reference_path, [], replacements)
        if not ok:
            print(f"Output mismatch for {actual_relpath} file.")
            print(diff)
            return 1

    if options.verbose:
        print(f"{launcher} test passed.")

    return 0
