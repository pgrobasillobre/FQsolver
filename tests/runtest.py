#!/usr/bin/env python3

import argparse
import difflib
import os
import re
import shlex
import subprocess


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


def _normalize_text(text):
    text = text.replace("\r\n", "\n")
    text = re.sub(
        r"Normal Termination of FQSolver program in date .*",
        "Normal Termination of FQSolver program in date <normalized>",
        text,
    )
    return text


def _compare_file(actual_path, reference_path, filter_specs):
    with open(actual_path, encoding="utf-8") as handle:
        actual_text = handle.read()
    with open(reference_path, encoding="utf-8") as handle:
        reference_text = handle.read()

    actual_text = _normalize_text(_remove_filtered_blocks(actual_text, filter_specs))
    reference_text = _normalize_text(_remove_filtered_blocks(reference_text, filter_specs))

    if actual_text == reference_text:
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


def run(options, configure, input_files=None, filters=None):
    input_files = input_files or []
    filters = filters or {}

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

    if options.verbose:
        print(f"{launcher} test passed.")

    return 0
