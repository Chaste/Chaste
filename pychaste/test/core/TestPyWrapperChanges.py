"""Copyright (c) 2005-2024, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import os
import re
import unittest
from difflib import context_diff


def file_diff(file_a: str, file_b: str) -> bool:
    """Check if two files have the same content.

    :param file_a: The path to the first file
    :param file_b: The path to the second file
    :return: A diff of the two files
    """
    # Read files and remove excess whitespace
    with open(file_a, "r") as fa:
        a = [line.strip() for line in fa]
        a = [line for line in a if line]

    with open(file_b, "r") as fb:
        b = [line.strip() for line in fb]
        b = [line for line in b if line]

    return "\n".join(context_diff(a, b))


class TestPyWrapperChanges(unittest.TestCase):
    """Test that the generated wrappers are up to date."""

    def test_wrapper_changes(self) -> None:
        """Check for wrapper changes."""
        src_wrapper_dir = os.path.abspath("wrappers")
        gen_wrapper_dir = os.path.abspath("wrappers.gen")

        self.assertTrue(
            os.path.isdir(src_wrapper_dir),
            "Cannot find original wrappers: " + src_wrapper_dir,
        )

        self.assertTrue(
            os.path.isdir(gen_wrapper_dir),
            "Cannot find generated wrappers: " + gen_wrapper_dir,
        )

        for dirpath, _, filenames in os.walk(src_wrapper_dir):
            for filename in filenames:
                if filename.endswith(".cppwg.cpp") or filename.endswith(".cppwg.hpp"):
                    src_file = os.path.join(dirpath, filename)
                    gen_file = os.path.join(
                        gen_wrapper_dir,
                        os.path.relpath(src_file, src_wrapper_dir),
                    )

                    self.assertTrue(
                        os.path.isfile(gen_file),
                        f"\n[{os.path.relpath(src_file, src_wrapper_dir)}]"
                        "\n-- Found ungenerated wrapper in src"
                        "\n-- if it has been deliberately removed from"
                        " config.yaml, the wrapper file should be removed as well"
                        "\n-- if it is not to be removed, the relevant class"
                        " entry should be added in config.yaml.",
                    )

                    diff = file_diff(src_file, gen_file)
                    self.assertEqual(
                        diff,
                        "",
                        f"\n[{os.path.relpath(src_file, src_wrapper_dir)}]"
                        "\n-- Generated wrapper differs from existing"
                        "\n-- to update, run `make pychaste_wrappers`.",
                    )

        for dirpath, _, filenames in os.walk(gen_wrapper_dir):
            for filename in filenames:
                if filename.endswith(".cppwg.cpp") or filename.endswith(".cppwg.hpp"):
                    gen_file = os.path.join(dirpath, filename)
                    src_file = os.path.join(
                        src_wrapper_dir,
                        os.path.relpath(gen_file, gen_wrapper_dir),
                    )
                    self.assertTrue(
                        os.path.isfile(src_file),
                        f"\n[{os.path.relpath(gen_file, gen_wrapper_dir)}]"
                        "\n-- New wrapper generated"
                        "\n-- to add to src, run `make pychaste_wrappers`"
                        "\n-- to exclude, add `exclude` option in config.yaml.",
                    )

    def test_unwrapped_classes(self) -> None:
        """Find unwrapped classes."""

        log_file = os.path.abspath("cppwg.log")
        self.assertTrue(
            os.path.isfile(log_file),
            "Cannot find wrapper generator logs: " + log_file,
        )

        unknown_classes = []
        with open(log_file, "r") as lf:
            for line in lf:
                if (
                    "Unknown class" in line
                    and "/cell_based/src/" in line
                    and not "/cell_based/src/fortests/" in line
                    and not re.search(r"Unknown class guid_defined<.*>", line)
                    and not re.search(r"Unknown class [\w]*Iterator\b", line)
                ):
                    unknown_classes.append(line)
        self.assertEqual(
            len(unknown_classes),
            0,
            "\n" + "".join(unknown_classes) + "Found unknown classes"
            "\n-- to wrap, add relevant entry to config.yaml. "
            "\n-- to exclude from wrapping, add to config.yaml with `exclude` option.",
        )


if __name__ == "__main__":
    unittest.main()
