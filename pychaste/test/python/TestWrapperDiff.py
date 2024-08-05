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
import unittest

from typing import List


def get_file_lines(file_path: str) -> List[str]:
    """ Load a file into a list of lines.

    :param file_path: The path to the file to load
    :return: File lines with excess whitespace and empty lines removed
    """

    with open(file_path, "r") as in_file:
        # remove excess whitespace
        lines = [line.rstrip().lstrip() for line in in_file]
        # remove empty lines
        lines = [line for line in lines if line]

    return lines


def compare_files(file_path_a: str, file_path_b: str) -> bool:
    """Check if two files have the same content.

    :param file_path_a: The path to the first file
    :param file_path_b: The path to the second file
    :return: True if the files have the same content else False
    """
    
    file_lines_a = get_file_lines(file_path_a)
    file_lines_b = get_file_lines(file_path_b)

    return file_lines_a == file_lines_b


class TestWrapperDiff(unittest.TestCase):

    def test_wrapper_diff(self) -> None:
        """Compare generated wrappers with existing wrappers."""

        old_wrappers = os.path.abspath("pychaste/dynamic/wrappers")
        new_wrappers = os.path.abspath("pychaste/dynamic/wrappers_diff")

        self.assertTrue(os.path.isdir(old_wrappers))
        self.assertTrue(os.path.isdir(new_wrappers))

        for dirpath, _, filenames in os.walk(old_wrappers):
            for filename in filenames:
                if filename.endswith(".cppwg.cpp") or filename.endswith(".cppwg.hpp"):
                    old_file = os.path.join(dirpath, filename)
                    new_file = old_file.replace(old_wrappers, new_wrappers, 1)

                    self.assertTrue(os.path.isfile(old_file))
                    self.assertTrue(os.path.isfile(new_file))

                    self.assertTrue(compare_files(old_file, new_file))


if __name__ == "__main__":
    unittest.main()
