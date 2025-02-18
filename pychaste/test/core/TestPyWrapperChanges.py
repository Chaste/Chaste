"""Copyright (c) 2005-2025, University of Oxford.
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


class TestPyWrapperChanges(unittest.TestCase):
    """Test that the generated wrappers are up to date."""

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
                    and not re.search(r"Unknown class pack<.*>", line)
                    and not re.search(r"Unknown class [\w]*Iterator\b", line)
                ):
                    unknown_classes.append(line)
        self.assertEqual(
            len(unknown_classes),
            0,
            "\n" + "".join(unknown_classes) + "Found unknown classes"
            "\n-- to wrap, add relevant entries to config.yaml. "
            "\n-- to exclude from wrapping, add to config.yaml with the `exclude` option.",
        )


if __name__ == "__main__":
    unittest.main()
