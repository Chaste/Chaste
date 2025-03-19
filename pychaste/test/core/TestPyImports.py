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

import chaste
import chaste.cell_based
import chaste.core
import chaste.mesh
import chaste.ode
import chaste.pde
import chaste.visualization
from chaste import *
from chaste._syntax import TemplateClassDict
from chaste.cell_based import *
from chaste.core import *
from chaste.mesh import *
from chaste.ode import *
from chaste.pde import *
from chaste.visualization import *


class TestPyImports(unittest.TestCase):

    def test_imports(self):
        module_wrapper = os.path.abspath("wrappers/all/_pychaste_all.main.cppwg.cpp")
        class_names = []

        class_regex = re.compile(r"^\s*register_(\w+)_class\(m\);\s*$")
        with open(module_wrapper, "r") as f:
            lines = f.readlines()

        for line in lines:
            class_match = class_regex.match(line)
            if class_match:
                class_name = class_match.groups(0)[0]
                if not class_name.startswith("Abstract"):
                    class_names.append(class_name)

        for class_name in class_names:
            # Check that wrapped classes are imported in PyChaste
            clas = globals().get(class_name, None)
            self.assertTrue(
                clas is not None,
                f"\nClass {class_name} is wrapped but not imported in PyChaste"
                "\n-- to fix, add an import in the relevant PyChaste module"
                " e.g. to add to pychaste core, update"
                " pychaste/src/py/chaste/core/__init__.py",
            )

            # Extra checks for templated classes
            if "_" in class_name:
                template_name, *args = class_name.split("_")
                template = globals().get(template_name, None)

                # Templated classes should have a TemplateClassDict entry
                defined = template and isinstance(template, TemplateClassDict)
                self.assertTrue(
                    defined,
                    f"\nClass {class_name} does not have a corresponding TemplateClassDict entry"
                    "\n-- to fix, add an entry in the relevant PyChaste module"
                    " e.g. to add to pychaste core, update"
                    " pychaste/src/py/chaste/core/__init__.py",
                )

                # Check that the TemplateClassDict returns the correct type
                self.assertEqual(
                    template[tuple(args)],
                    clas,
                    f"\nTemplatedClass {template_name}[{args}] does not return {class_name}"
                    "\n-- to fix, check the TemplateClassDict initialisation for {class_name}"
                    " in the relevant PyChaste module e.g. if the class is in pychaste"
                    " core, check pychaste/src/py/chaste/core/__init__.py",
                )


if __name__ == "__main__":
    unittest.main()
