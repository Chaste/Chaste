"""Visualization Module"""

__copyright__ = """Copyright (c) 2005-2025, University of Oxford.
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

import importlib.util
import warnings

from chaste._pychaste_all import (
    CellPopulationPyChasteActorGenerator_2,
    CellPopulationPyChasteActorGenerator_3,
    VtkScene_2,
    VtkScene_3,
)
from chaste._syntax import DeprecatedClass, TemplateClassDict

ipython_spec = importlib.util.find_spec("IPython")
if ipython_spec is None:
    warnings.warn("IPython not found... skipping Jupyter imports.")

else:
    from chaste.visualization._jupyter import (
        JupyterNotebookManager,
        JupyterSceneModifier_2,
        JupyterSceneModifier_3,
    )

# Template Class Syntax
CellPopulationPyChasteActorGenerator = TemplateClassDict(
    {
        ("2",): CellPopulationPyChasteActorGenerator_2,
        ("3",): CellPopulationPyChasteActorGenerator_3,
    }
)

VtkScene = TemplateClassDict(
    {
        ("2",): VtkScene_2,
        ("3",): VtkScene_3,
    }
)

if ipython_spec is not None:
    JupyterSceneModifier = TemplateClassDict(
        {
            ("2",): JupyterSceneModifier_2,
            ("3",): JupyterSceneModifier_3,
        }
    )

# Deprecated Class Syntax
CellPopulationPyChasteActorGenerator2 = DeprecatedClass("CellPopulationPyChasteActorGenerator2", CellPopulationPyChasteActorGenerator_2)
CellPopulationPyChasteActorGenerator3 = DeprecatedClass("CellPopulationPyChasteActorGenerator3", CellPopulationPyChasteActorGenerator_3)
VtkScene2 = DeprecatedClass("VtkScene2", VtkScene_2)
VtkScene3 = DeprecatedClass("VtkScene3", VtkScene_3)

if ipython_spec is not None:
    JupyterSceneModifier2 = DeprecatedClass("JupyterSceneModifier2", JupyterSceneModifier_2)
    JupyterSceneModifier3 = DeprecatedClass("JupyterSceneModifier3", JupyterSceneModifier_3)

del ipython_spec
