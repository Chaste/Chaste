"""Helper classes for running visualization tests"""

__copyright__ = """Copyright (c) 2005-2024, University of Oxford.
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

from chaste.cell_based import VtkSceneModifier_2, VtkSceneModifier_3

try:
    import IPython
except ImportError as e:
    raise ImportError("Cannot use JupyterSceneModifier without Jupyter.") from e


def JupyterSceneModifierFactory(VtkSceneModifier):

    class JupyterSceneModifier(VtkSceneModifier):
        """Class for real time plotting of output"""

        def __init__(self, plotting_manager):
            self.output_format = "png"
            self.plotting_manager = plotting_manager
            super().__init__()

        def UpdateAtEndOfTimeStep(self, cell_population):
            """Update the Jupyter notebook plot with the new scene"""

            super().UpdateAtEndOfTimeStep(cell_population)

            IPython.display.clear_output(wait=True)

            if self.output_format == "png":
                IPython.display.display(
                    self.plotting_manager.vtk_show(
                        self.GetVtkScene(), output_format=self.output_format
                    )
                )
            else:
                self.plotting_manager.vtk_show(
                    self.GetVtkScene(), output_format=self.output_format
                )

    return JupyterSceneModifier


JupyterSceneModifier_2 = JupyterSceneModifierFactory(VtkSceneModifier_2)
JupyterSceneModifier_3 = JupyterSceneModifierFactory(VtkSceneModifier_3)
