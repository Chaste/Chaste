"""Helper classes for visualising Chaste cell populations in Jupyter notebooks."""

from __future__ import annotations

__copyright__ = """Copyright (c) 2005-2026, University of Oxford.
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

try:
    import IPython
except ImportError as e:
    raise ImportError("Cannot use Jupyter imports without Jupyter.") from e

import os
import warnings

from importlib import resources
from io import StringIO
from typing import Any

import vtk
import xvfbwrapper

from chaste.cell_based import SimulationTime, VtkSceneModifier_2, VtkSceneModifier_3


class JupyterNotebookManager:
    """
    Singleton that manages cell-population rendering in a Jupyter notebook.

    On construction it starts a headless virtual display (Xvfb) so VTK can
    render off-screen. Use vtk_show() to display a VtkScene, as a static image
    or an interactive three.js plot. It is a singleton so the virtual display
    is started only once per kernel.
    """

    _instance = None

    def __new__(cls) -> JupyterNotebookManager:
        """Return the singleton instance, creating it on the first call."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self) -> None:
        """Start the virtual display and locate the JavaScript assets (once)."""
        # __init__ fires on every construction (singleton), so guard the
        # one-time setup below.
        if getattr(self, "_initialized", False):
            return
        self._initialized = True

        # Xvfb gives VTK a headless X display for its off-screen GL rendering.
        try:
            self.vdisplay = xvfbwrapper.Xvfb()
            self.vdisplay.start()
        except OSError:
            self.vdisplay = None
            warnings.warn("Could not start Xvfb.")

        self.js_resource_dir = resources.files("chaste").joinpath("external")

        # Counter for the unique DOM id given to each interactive plot's container.
        self.container_id = 0

    def _display_image(self, scene: Any, width: int, height: int) -> None:
        """Render the scene to a width x height PNG and display it."""
        # The render and PNG capture are done in C++ (GetSceneAsCharBuffer)
        # Size the render window first so the image honours the requested width/height.
        scene.GetRenderWindow().SetSize(width, height)

        # memoryview exposes the wrapped vtkUnsignedCharArray as a bytes buffer.
        data = memoryview(scene.GetSceneAsCharBuffer())
        IPython.display.display(IPython.display.Image(data))

    def _display_interactive(
        self, scene: Any, width: int, height: int, increment: bool = True
    ) -> None:
        """
        Export the scene to VRML and display it as an interactive three.js plot.
        """
        # Build the scene before rendering and exporting it (0 = time step, used
        # only for naming saved output, which this path does not write).
        scene.RenderFrame(0)

        render_window = scene.GetRenderWindow()
        render_window.SetOffScreenRendering(1)
        render_window.SetSize(width, height)
        render_window.Render()

        # Written to the working directory so the notebook can serve it back to
        # the interactive plot via its "files/" URL.
        wrl_filename = "temp_scene.wrl"
        exporter = vtk.vtkVRMLExporter()
        exporter.SetInput(render_window)
        exporter.SetFileName(os.path.join(os.getcwd(), wrl_filename))
        exporter.Write()

        if increment:
            self.container_id += 1

        # A unique container id per plot, so each output's script
        # targets its own <div> rather than another plot's.
        container = f"pychaste_plotting_container_{self.container_id}"
        
        # Each plot's output is self-contained. It loads the plotting libraries
        # itself if they are not already on the page, then renders.
        head = (
            f'<div id="{container}" style="width:{width}px; height:{height}px"></div>\n'
            '<script type="text/javascript">\n'
            "(function () {\n"
            "    function plot() {\n"
            f'        var container = document.getElementById("{container}");\n'
            f'        pychaste_plot(container, "files/{wrl_filename}", {width}, {height});\n'
            "    }\n"
            '    if (typeof pychaste_plot === "function") {\n'
            "        plot();\n"
            "    } else {\n"
        )
        tail = "\n        plot();\n    }\n})();\n</script>\n"
        # head ends at "} else {", so the libraries form the else-branch bootstrap
        # that runs only when pychaste_plot is not already defined on the page.
        html_source = head + self._read_plotting_libraries() + tail

        IPython.display.display(IPython.display.HTML(html_source))

    def _read_plotting_libraries(self) -> str:
        """
        Return the bundled three.js libraries and plotting script concatenated
        into a single block of JavaScript, ready to inject into a notebook cell.
        """
        js = StringIO()
        # Hide RequireJS's define so the three.js UMD builds set the global THREE
        # instead of registering as anonymous AMD modules; restore it afterwards.
        js.write("var __pychaste_define = window.define;\n" "window.define = undefined;\n")
        for library in ("three.min.js", "OrbitControls.js", "VRMLLoader.js"):
            with self.js_resource_dir.joinpath(library).open() as infile:
                js.write(infile.read())

        js.write("window.define = __pychaste_define;\n")

        with self.js_resource_dir.joinpath("plotting_script.js").open() as infile:
            js.write(infile.read())
        return js.getvalue()

    def vtk_show(
        self,
        scene: Any,
        width: int = 400,
        height: int = 300,
        output_format: str = "png",
        increment: bool = True,
    ) -> None:
        """
        Display a VtkScene in the current Jupyter notebook cell, sized to
        width x height pixels.

        For output_format "png" (the default) the scene is rendered to a static
        image; for "wrl" it is exported to VRML and shown as an interactive
        three.js plot (drag to rotate, scroll to zoom). Pass increment=False to
        overwrite the previous interactive plot in place instead of adding one.
        """
        if output_format == "wrl":
            self._display_interactive(scene, width, height, increment)
        else:
            self._display_image(scene, width, height)


def JupyterSceneModifierFactory(VtkSceneModifier: type) -> type:
    """
    Build a JupyterSceneModifier subclass of the given dimension-specific
    VtkSceneModifier (used to create the 2D and 3D variants below).
    """

    class JupyterSceneModifier(VtkSceneModifier):
        """
        Simulation modifier that renders the scene into the active Jupyter
        notebook cell at the end of each time step, for real-time plotting.
        """

        def __init__(
            self, plotting_manager: JupyterNotebookManager, output_format: str = "png"
        ) -> None:
            """
            Create the modifier with the manager that displays each frame and
            the output format ("png" for a static image, "wrl" for interactive).
            """
            super().__init__()
            self.plotting_manager = plotting_manager
            self.output_format = output_format

        def UpdateAtEndOfTimeStep(self, cell_population: Any) -> None:
            """
            Re-render the cell population into the notebook cell at the end of a
            time step (only on steps that are a multiple of the update frequency).
            """

            # Update the population directly rather than via the base class: its
            # render-if-due path only writes file output (unused here) and would
            # rebuild the scene redundantly before vtk_show renders it below.
            self.UpdateCellData(cell_population)

            # Only refresh the plot every GetUpdateFrequency() time steps.
            time_step = SimulationTime.Instance().GetTimeStepsElapsed()
            if time_step % self.GetUpdateFrequency() != 0:
                return

            # Clear the previous frame so this one replaces it in the same cell.
            IPython.display.clear_output(wait=True)
            self.plotting_manager.vtk_show(
                self.GetVtkScene(), output_format=self.output_format
            )

    # Give each dimension's class a distinct, informative name.
    JupyterSceneModifier.__name__ = VtkSceneModifier.__name__.replace("VtkSceneModifier", "JupyterSceneModifier")
    JupyterSceneModifier.__qualname__ = JupyterSceneModifier.__name__

    return JupyterSceneModifier


JupyterSceneModifier_2 = JupyterSceneModifierFactory(VtkSceneModifier_2)
JupyterSceneModifier_3 = JupyterSceneModifierFactory(VtkSceneModifier_3)
