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

import os.path

from io import StringIO
from pkg_resources import resource_filename

import vtk
import xvfbwrapper

try:
    import IPython
except ImportError as e:
    raise ImportError("Cannot use JupyterNotebookManager without Jupyter.") from e


class JupyterNotebookManager:
    """Singleton class for managing plotting in a Jupyter notebook"""

    def __new__(cls, *args, **kwds):
        """Singleton pattern"""
        # https://www.python.org/download/releases/2.2/descrintro/#__new__
        it = cls.__dict__.get("__it__", None)
        if it is not None:
            return it
        it = object.__new__(cls)
        cls.__it__ = it
        it.init(*args, **kwds)
        return it

    def init(self, *args, **kwds):
        self.interactive_plotting_loaded = False
        self.three_js_dir = resource_filename("chaste", os.path.join("external"))
        self.container_id = 0

        try:
            self.vdisplay = xvfbwrapper.Xvfb()
            self.vdisplay.start()
        except OSError:
            self.vdisplay = None

        self.renderWindow = vtk.vtkRenderWindow()

    def _interactive_plot_init(self):
        if not self.interactive_plotting_loaded:

            library_javascript = StringIO()
            library_javascript.write(
                """
            <script type="text/javascript">
            var pychaste_javascript_injected = true;
            """
            )

            three_js_libraries = (
                "three.min.js",
                "OrbitControls.js",
                "VRMLLoader.js",
                "Detector.js",
            )

            # Include three.js
            for library in three_js_libraries:
                with open(os.path.join(self.three_js_dir, library)) as infile:
                    library_javascript.write(infile.read())

            # Include internal plotting functions
            with open(os.path.join(self.three_js_dir, "plotting_script.js")) as infile:
                library_javascript.write(infile.read())

            self.interactive_plotting_loaded = True
            IPython.display(IPython.HTML(library_javascript.getvalue()))

    def _interactive_plot_show(
        self, width, height, file_name="temp_scene.wrl", increment=True
    ):

        self._interactive_plot_init()
        if increment:
            self.container_id += 1

        html_source = f"""
        <div id="pychaste_plotting_container_{self.container_id}" style="width:{width}px; height:{height}px">
            &nbsp;
        </div>
        <script type="text/javascript">
            (function(){{
            var three_container = document.getElementById("pychaste_plotting_container_{self.container_id}");
            pychaste_plot(three_container, files/{file_name}, {width}, {height});
            }})();
        </script>
        """

        IPython.display(IPython.HTML(html_source))

    def vtk_show(
        self, scene, width=400, height=300, output_format="png", increment=True
    ):
        """
        Takes vtkRenderer instance and returns an IPython Image with the rendering.
        """

        scene.ResetRenderer(0)
        renderer = scene.GetRenderer()

        self.renderWindow.SetOffScreenRendering(1)
        self.renderWindow.AddRenderer(renderer)
        self.renderWindow.SetSize(width, height)
        self.renderWindow.Render()

        if output_format == "wrl":
            exporter = vtk.vtkVRMLExporter()
            exporter.SetInput(self.renderWindow)
            exporter.SetFileName(os.getcwd() + "/temp_scene.wrl")
            exporter.Write()
            self._interactive_plot_show(width, height, "temp_scene.wrl", increment)
            self.renderWindow.RemoveRenderer(renderer)

        else:
            windowToImageFilter = vtk.vtkWindowToImageFilter()
            windowToImageFilter.SetInput(self.renderWindow)
            windowToImageFilter.Update()

            writer = vtk.vtkPNGWriter()
            writer.SetWriteToMemory(1)
            writer.SetInputConnection(windowToImageFilter.GetOutputPort())
            writer.Write()
            data = memoryview(writer.GetResult())
            self.renderWindow.RemoveRenderer(renderer)

            return IPython.Image(data)
