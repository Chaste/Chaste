#!/usr/bin/env python
"""Helper classes for running tests
"""

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

PYCHASTE_CAN_IMPORT_IPYTHON = True
PYCHASTE_CAN_IMPORT_VTK = True

from six.moves import StringIO
import os.path
from pkg_resources import resource_filename
from xvfbwrapper import Xvfb

import chaste.cell_based

try:
    import IPython                    
except ImportError:
    PYCHASTE_CAN_IMPORT_IPYTHON = False  
    
try:
    import vtk                    
except ImportError:
    PYCHASTE_CAN_IMPORT_VTK = False      

if PYCHASTE_CAN_IMPORT_IPYTHON:
    #from IPython import display
    from IPython.display import Image, HTML, display
    
    
    def vtk_show_basic(scene, width=800, height=600):
        
        """
        Takes a Scene instance and returns an IPython Image with the rendering.
        """
        
        data = memoryview(scene.GetSceneAsCharBuffer())
        
        return Image(data)
    
    if PYCHASTE_CAN_IMPORT_VTK:
        
        class JupyterNotebookManager(object):
            
            interactive_plotting_loaded = False
            three_js_dir = resource_filename('chaste', os.path.join('external'))
            container_id = 0

            def __new__(cls, *args, **kwds):
                it = cls.__dict__.get("__it__")
                if it is not None:
                    return it
                cls.__it__ = it = object.__new__(cls)
                it.init(*args, **kwds)
                return it

            def init(self, *args, **kwds):
                try:
                    self.vdisplay = Xvfb()
                    self.vdisplay.start()
                except OSError:
                    self.vdisplay = None

                self.renderWindow = vtk.vtkRenderWindow()
            
            def interactive_plot_init(self):
                
                if not JupyterNotebookManager.interactive_plotting_loaded:
                
                    library_javascript = StringIO()
                    library_javascript.write("""
                    <script type="text/javascript">
                    var pychaste_javascript_injected = true;
                    """)
                    
                    three_js_libraries = ("three.min.js", "OrbitControls.js", "VRMLLoader.js", "Detector.js")
                    
                    # Include three.js
                    for library in three_js_libraries:
                        with open(os.path.join(JupyterNotebookManager.three_js_dir, library)) as infile:
                            library_javascript.write(infile.read())
                            
                    # Include internal plotting functions
                    with open(os.path.join(JupyterNotebookManager.three_js_dir, "plotting_script.js")) as infile:
                        library_javascript.write(infile.read())
                        
                    JupyterNotebookManager.interactive_plotting_loaded = True
                    display(HTML(library_javascript.getvalue()))
                    
            def interactive_plot_show(self, width, height, file_name = "temp_scene.wrl", increment = True):
                
                self.interactive_plot_init()
                if increment:
                    JupyterNotebookManager.container_id +=1 
                
                html_source = """
                <div id="pychaste_plotting_container_{container_id}" style="width:{width}px; height:{height}px">
                <script type="text/javascript">
                    (function(){{
                    var three_container = document.getElementById("pychaste_plotting_container_{container_id}");
                    pychaste_plot(three_container, {file}, {width}, {height});
                    }})();
                </script>
                """.format(width=str(width), height=str(height), 
                           container_id=str(JupyterNotebookManager.container_id), file='"files/' + file_name + '"')

                display(HTML(html_source))
        
            def vtk_show(self, scene, width=400, height=300, output_format = "png", increment = True):
                
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
                    self.interactive_plot_show(width, height, "temp_scene.wrl", increment)
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
                     
                    return Image(data)
                
    class JupyterSceneModifier2(chaste.cell_based.VtkSceneModifier2):
        
        """ Class for real time plotting of output
        """

        def __init__(self, plotting_manager):
            
            self.output_format = 'png'
            self.plotting_manager = plotting_manager
            super(JupyterSceneModifier2, self).__init__()
            
            
        def UpdateAtEndOfTimeStep(self, cell_population):
            
            """ Update the Jupyter notebook plot with the new scene
            
            """
            
            super(JupyterSceneModifier2, self).UpdateAtEndOfTimeStep(cell_population)
            
            IPython.display.clear_output(wait=True)
            
            if self.output_format == 'png':
                IPython.display.display(self.plotting_manager.vtk_show(self.GetVtkScene(), output_format=self.output_format))
            else:
                self.plotting_manager.vtk_show(self.GetVtkScene(), output_format=self.output_format)
            
            
    class JupyterSceneModifier3(chaste.cell_based.VtkSceneModifier3):
        
        """ Class for real time plotting of output
        """

        def __init__(self, plotting_manager):
            
            self.output_format = 'png'
            self.plotting_manager = plotting_manager
            super(JupyterSceneModifier3, self).__init__()
            
            
        def UpdateAtEndOfTimeStep(self, cell_population):
            
            """ Update the Jupyter notebook plot with the new scene
            
            """
            
            super(JupyterSceneModifier3, self).UpdateAtEndOfTimeStep(cell_population)
            
            IPython.display.clear_output(wait=True)
            
            if self.output_format == 'png':
                IPython.display.display(self.plotting_manager.vtk_show(self.GetVtkScene(), output_format=self.output_format))
            else:
                self.plotting_manager.vtk_show(self.GetVtkScene(), output_format=self.output_format)

