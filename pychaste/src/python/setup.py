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

from setuptools import setup, Distribution,find_packages

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False
    
    def has_ext_modules(self):
        return True

setup(
    name = "chaste",
    version = "2024.1",
    packages = find_packages(),
    
    package_data={
        'chaste': ['_chaste_project_PyChaste_preload.so', 
                   'cell_based/_chaste_project_PyChaste_cell_based.so',
                   'core/_chaste_project_PyChaste_core.so', 
                   'mesh/_chaste_project_PyChaste_mesh.so', 
                   'ode/_chaste_project_PyChaste_ode.so', 
                   'pde/_chaste_project_PyChaste_pde.so', 
                   'visualization/_chaste_project_PyChaste_visualization.so', 
                   'tutorial/_chaste_project_PyChaste_tutorial.so',],},

    data_files = [('tutorials', ['doc/tutorials/TestMeshBasedCellSimulationsPythonTutorial.ipynb', 
                                'doc/tutorials/TestNodeBasedCellSimulationsPythonTutorial.ipynb',
                                'doc/tutorials/TestPottsBasedCellSimulationsPythonTutorial.ipynb',
                                'doc/tutorials/TestVertexBasedCellSimulationsPythonTutorial.ipynb',
                                'doc/tutorials/TestScratchAssayTutorial.ipynb',
                                'doc/tutorials/TestSpheroidTutorial.ipynb',
                                'doc/tutorials/TestTensileTestTutorial.ipynb',
                                'doc/tutorials/TestCellSortingTutorial.ipynb'])],
      
    include_package_data=True,
    zip_safe = False,

    # Project Metadata
    author = "Chaste Developers",
    author_email = "chaste-users@maillist.ox.ac.uk",
    description = "Python bindings for Chaste, a general purpose simulation package for computational biology.",
    license = "BSD-3-Clause",
    keywords = "cancer developmental biology electrophysiology scientific",

    classifiers=[
        "Development Status :: 4 - Beta",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        "Operating System :: POSIX",
        'License :: OSI Approved :: BSD License',
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: Implementation :: CPython",
    ],

    distclass=BinaryDistribution
)
