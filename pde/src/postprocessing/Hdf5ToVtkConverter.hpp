/*

Copyright (c) 2005-2019, University of Oxford.
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

*/

#ifndef HDF5TOVTKCONVERTER_HPP_
#define HDF5TOVTKCONVERTER_HPP_

#include "AbstractHdf5Converter.hpp"

/**
 * This class converts from Hdf5 format to Vtk format.
 * The output will be one .vtu file with separate vtkPointData for each time step.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class Hdf5ToVtkConverter : public AbstractHdf5Converter<ELEMENT_DIM, SPACE_DIM>
{
public:
    /**
     * Constructor, which does the conversion and writes the .vtu file.
     *
     * @note This method is collective, and hence must be called by all processes.
     *
     * @param rInputDirectory The input directory, relative to CHASTE_TEST_OUTPUT, where the .h5 file has been written
     * @param rFileBaseName The base name of the data file.
     * @param pMesh Pointer to the mesh.
     * @param parallelVtk When true, write with .pvtu and fragment meshes (only works for DistributedTetrahedralMesh)
     * @param usingOriginalNodeOrdering Whether HDF5 output was written using the original node ordering
     */
    Hdf5ToVtkConverter(const FileFinder& rInputDirectory,
                       const std::string& rFileBaseName,
                       AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                       bool parallelVtk,
                       bool usingOriginalNodeOrdering);
};

#endif /*HDF5TOVTKCONVERTER_HPP_*/
