/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
     * @param inputDirectory The input directory, relative to CHASTE_TEST_OUTPUT, where the .h5 file has been written
     * @param fileBaseName The base name of the data file.
     * @param pMesh Pointer to the mesh.
     * @param parallelVtk When true, write with .pvtu and fragment meshes (only works for DistributedTetrahedralMesh)
     * @param usingOriginalNodeOrdering Whether HDF5 output was written using the original node ordering
     */
    Hdf5ToVtkConverter(std::string inputDirectory,
                       std::string fileBaseName,
                       AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                       bool parallelVtk,
                       bool usingOriginalNodeOrdering);
};

#endif /*HDF5TOVTKCONVERTER_HPP_*/
