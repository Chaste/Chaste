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

#ifndef HDF5TOCMGUICONVERTER_HPP_
#define HDF5TOCMGUICONVERTER_HPP_

#include "AbstractHdf5Converter.hpp"

/**
 *  This class converts from Hdf5 format to Cmgui format.
 *  The output will be one .exnode file per time step.
 *  The format that cmgui accepts is (after the headers):
 *
 *   Node: 1
 *   Value_at_node_1
 *   Node:2
 *   Value_at_node_2
 *   .....
 *
 *   For bidomain simulations, we will have two fields, one for Vm and one for Phie.
 *   The Cmgui format for two fields is as follows:
 *
 *   Node: 1
 *   Vm_node_1
 *   Phie_at_node_1
 *   Node:2
 *   Vm_at_node_2
 *   Phie_at_node_2
 *   .....
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class Hdf5ToCmguiConverter : public AbstractHdf5Converter<ELEMENT_DIM, SPACE_DIM>
{
private:

    /**
     * A helper method which takes in a string, which must be 'Mono' or 'Bi'
     * and reads the data from the hdf5 file, writing it out in Cmgui format.
     *
     * @param type the type of simulation (Mono  or Bi)
     */
    void Write(std::string type);

    /**
     * Writes a basic script for visualization of mesh and data.
     *
     * It loads the nodes and elements. It also asks cmgui to create faces and
     * lines (because our output files don't have the information). Data are
     * loaded by means of a 'for' loop.
     *
     * After loading the script, nodes (seen as small dots) and lines connecting
     * them will be displayed. Both nodes  and lines will be coloured according
     * to the first variable in the Hdf5 file (usually Vm). The Cmgui default
     * spectrum (0 to 1, blue to red) is used. This can be changed by clicking
     * on 'graphics->spectrum editor'. Other manual modification can be made by
     * clicking on 'graphics->scene editor'.
     */
    void WriteCmguiScript();

public:

    /**
     * Constructor, which does the conversion.
     *
     * @note This method is collective, and hence must be called by all processes.
     *
     * @param inputDirectory The input directory, relative to CHASTE_TEST_OUTPUT, where the .h5 file has been written
     * @param fileBaseName The base name of the data file.
     * @param pMesh Pointer to the mesh.
     * @param hasBath whether the mesh has a bath or not. Defaults to false.
     */
    Hdf5ToCmguiConverter(std::string inputDirectory,
                         std::string fileBaseName,
                         AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                         bool hasBath = false);
};

#endif /*HDF5TOCMGUICONVERTER_HPP_*/
