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
     * @param rInputDirectory The input directory, relative to CHASTE_TEST_OUTPUT, where the .h5 file has been written
     * @param rFileBaseName The base name of the data file.
     * @param pMesh Pointer to the mesh.
     * @param hasBath whether the mesh has a bath or not. Defaults to false.
     * @param precision  The number of digits to use when printing to file (0u = default).
     */
    Hdf5ToCmguiConverter(const FileFinder& rInputDirectory,
                         const std::string& rFileBaseName,
                         AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                         bool hasBath = false,
                         unsigned precision = 0u);
};

#endif /*HDF5TOCMGUICONVERTER_HPP_*/
