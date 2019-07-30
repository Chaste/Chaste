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

#ifndef _GENERICMESHREADER_HPP_
#define _GENERICMESHREADER_HPP_

#include <memory>
#include <string>

#include "AbstractMeshReader.hpp"

// Possible mesh reader classes to create
#include "TrianglesMeshReader.hpp"
#include "MemfemMeshReader.hpp"
#include "VtkMeshReader.hpp"

/**
 * This function creates a mesh reader of a suitable type to read the mesh file given.
 * It can use any of the following readers:
 *  - TrianglesMeshReader
 *  - MemfemMeshReader
 *  - VtkMeshReader
 *
 * The created mesh reader is returned as a std::shared_ptr to ease memory management.
 *
 * @param rPathBaseName  the base name of the files from which to read the mesh data
 *    (either absolute, or relative to the current directory)
 * @param orderOfElements  the order of each element: 1 for linear, 2 for quadratic (defaults to 1)
 * @param orderOfBoundaryElements the order of each boundary element: 1 for linear, 2 for quadratic (defaults to 1. May
 *    or may not be different to orderOfElements (Note tetgen with the -o2 flag creates quadratic elements but doesn't
 *    create quadratic faces, hence the need for this third parameter)
 * @param readContainingElementsForBoundaryElements Whether to read in the containing element information
 *    for each boundary element (in the .face file if tetgen was run with '-nn').
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::shared_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > GenericMeshReader(const std::string& rPathBaseName,
                                                                             unsigned orderOfElements=1,
                                                                             unsigned orderOfBoundaryElements=1,
                                                                             bool readContainingElementsForBoundaryElements=false)
{
    std::shared_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > p_reader;
    try
    {
        p_reader.reset(new TrianglesMeshReader<ELEMENT_DIM, SPACE_DIM>(rPathBaseName,
                                                                       orderOfElements,
                                                                       orderOfBoundaryElements,
                                                                       readContainingElementsForBoundaryElements));
    }
    catch (const Exception& r_triangles_exception)
    {
        if (orderOfElements!=1 || orderOfBoundaryElements!=1 || readContainingElementsForBoundaryElements)
        {
            EXCEPTION("Quadratic meshes are only supported in Triangles format.");
        }

        try
        {
            p_reader.reset(new MemfemMeshReader<ELEMENT_DIM, SPACE_DIM>(rPathBaseName));
        }
        catch (const Exception& r_memfem_exception)
        {
#ifdef CHASTE_VTK
            try
            {
                p_reader.reset(new VtkMeshReader<ELEMENT_DIM, SPACE_DIM>(rPathBaseName));
            }
            catch (const Exception& r_vtk_exception)
            {
#endif // CHASTE_VTK
                std::string eol("\n");
                std::string combined_message = "Could not open appropriate mesh files for " + rPathBaseName + eol;
                combined_message += "Triangle format: " + r_triangles_exception.GetShortMessage() + eol;
                combined_message += "Memfem format: " + r_memfem_exception.GetShortMessage() + eol;
#ifdef CHASTE_VTK
                combined_message += "Vtk format: " + r_vtk_exception.GetShortMessage() + eol;
#endif // CHASTE_VTK
                EXCEPTION(combined_message);
#ifdef CHASTE_VTK
            }
#endif // CHASTE_VTK
        }
    }
    return p_reader;
}

#endif //_GENERICMESHREADER_HPP_
