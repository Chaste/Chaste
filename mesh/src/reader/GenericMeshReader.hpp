/*

Copyright (C) University of Oxford, 2005-2012

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
 * The created mesh reader is returned as a std::auto_ptr to ease memory management.
 *
 * @param rPathBaseName  the base name of the files from which to read the mesh data
 *    (either absolute, or relative to the current directory)
 * @param orderOfElements  the order of each element: 1 for linear, 2 for quadratic (defaults to 1)
 * @param orderOfBoundaryElements the order of each boundary element: 1 for linear, 2 for quadratic (defaults to 1. May
 *    or may not be different to orderOfElements (Note tetgen with the -o2 flag creates quadratic elements but doesn't
 *    create quadratic faces, hence the need for this third parameter)
 * @param readContainingElementsForBoundaryElements Whether to read in the containing element infomation
 *    for each boundary element (in the .face file if tetgen was run with '-nn').
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::auto_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > GenericMeshReader(const std::string& rPathBaseName,
                                                                             unsigned orderOfElements=1,
                                                                             unsigned orderOfBoundaryElements=1,
                                                                             bool readContainingElementsForBoundaryElements=false)
{
    std::auto_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > p_reader;
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
