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

#ifndef POTTSMESHGENERATOR_HPP_
#define POTTSMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "PottsMesh.hpp"

/**
 * Generator of regular Potts meshes, used as starting points for many simulations.
 *
 * This class takes in options such as width, height,
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
template<unsigned DIM>
class PottsMeshGenerator
{
protected:

    /** A pointer to the mesh this class creates */
    PottsMesh<DIM>* mpMesh;

public:

    /**
     * Constructor.
     *
     * @param numNodesAcross  The number of columns of nodes in the mesh
     * @param numElementsAcross  The number of columns of elements in the mesh
     * @param elementWidth  Width of the elements
     * @param numNodesUp  The number of rows of nodes in the mesh (defaults to 1)
     * @param numElementsUp  The number of rows of elements in the mesh (defaults to 1)
     * @param elementHeight  Height of the elements (defaults to 1)
     * @param numNodesDeep  The number nodes deep for this mesh (defaults to 1)
     * @param numElementsDeep  The number of elements deep for this mesh (defaults to 1)
     * @param elementDepth  The number of rows of nodes in each element (defaults to 1)
     * @param startAtBottomLeft  If true then the mesh starts in the bottom left corner
     *     of the domain rather than the centre, used for simple tests (defaults to false)
     * @param isPeriodicInX  If true then the mesh is periodic in the x dimension (defaults to false)
     * @param isPeriodicInY  If true then the mesh is periodic in the y dimension (defaults to false)
     * @param isPeriodicInZ  If true then the mesh is periodic in the y dimension (defaults to false)
     */
    PottsMeshGenerator(unsigned numNodesAcross,
                       unsigned numElementsAcross,
                       unsigned elementWidth,
                       unsigned numNodesUp=1u,
                       unsigned numElementsUp=1u,
                       unsigned elementHeight=1u,
                       unsigned numNodesDeep=1u,
                       unsigned numElementsDeep=1u,
                       unsigned elementDepth=1u,
                       bool startAtBottomLeft = false,
                       bool isPeriodicInX = false,
                       bool isPeriodicInY = false,
                       bool isPeriodicInZ = false);

    /**
     * Null constructor for derived classes to call.
     */
    PottsMeshGenerator()
    {
    }

    /**
     * Destructor - deletes the mesh object and pointer
     */
    virtual ~PottsMeshGenerator();

    /**
     * Helper method to calculate the Noore and Von Neumann Neighbourhoods of all nodes
     *
     * @param isPeriodicInX  If true then the mesh is periodic in the x dimnension
     */
    void CaclulateNeighbouringNodeIndices(bool isPeriodicInX);


    /**
     * @return a Cubouid or rectangulr Potts mesh.
     */
    virtual PottsMesh<DIM>* GetMesh();
};

#endif /*POTTSMESHGENERATOR_HPP_*/
