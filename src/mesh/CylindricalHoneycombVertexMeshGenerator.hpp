/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef CYLINDRICALHONEYCOMBVERTEXMESHGENERATOR_HPP_
#define CYLINDRICALHONEYCOMBVERTEXMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "Cylindrical2dVertexMesh.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

/**
 *  Generator of honeycomb meshes, used as starting points for many simulations.
 *
 *  This class takes in options such as width, height, number of ghost nodes
 *  and generates a honeycomb mesh (with equal distance between nodes), and ghost
 *  node information when requested.
 *
 *  NOTE: the user should delete the mesh after use to manage memory.
 */
class CylindricalHoneycombVertexMeshGenerator : HoneycombVertexMeshGenerator
{

public:

    /**
     * Constructor.
     *
     * @param numElementsAcross  The number of columns of elements in the mesh
     * @param numElementsUp  The number of rows of elements in the mesh
     * @param isFlatBottom  Whether to enforce a flat bottom to the mesh (defaults to false; only used if isCylindrical is true)
     * @param cellRearrangementThreshold the minimum threshold distance for element rearrangment (defaults to 0.01)
     * @param edgeDivisionThreshold the maximum threshold distance for edge division (defaults to DBL_MAX)
     * @param t2Threshold the maximum threshold distance for Type 2 swaps (defaults to 0.001)
     */
    CylindricalHoneycombVertexMeshGenerator(unsigned numElementsAcross,
                                 unsigned numElementsUp,
                                 bool isFlatBottom=false,
                                 double cellRearrangementThreshold=0.01,
                                 double edgeDivisionThreshold=DBL_MAX,
                                 double t2Threshold=0.001);
    /**
     * Destructor - deletes the mesh object and pointer
     */
//    ~CylindricalHoneycombVertexMeshGenerator();

    /**
     * @return a honeycomb mesh based on a 2D plane.
     */
    MutableVertexMesh<2,2>* GetMesh();

    /**
     * @return a honeycomb mesh with cylindrical boundaries.
     */
    Cylindrical2dVertexMesh* GetCylindricalMesh();
};

#endif /*CYLINDRICALHONEYCOMBVERTEXMESHGENERATOR_HPP_*/
