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

#ifndef CYLINDRICALHONEYCOMBMESHGENERATOR_HPP_
#define CYLINDRICALHONEYCOMBMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "HoneycombMeshGenerator.hpp"
#include "Cylindrical2dMesh.hpp"

/**
 * Honeycomb mesh generator that creates a 2D "cylindrical" mesh (one in which
 * periodicity is imposed on the left and right boundaries) for use in cell-centre
 * simulations.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 */
class CylindricalHoneycombMeshGenerator : public HoneycombMeshGenerator
{
public:

    /**
     * Default constructor.
     *
     * @param numNodesAlongWidth  The number of cells you want alopng the bottom of the domain
     * @param numNodesAlongLength  The number of cells you want sides of the domain
     * @param ghosts  The thickness of ghost nodes to put around the edge (defaults to 3)
     * @param scaleFactor  The scale factor for the width (circumference) of the cells (defaults to 1.0)
     */
    CylindricalHoneycombMeshGenerator(unsigned numNodesAlongWidth, unsigned numNodesAlongLength, unsigned ghosts=3, double scaleFactor=1.0);

    /**
     * @return a 2D honeycomb mesh
     */
    MutableMesh<2,2>* GetMesh();

    /**
     * @return a 2D honeycomb mesh with periodic left/right boundaries
     */
    Cylindrical2dMesh* GetCylindricalMesh();
};

#endif /*CYLINDRICALHONEYCOMBMESHGENERATOR_HPP_*/
