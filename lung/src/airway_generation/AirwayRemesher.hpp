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


#ifndef AIRWAY_REMESHER_HPP_
#define AIRWAY_REMESHER_HPP_

#include "MutableMesh.hpp"
#include "AirwayTreeWalker.hpp"
#include "AirwayPropertiesCalculator.hpp"

/**
 * Class to allow rebalancing of airway tree meshes to improve condition number.
 */
class AirwayRemesher
{
public:
    /**
     * Constructor
     *
     * @param rMesh The centrelines mesh to remesh.
     * @param rootIndex The root node index corresponding to the trachea.
     */
    AirwayRemesher(TetrahedralMesh<1,3>& rMesh,
                   unsigned rootIndex);

    /**
     * Creates a remeshed version of the underlying mesh.
     *
     * This version attempts to balance the mesh in such a way as to minimise the
     * condition number matrices resulting from Poiseuille based network flow problems.
     *
     * @param rOutputMesh The mesh object to be written to.
     * @param maximumResistance The maximum allowed resistance of an element
     */
    void Remesh(MutableMesh<1,3>& rOutputMesh, double maximumResistance);

    /**
     * Creates a remeshed version of the underlying mesh.
     *
     * This version simply removes all intermediate, non bifurcation, nodes.
     *
     * @param rOutputMesh The mesh object to be written to.
     * @param maximumResistance The maximum allowed resistance of an element
     */
    void Remesh(MutableMesh<1,3>& rOutputMesh);

private:
    /**
     * A mesh containing the major airways
     */
    TetrahedralMesh<1,3>& mrMesh;

    /**
     * The index of the root of the airway tree (trachea)
     */
    unsigned mOutletNodeIndex;

    /**
     * Used to navigate through the airways mesh.
     */
    AirwayTreeWalker mWalker;

    /**
     * Used to calculate order on the airways mesh.
     */
    AirwayPropertiesCalculator mCalculator;
};

#endif //  AIRWAY_REMESHER_HPP_
