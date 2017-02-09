/*

Copyright (c) 2005-2017, University of Oxford.
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

// Created by bartosz on 09/02/17.

#ifndef CHASTE_IMMERSEDBOUNDARYDROSOPHILAEMBRYOCROSSSECTIONMESHGENERATOR_HPP
#define CHASTE_IMMERSEDBOUNDARYDROSOPHILAEMBRYOCROSSSECTIONMESHGENERATOR_HPP

#include <cmath>
#include <vector>

#include "ImmersedBoundaryMesh.hpp"
#include "SuperellipseGenerator.hpp"
#include "UblasCustomFunctions.hpp"

/**
 * Creates a collection of immersed boundary elements to model a cross section of a drosophila embryo as an array of cells
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 * NOTE: the user should change mesh parameters like fluid grid spacing, as these are not altered from defaults here.
 */

class ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator
{
protected:
    /** A pointer to the mesh this class creates. */
    ImmersedBoundaryMesh<2,2>* mpMesh;

    /** Calculates the perimeter of the membrane */
    double CalculateTotalArcLengthOfSuperellipse(double membraneWidth, double membraneHeight, double ellipseExponent);

    /** Calculates the normal to the two points specified */
    c_vector<double, 2> CalculateNormal(c_vector<double, 2> point_1, c_vector<double, 2> point_2);

    /** @returns a set of points that represent a cell */
    const std::vector<c_vector<double, 2> > BuildACell(double targetNodeSpacing,
                                                       std::vector<c_vector<double, 2> > base,
                                                       c_vector<double, 2> normal_1,
                                                       c_vector<double, 2> normal_2,
                                                       double cell_height) const;

public:

    /**
     * Default constructor
     * @param numCells - number of cells
     * @param numNodesPerCell - number of nodes per cell
     * @param radius - inner radius of the ring of cells
     * @param cellHeight - height of the cells
     */
    ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator(unsigned numCells=10,
                                                              unsigned numNodesPerCell=100,
                                                              double cellHeight=0.1,
                                                              double radius=0.2,
                                                              double aspectRatio=1.0,
                                                              double ellipseExponent=1.0);

    /** Null constructor - for derived classes to call */
    ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator()
    {
    };

    /**
     * Destructor
     *
     * Deletes the mesh and the pointer to it
     */
    virtual ~ImmersedBoundaryDrosophilaEmbryoCrossSectionMeshGenerator();

    /**
     * @return a 2D honeycomb mesh based on a 2D plane
     */
    ImmersedBoundaryMesh<2,2>* GetMesh();

};


#endif //CHASTE_IMMERSEDBOUNDARYDROSOPHILAEMBRYOCROSSSECTIONMESHGENERATOR_HPP
