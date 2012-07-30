/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef SEMMESHGENERATOR_HPP_
#define SEMMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "SemMesh.hpp"

/**
 * Generator of regular SemMesh.
 *
 * This takes in a request for a mesh with a given number of
 * cells in each direction, and an inital number of subcellular elements
 * per cell, and constructs cuboid cells in that configuration.
 */
template<unsigned DIM>
class SemMeshGenerator
{
protected:

    /** A pointer to the mesh this class creates */
    SemMesh<DIM>* mpMesh;

    /** The total number of subcellular elements per cell */
    unsigned mTotalSubcellularElementsPerCell;

    /** The equilibrium distance between nodes */
    double mNodeEquilibriumDistance;

public:

    /**
     * Constructor.
     *
     * Some garb here about how the size of each cell is determined.
     *
     * @param numCellsAcross the number of cells in the x-direction
     * @param numCellsUp the number of cells in the y-direction defaults to one for a string of cells.
     * @param numCellDeep the number of cells in the z-direction. defaults to one for a flat plane of cells.
     * @param numSubCellularElementsPerCellAcross the number of subcellular elements per cell along the x-axis, represented by Nodes defaults to 10
     * @param numSubCellularElementsPerCellUp the number of subcellular elements per cell along the y-axis, represented by Nodes defaults to 10
     * @param numSubCellularElementsPerCellDeep the number of subcellular elements per cell along the z-axis, represented by Nodes defaults to 0 for 2D.
     */
    SemMeshGenerator(   unsigned numCellsAcross,
                        unsigned numCellsUp = 1,
                        unsigned numCellsDeep = 1,
                        unsigned numSubCellularElementsPerCellAcross = 10,
                        unsigned numSubCellularElementsPerCellUp = 10,
                        unsigned numSubCellularElementsPerCellDeep = 1);

    /**
     * Null constructor for derived classes to call.
     */
    SemMeshGenerator()
    {
    }

    /**
     * Return the calculated equilibrium distance between two nodes.
     *
     * @return mNodeEquilibriumDistance
     */
    double GetEquilibriumDistance();

    /**
     * Destructor - deletes the mesh object and pointer
     */
    virtual ~SemMeshGenerator();

    /**
     * @return a Cubouid or rectangulr Potts mesh.
     */
    virtual SemMesh<DIM>* GetMesh();
};

#endif /*SEMMESHGENERATOR_HPP_*/
