/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef IMMERSEDBOUNDARYPALISADEMESHGENERATOR_HPP_
#define IMMERSEDBOUNDARYPALISADEMESHGENERATOR_HPP_

#include <cmath>
#include <vector>

#include "ImmersedBoundaryMesh.hpp"
#include "SuperellipseGenerator.hpp"

/**
 * Creates a palisade of immersed boundary elements.
 *
 * NOTE: the user should delete the mesh after use to manage memory.
 * NOTE: the user should change mesh parameters like fluid grid spacing, as these are not altered from defaults here.
 */
class ImmersedBoundaryPalisadeMeshGenerator
{
protected:

    /** A pointer to the mesh this class creates. */
    ImmersedBoundaryMesh<2,2>* mpMesh;

    /** The number of cells from left to right. */
    unsigned mNumCellsWide;

    /** The number of nodes per cell. */
    unsigned mNumNodesPerCell;

    /** The exponent of the superellipse. */
    double mEllipseExponent;

    /** The aspect ratio of each cell (height / width). */
    double mCellAspectRatio;

    /** The random variation in y_pos. */
    double mRandomYMult;

    /** Whether the mesh has a basement membrane. */
    bool mMembrane;

public:

    /**
     * Default constructor.
     *
     * @param numCellsWide  the number of cells from left to right along the domain
     * @param numNodesPerCell  the number of nodes per cell (defaults to 100)
     * @param ellipseExponent  the exponent of the superellipse (defaults to 0.2)
     * @param cellAspectRatio  the aspect ratio of each cell (defaults to 2)
     * @param randomYMult  the random variation in y_pos (defaults to 0)
     * @param membrane  whether the mesh has a basement membrane (defaults to false)
     */
    ImmersedBoundaryPalisadeMeshGenerator(unsigned numCellsWide,
                                          unsigned numNodesPerCell=100,
                                          double ellipseExponent=0.2,
                                          double cellAspectRatio=2.0,
                                          double randomYMult=0.0,
                                          bool membrane=false);

    /**
     * Null constructor for derived classes to call.
     */
    ImmersedBoundaryPalisadeMeshGenerator()
    {
    }

    /**
     * Destructor.
     *
     * Deletes the mesh object and pointer.
     */
    virtual ~ImmersedBoundaryPalisadeMeshGenerator();

    /**
     * @return a 2D honeycomb mesh based on a 2D plane
     */
    ImmersedBoundaryMesh<2,2>* GetMesh();

    ///\todo document this method
    void SetRandomYMult(double mult);
};

#endif /*IMMERSEDBOUNDARYPALISADEMESHGENERATOR_HPP_*/
