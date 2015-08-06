/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef SUPERELLIPSEGENERATOR_HPP_
#define SUPERELLIPSEGENERATOR_HPP_

#include <cmath>
#include <vector>
#include "ChastePoint.hpp"
#include "Node.hpp"
#include "ImmersedBoundaryElement.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "Exception.hpp"

/**
 * Class to generate roughly equally spaced points around a 2D superellipse (x/a)^n + (y/b)^n = 1
 */
class SuperellipseGenerator
{

private:

    /** The target node spacing */
    double mTargetNodeSpacing;

    /** Vector to store the points */
    std::vector<c_vector<double, 2> > mPoints;

    /** Immersed boundary mesh that can optionally be returned */
    ImmersedBoundaryMesh<2,2>* mpMesh;

public:

    /**
     * Constructor.
     *
     * Default values produce a circle bounded by the unit square
     *
     * @param numPoints  The number of points to generate around the circumference
     * @param ellipseExponent  The exponent of the superellipse
     * @param width  The width of the superellipse
     * @param height The height of the superellipse
     * @param botLeftX The x-coordinate of the bottom-left of the superellipse bounding box
     * @param botLeftY The y-coordinate of the bottom-left of the superellipse bounding box
     */
    SuperellipseGenerator(unsigned numPoints=100,
                          double ellipseExponent=1.0,
                          double width=1.0,
                          double height=1.0,
                          double botLeftX=0.0,
                          double botLeftY=0.0);

    /**
     * Destructor - clears mPoints.
     */
    virtual ~SuperellipseGenerator();

    /**
     * @return the target node spacing
     */
    double GetTargetNodeSpacing();

    /**
     * @return vector of c_vector objects
     */
    const std::vector<c_vector<double, 2> > GetPointsAsVectors() const;

    /**
     * @return vector of ChastePoint objects
     */
    const std::vector<ChastePoint<2> > GetPointsAsChastePoints() const;

    /**
     * @return an ImmersedBoundaryMesh based on generated locations
     */
    ImmersedBoundaryMesh<2,2>* GetMesh();
};

#endif /*SUPERELLIPSEGENERATOR_HPP_*/
