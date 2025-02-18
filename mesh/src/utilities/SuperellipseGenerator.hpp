/*

Copyright (c) 2005-2025, University of Oxford.
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
#include "Exception.hpp"
#include "Node.hpp"

/**
 * Class to generate roughly equally spaced points around a 2D superellipse.
 * This shape is described by the equation (x/a)^n + (y/b)^n = 1, where n, a
 * and b are positive numbers.
 */
class SuperellipseGenerator
{
private:

    /** The target spacing between points. */
    double mTargetNodeSpacing;

    /**
     * The height at which the point of maximal curvature is achieved in the
     * top right-hand corner of the superellipse.
     *
     * If the superellipse exponent is 1.0 this definition does not make sense
     * and we choose 80% of the height of the cell.
     */
    double mHeightOfTopSurface;

    /** Vector to store the points. */
    std::vector<c_vector<double, 2> > mPoints;

public:

    /**
     * Constructor.
     *
     * Default values of the input arguments produce a circle bounded by the unit square.
     *
     * @param numPoints  The number of points to generate around the circumference (defaults to 100)
     * @param ellipseExponent  The exponent of the superellipse (defaults to 1)
     * @param width  The width of the superellipse (defaults to 1)
     * @param height The height of the superellipse (defaults to 1)
     * @param botLeftX The x-coordinate of the bottom-left of the superellipse bounding box (defaults to 0)
     * @param botLeftY The y-coordinate of the bottom-left of the superellipse bounding box (defaults to 0)
     */
    SuperellipseGenerator(unsigned numPoints=100,
                          double ellipseExponent=1.0,
                          double width=1.0,
                          double height=1.0,
                          double botLeftX=0.0,
                          double botLeftY=0.0);

    /**
     * Destructor.
     *
     * Clears mPoints.
     */
    virtual ~SuperellipseGenerator();

    /**
     * @return #mPoints
     */
    double GetTargetNodeSpacing();

    /**
     * @return #mHeightOfTopSurface
     */
    double GetHeightOfTopSurface();

    /**
     * @return #mPoints
     */
    const std::vector<c_vector<double, 2> > GetPointsAsVectors() const;

    /**
     * @return #mPoints as a vector of ChastePoint objects
     */
    const std::vector<ChastePoint<2> > GetPointsAsChastePoints() const;
};

#endif /*SUPERELLIPSEGENERATOR_HPP_*/
