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

/*
 * Quadrilateral.hpp
 *
 *  Created on: 29 May 2016
 *      Author: Bartosz Jan Bartmanski
 */

#ifndef PROJECTS_IMMERSEDBOUNDARY_BJB_TEST_QUADRILATERAL_HPP_
#define PROJECTS_IMMERSEDBOUNDARY_BJB_TEST_QUADRILATERAL_HPP_

#include <cmath>
#include <vector>

#include "ChastePoint.hpp"
#include "Exception.hpp"
#include "Node.hpp"

class Quadrilateral
{

private:

	/** The target spacing between points. */
	double mTargetNodeSpacing;

	/** Vector to store the points. */
	std::vector<c_vector<double, 2> > mPoints;

	/** The perimeter of the trapezoid */
	double mPerimeter;

	/** The indices of the corners of the quadrilateral */
	std::vector<unsigned> mCornerIndices;

public:

	/**
	 * Constructor
	 */
    Quadrilateral(double nodeNumOrSpacing,
                  c_vector<double, 2> start_point,
                  c_vector<double, 2> tangent,
                  c_vector<double, 2> normal_right,
                  c_vector<double, 2> normal_left,
                  bool num_or_spacing);

    /**
     * Destructor - clears mPoints
     */
    virtual ~Quadrilateral();

    /**
    * @return #mPoints
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
    * @return the perimeter of the trapezoid
    */
    double GetPerimeter();

    /**
    * @return the indices of the corners
    */
    std::vector<unsigned> GetCornerIndices();

};



#endif /* PROJECTS_IMMERSEDBOUNDARY_BJB_TEST_QUADRILATERAL_HPP_ */
