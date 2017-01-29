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
 * CellShape.cpp
 *
 *  Created on: 23 Jun 2016
 *      Author: Bartosz Jan Bartmanski
 */

#include "CellShape.hpp"

#include "Debug.hpp"

CellShape::CellShape(double nodeNumOrSpacing,
					 std::vector<c_vector<double, 2> > base,
					 double cellHeight,
					 bool num_or_spacing)
{
	// Check for correct input
	if ( nodeNumOrSpacing < 1 && num_or_spacing == false)
	{
		EXCEPTION("Number of nodes can't be a fraction");
	}
	if ( nodeNumOrSpacing > 1 && num_or_spacing == true)
	{
		EXCEPTION("Node spacing can't be an integer");
	}

	double node_spacing;
	unsigned numPoints;
	std::vector<unsigned> corner_indices;

	unsigned num_points = base.size();
	// Calculate the normals on both sides of the base, which second requires tangents

	// second side
	c_vector<double, 2> second_tangent = base[1] - base[0];
	second_tangent *= (1.0 / norm_2(second_tangent));
	c_vector<double, 2> second_normal;
	second_normal[0] = second_tangent[1];
	second_normal[1] = -1.0 * second_tangent[0];
	second_normal *= cellHeight;

	// first side
	c_vector<double, 2> first_tangent = base[num_points-1] - base[num_points - 2];
	first_tangent *= (1.0 / norm_2(first_tangent));
	c_vector<double, 2> first_normal;
	first_normal[0] = first_tangent[1];
	first_normal[1] = -1.0 * first_tangent[0];
	first_normal *= cellHeight;

	c_vector<double, 2> bisector = 0.5 * (second_normal + first_normal);
	bisector *= (1.0 / norm_2(bisector));

	double len_of_base = 0;
	for (unsigned i = 0; i < base.size()-1; i++)
	{
		len_of_base += norm_2(base[i+1] - base[i]);
	}

	// Calculate the vector of the top side of the polygon
	c_vector<double, 2> top = second_normal + (base[0] - base[num_points-1]) - first_normal;
	double perimeter = len_of_base + norm_2(second_normal) + norm_2(top) + norm_2(first_normal);

	if ( num_or_spacing == false )
	{
		numPoints = nodeNumOrSpacing;
		// Need to calculate the node spacing, so need the perimeter
		node_spacing = perimeter / numPoints;
	}
	else
	{
		node_spacing = nodeNumOrSpacing;
	}

	// Create the base of the cell, using the @param tangent
	// Add the corner index to the list of corner indices
	corner_indices.push_back(0);
	for ( unsigned i = 0; i < base.size(); i++)
	{
		mPoints.push_back(base[i] + 0.001 * bisector);
	}
	// Add the corner index to the list of corner indices
	corner_indices.push_back(mPoints.size()-1);
	// Now for the rest
	unsigned numPointsLeft = floor( norm_2(first_normal) / node_spacing);
	for ( unsigned i = 0; i < numPointsLeft; i++)
	{
		mPoints.push_back( mPoints.back() + (node_spacing / norm_2(first_normal)) * first_normal);
	}
	// Add the corner index to the list of corner indices
	corner_indices.push_back(mPoints.size()-1);

	unsigned numPointsTop = ceil( norm_2(top) / node_spacing );
	for ( unsigned i = 0; i < numPointsTop; i++)
	{
		mPoints.push_back( mPoints.back() + (node_spacing / norm_2(top)) * top);
	}
	// Add the corner index to the list of corner indices
	corner_indices.push_back(mPoints.size()-1);

	unsigned numPointsRight = floor( norm_2(second_normal) / node_spacing ) - 1;
	for ( unsigned i = 0; i < numPointsRight; i++)
	{
		mPoints.push_back( mPoints.back() - (node_spacing / norm_2(second_normal)) * second_normal);
	}

	mPerimeter = perimeter;
	mCornerIndices = corner_indices;
}

CellShape::~CellShape()
{
	mPoints.clear();
}

const std::vector<c_vector<double, 2> > CellShape::GetPointsAsVectors() const
{
	return mPoints;
}

const std::vector<ChastePoint<2> > CellShape::GetPointsAsChastePoints() const
{
	std::vector<ChastePoint<2> > chaste_points;

	for ( unsigned point = 0; point < mPoints.size(); point++ )
	{
		chaste_points.push_back(ChastePoint<2>(mPoints[point]));
	}

	return chaste_points;

}

double CellShape::GetPerimeter()
{
	return mPerimeter;
}

std::vector<unsigned> CellShape::GetCornerIndices()
{
	return mCornerIndices;
}


