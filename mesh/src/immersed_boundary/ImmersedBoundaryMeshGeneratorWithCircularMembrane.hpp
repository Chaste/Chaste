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
 * ImmersedBoundaryCircularBoundaryGenerator.hpp
 *
 *  Created on: 2 May 2016
 *      Author: bartosz
 */

#ifndef IMMERSEDBOUNDARYMESHGENERATORWITHCIRCULARMEMBRANE_HPP_
#define IMMERSEDBOUNDARYMESHGENERATORWITHCIRCULARMEMBRANE_HPP_

#include <cmath>
#include <vector>

#include "ImmersedBoundaryMesh.hpp"
// Might as well use Fergus' code to generate superellipse to get a circular boundary for my simulation
#include "SuperellipseGenerator.hpp"


// Creates a collection of cells on a circular membrane within the Immersed Boundary framework.

class ImmersedBoundaryMeshGeneratorWithCircularMembrane
{
protected:

	// A pointer to the mesh this class creates
	ImmersedBoundaryMesh<2,2>* mpMesh;

public:

	/* The default constructor
	 *
	 * All the parameters are as described above
	 */
	ImmersedBoundaryMeshGeneratorWithCircularMembrane(unsigned numCellWide,
													  unsigned numNodesPerCell=100,
													  double ellipseExponent=0.2,
													  double cellAspectRatio=2.0,
													  double radius=0.75);


	// Null constructor for the derived classes to call
	ImmersedBoundaryMeshGeneratorWithCircularMembrane()
	{
	}

	// Destructor - deletes the mesh object and pointer
	virtual ~ImmersedBoundaryMeshGeneratorWithCircularMembrane();

	// Return a 2D honeycomb mesh based on a 2D plane
	ImmersedBoundaryMesh<2,2>* GetMesh();

};

#endif /* PROJECTS_IMMERSEDBOUNDARY_BJB_SRC_IMMERSEDBOUNDARYMESHGENERATORWITHCIRCULARMEMBRANE_HPP_ */
