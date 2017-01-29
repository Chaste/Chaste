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
 * ImmersedBoundaryCircularBoundaryGenerator.cpp
 *
 *  Created on: 2 May 2016
 *      Author: bartosz
 */

#include "ImmersedBoundaryMeshGeneratorWithCircularMembrane.hpp"
#include "ImmersedBoundaryEnumerations.hpp"
#include "RandomNumberGenerator.hpp"

#include "Debug.hpp"

ImmersedBoundaryMeshGeneratorWithCircularMembrane::ImmersedBoundaryMeshGeneratorWithCircularMembrane(unsigned numCells,
																									 unsigned numNodesPerCell,
																									 double ellipseExponent,
																									 double cellAspectRatio,
																									 double radius)
	: mpMesh(NULL)
{
	// Check for sensible input
	assert(numCells > 0);
	assert(numNodesPerCell > 3);
	assert(ellipseExponent > 0.0);
	assert(cellAspectRatio > 0.0);
	assert(radius < 0.5);

	// Helper vectors
	unit_vector<double> x_unit(2,0);
	unit_vector<double> y_unit(2,1);
	c_vector<double, 2> center = 0.5 * x_unit + 0.5 * y_unit;
	c_vector<double, 2> zero = 0.0 * x_unit + 0.0 * y_unit;

	// Calculating cell width and cell height based on the radius of the circular membrane
	double cell_width = M_PI * radius / numCells;
	double cell_height = cellAspectRatio * cell_width;

	// Make sure that the cells are not too large
	if (cell_height > 0.2)
	{
		cell_height = 0.2;
		cell_width = cell_height / cellAspectRatio;
	}

	// Checking that the membrane and the cells will fit into the simulation box
	if ( radius + cell_height > 0.5)
	{
		EXCEPTION("Cells don't fit in");
	}

	// Calculate the number of nodes in the base of a cell
	double node_spacing = 2 * ( cell_width + cell_height ) / numNodesPerCell;
	unsigned num_nodes_base = ceil(cell_width / node_spacing) + 1;

	// This is to ensure that the nodes spacing in the membrane and the cells is the same
	unsigned mem_num_nodes =  ceil( (numCells * numNodesPerCell) / ( 1.0 + cellAspectRatio) );

	// Generating a super-ellipse to get the locations of the nodes of the membrane
	// TODO: Check that the Superellipse generator still works or use Fergus' one!
	SuperellipseGenerator* p_circle_gen = new SuperellipseGenerator(mem_num_nodes, 1.0, 2*radius, 2*radius, 0.5-radius, 0.5-radius);
	std::vector<c_vector<double, 2> > circle_locations = p_circle_gen->GetPointsAsVectors();
	/*double length = p_circle_gen->GetTotalArcLength();*/ // May need to reimplement this function
	/*delete p_circle_gen;*/ // TODO: not sure whether this is necessary

	/* Create a vector of IB elements and a vector of nodes
	 * ib_elements stores the elements of the simulation - membrane and cells
	 * nodes stores all of the nodes for this simulation
	 */
	std::vector<ImmersedBoundaryElement<2,2>*> ib_elements;
	std::vector<Node<2>*> nodes;

	// Creating an IB element to represent the membrane - for some reason it's not closed, so will have to ask Fergus about it
	std::vector<Node<2>*> nodes_this_elem;
	for (unsigned mem_node_idx=0; mem_node_idx<circle_locations.size(); mem_node_idx++)
	{
		nodes.push_back(new Node<2>(mem_node_idx, circle_locations[mem_node_idx], true));
		nodes_this_elem.push_back(nodes.back());
	}
	ib_elements.push_back(new ImmersedBoundaryElement<2,2>(0, nodes_this_elem));

	// Generating a super-ellipse to get the locations of the nodes of the reference cell
	SuperellipseGenerator* p_gen = new SuperellipseGenerator(numNodesPerCell, ellipseExponent, cell_width, cell_height, -0.5*cell_width, -0.5*cell_height);
	std::vector<c_vector<double, 2> > locations = p_gen->GetPointsAsVectors();
	delete p_gen;

	// Need an angle by which this cell will be rotated by
	double angle = M_PI_2;
	// Angle corresponding to the arc on which a cell lies
	// double theta = M_PI / numCells;
	double theta = num_nodes_base * (2.0 * M_PI / circle_locations.size());
	// A small quantity, which will be used to move
	double epsilon=0.01;

	angle -= 0.5*theta;

	// Place each cell to a correct location around the membrane
	for ( unsigned cell_index = 0; cell_index < numCells; cell_index++)
	{
		// Set things up for the spatial transformations

		// Calculate the index of the membrane node next to which the cell will be placed
		unsigned mem_idx = floor(0.5 * circle_locations.size()) + floor((cell_index-0.5) * num_nodes_base);
		// Check that it hasn't gone out of scope and if it has jump back to the beginning
		if (mem_idx > circle_locations.size())
		{
			mem_idx = fmod(mem_idx, circle_locations.size());
		}
		// Get the location of this membrane node
		c_vector<double,2> mem_node_loc = circle_locations[mem_idx];
		// Helper vector dependent on the membrane node
		c_vector<double,2> correction =(mem_node_loc - center);
		// Normalise this vector
		double norm = sqrt(correction[0]*correction[0] + correction[1]*correction[1]);
		correction *= 1.0 / norm;
		// Collect all the necessary vectors into one
		c_vector<double,2> translate = mem_node_loc+ ((0.5 * cell_height) + epsilon)*correction;
		// Vector with the locations of nodes of a transformed cell
		std::vector<c_vector<double,2> > new_locations;
		new_locations.resize(locations.size());
		// Do all the transformations on each cell
		for ( unsigned node_idx = 0; node_idx < numNodesPerCell; node_idx++)
		{
			new_locations[node_idx] = zero;
			// Rotate the angle by a small amount
			new_locations[node_idx][0] = ((double)cos(angle)) * locations[node_idx][0] - ((double)sin(angle)) * locations[node_idx][1];
			new_locations[node_idx][1] = ((double)sin(angle)) * locations[node_idx][0] + ((double)cos(angle)) * locations[node_idx][1];
			// Shrink the cell by a small amount - may not be necessary, but just to be safe
			new_locations[node_idx] *= 0.95;
			// Translate each node of the cell to the correct position
			new_locations[node_idx] += translate;
		}
		// Need a vector to store all the nodes, hence:
		std::vector<Node<2>*> nodes_this_elem;
		// Pass on the node locations to the ib_elements and nodes containers
		for ( unsigned node_idx = 0; node_idx < numNodesPerCell; node_idx++)
		{
			// Size of the nodes vector, which becomes the index of the new node being added
			unsigned node_index = nodes.size();
			nodes.push_back(new Node<2>(node_index, new_locations[node_idx], true));
			nodes_this_elem.push_back(nodes.back());
		}
		ib_elements.push_back(new ImmersedBoundaryElement<2,2>(ib_elements.size(), nodes_this_elem));
		// Update the angle of the rotation
		angle += theta;
	}

	mpMesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_elements);
}

ImmersedBoundaryMeshGeneratorWithCircularMembrane::~ImmersedBoundaryMeshGeneratorWithCircularMembrane()
{
	delete mpMesh;
}
ImmersedBoundaryMesh<2,2>* ImmersedBoundaryMeshGeneratorWithCircularMembrane::GetMesh()
{
	return mpMesh;
}

