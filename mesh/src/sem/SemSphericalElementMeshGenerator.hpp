/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef SEMSPHERICALELEMENTMESHGENERATOR_HPP_
#define SEMSPHERICALELEMENTMESHGENERATOR_HPP_


#include <vector>
#include <memory>

#include "SemEnumerations.hpp"
#include "SemMesh.hpp"
#include "UblasVectorInclude.hpp"


/**
 * Sem mesh generator that creates a single approximately spherical SEM element, in either 1D,
 * 2D, or 3D. This is the geometry of the equilibrated single cell of Sandersius & Newman (2008)
 * Section 1.2, in which an aggregate of N subcellular elements forms a roughly spherical cell of
 * radius R_cell.
 *
 * Where SemSingleElementMeshGenerator takes a node count per dimension and produces a cuboid,
 * this generator takes a total node count and a cell radius, which is the parameterisation the
 * paper's scaling laws are written in (see SemParameterScaling).
 *
 * The cell is carved out of a close-packed lattice by keeping the requested number of lattice
 * sites nearest the centre, then scaling so that the outermost node lies exactly at the requested
 * radius. This gives exactly the requested number of nodes with a nearest-neighbour distance
 * equal to the node spacing everywhere.
 *
 * Nodes are labelled SEM_BOUNDARY_REGION when they lie within a given thickness of the cell
 * surface, and SEM_INTERIOR_REGION otherwise; see IsBoundaryNode(). A thickness of one node
 * spacing, the default, gives a single-node-thick cortex, matching the outermost grid layer that
 * SemSingleElementMeshGenerator labels.
 */
template <unsigned DIM> class SemSphericalElementMeshGenerator
{
private:
    /** A pointer to the mesh this class creates */
    std::shared_ptr<SemMesh<DIM>> mpMesh;

    /** Total number of subcellular nodes in the element */
    unsigned mNumNodes;

    /** The radius of the element, i.e. the distance from its centre to its outermost node */
    double mCellRadius;

    /** Thickness of the boundary region, in units of the node spacing */
    double mBoundaryThickness;

    /** The node spacing resulting from fitting mNumNodes nodes into a cell of radius mCellRadius */
    double mNodeSpacing;

    /**
     * Generate node positions filling a ball centred on the origin, and set mNodeSpacing to the
     * resulting nearest-neighbour distance.
     *
     * @return vector of node coordinates (each as c_vector<double, DIM>)
     */
    std::vector<c_vector<double, DIM>> GenerateNodePositions();

    /**
     * Generate the mesh from the node positions.
     *
     * @param positions vector of node coordinates produced by GenerateNodePositions()
     */
    void GenerateMesh(std::vector<c_vector<double, DIM>> positions);

    /**
     * Return whether the node at the given position is in the boundary region, i.e. whether it
     * lies within mBoundaryThickness node spacings of the cell surface. Region 0 = interior,
     * 1 = boundary.
     *
     * The generator knows the cell's shape analytically, so the surface is found by a radial
     * test rather than by reconstructing a hull from the point cloud.
     *
     * @param rPosition the position of the node, relative to the centre of the element
     * @return true if the node lies in the boundary region, false otherwise
     */
    bool IsBoundaryNode(const c_vector<double, DIM>& rPosition) const;

public:
    /**
     * Constructor.
     *
     * @param numNodes the total number of subcellular nodes, which must be at least two
     * @param cellRadius the distance from the centre of the element to its outermost node
     * @param boundaryThickness thickness of the boundary region in units of the node spacing
     */
    SemSphericalElementMeshGenerator(unsigned numNodes, double cellRadius = 0.5,
        double boundaryThickness = 1.0);

    /**
     * Default constructor
     */
    SemSphericalElementMeshGenerator() = default;

    /**
     * Default destructor
     */
    virtual ~SemSphericalElementMeshGenerator() = default;

    /**
     * @return a shared view into the mesh
     */
    virtual std::shared_ptr<SemMesh<DIM>> GetMesh();

    /**
     * The distance between neighbouring subcellular nodes, which follows from fitting the
     * requested number of nodes into a cell of the requested radius. This is the r_eq of
     * Sandersius & Newman (2008) Section 1.2.
     *
     * @return the node spacing
     */
    double GetNodeSpacing() const;
};

#endif /*SEMSPHERICALELEMENTMESHGENERATOR_HPP_*/
