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

#ifndef GEODESICSPHERE23GENERATOR_HPP_
#define GEODESICSPHERE23GENERATOR_HPP_

#include <vector>
#include "MutableVertexMesh.hpp"
#include "Node.hpp"
#include "VertexElement.hpp"

/**
 * A generator that creates a geodesic sphere as a two-dimensional triangulated
 * surface embedded in 3D space (a mesh of VertexElement<2,3> faces).
 *
 * The generator is initialised (in the constructor) with a regular icosahedron
 * inscribed in the unit sphere: 12 nodes, 30 edges and 20 triangular faces. Each
 * call to SubDivide() refines the surface by splitting every triangular face into
 * four smaller triangles, with the newly created edge-midpoint nodes projected back
 * onto the unit sphere. Repeated subdivision therefore yields successively finer,
 * approximately uniform, triangulations of the sphere.
 *
 * GetDual() returns the dual polyhedron of the current triangulation (a Goldberg
 * polyhedron, whose faces are hexagons and twelve pentagons) as a MutableVertexMesh.
 * This dual is useful as the apical (or basal) surface of a closed monolayer of
 * cells; see MonolayerVertexMeshGenerator.
 */
class GeodesicSphere23Generator
{
    friend class TestMonolayerVertexMeshGenerator;

protected:
    /** The nodes (vertices) of the triangulated sphere, lying on the unit sphere. */
    std::vector<Node<3>*> mNodes;

    /** The edges of the triangulated sphere, each a VertexElement<1,3>. */
    std::vector<VertexElement<1, 3>*> mEdges;

    /** The triangular faces of the sphere, each a VertexElement<2,3>. */
    std::vector<VertexElement<2, 3>*> mFaces;

public:
    /**
     * Constructor. Populates mNodes, mEdges and mFaces with a regular icosahedron
     * (12 nodes, 30 edges, 20 triangular faces) inscribed in the unit sphere.
     */
    GeodesicSphere23Generator();

    /**
     * Refine the current triangulation by subdividing each triangular face into four,
     * inserting a new node at the midpoint of each edge and projecting it onto the unit
     * sphere. May be called repeatedly to obtain finer triangulations.
     */
    void SubDivide();

    /**
     * @return a pointer to the dual of the current triangulation (a Goldberg polyhedron
     *         of hexagonal and pentagonal faces), as a new MutableVertexMesh<2,3>. The
     *         caller is responsible for deleting the returned mesh.
     */
    MutableVertexMesh<2, 3>* GetDual();
};

#endif /* GEODESICSPHERE23GENERATOR_HPP_ */
