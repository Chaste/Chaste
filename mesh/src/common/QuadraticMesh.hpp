/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef QUADRATICMESH_HPP_
#define QUADRATICMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <map>
#include <vector>

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

/**
 * A concrete quadratic mesh class that inherits from TetrahedralMesh.
 */
template<unsigned DIM>
class QuadraticMesh : public TetrahedralMesh<DIM, DIM>
{
    friend class TestQuadraticMesh;
protected:

    /** Number of vertices, ie non-internal (non-quadratic), nodes. */
    unsigned mNumVertices;

    /**
     * Count nodes which are vertices (not marked as internal)
     */
    void CountVertices();

    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the mesh.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<TetrahedralMesh<DIM, DIM> >(*this);
    }


    /**
     * Create a quadratic mesh on the interval [0,numElemX] with numElemX elements in each
     * direction.  This is private, users should call ConstructRegularSlabMesh();
     *
     * Very badly named. This creates a QUADRATIC mesh, the linear just refers to the fact the
     * mesh is a line in 1D. The name is inherited from the parent class.
     *
     * @param numElemX Number of elements in x-direction (also, the width of the final mesh)
     */
    void ConstructLinearMesh(unsigned numElemX);


    /**
     * Create a quadratic mesh on a rectangle from (0,0) to (numElemX,numElemY)
     * with that number of elements in each direction.
     *
     * The method overloads the equivalent method in
     * AbstractTetrahedralMesh. This is private, users should call ConstructRegularSlabMesh();
     *
     * @param numElemX Number of elements in x-direction (also, the width of the final mesh)
     * @param numElemY Number of elements in y-direction (also, the height of the final mesh)
     * @param stagger is the same as the over-ridden method with same name in AbstractTetrahedralMesh.
     *        stagger = false gives all back-slash diagonals
     *        stagger = true check-boards the diagonals as forward and back slashes
     */
    void ConstructRectangularMesh(unsigned numElemX, unsigned numElemY, bool stagger=true);


    /**
     * Create a quadratic mesh on a cuobid from (0,0,0) to (numElemX,numElemY,numElemZ)
     * with that number of elements in each direction.
     * The method overloads the equivalent method in
     * AbstractTetrahedralMesh. This is private, users should call ConstructRegularSlabMesh();
     *
     * @param numElemX Number of elements in x-direction (also, the width of the final mesh)
     * @param numElemY Number of elements in y-direction (also, the height of the final mesh)
     * @param numElemZ Number of elements in y-direction (also, the depth of the final mesh)
     */
    void ConstructCuboid(unsigned numElemX, unsigned numElemY, unsigned numElemZ);

    /**
     * A helper method used in the private structured mesh constructors (ConstructRectangularMesh etc).
     *
     * Method uses top (and the zero vector) to calculate if the node should be designated as a boundary node.
     * Method uses top to determine if the node is outside the cuboid -- this allows for simpler loops in the caller
     * Method creates node, pushes node onto mNodes and marks it as an internal node.
     * @return pointer to new node
     * @param rIndex  is the index in the mesh which this node should take. Note: this method increments rIndex
     * @param rLocation  the position of the node in space (coordinates should be integers or multiples of 1/2)
     * @param rTop  the position of top-most node in the line/slab/cuboid
     */
    Node<DIM>* MakeNewInternalNode(unsigned& rIndex, c_vector<double, DIM>& rLocation, c_vector<double, DIM>& rTop);

    /**
     * A helper method used in the private structured mesh constructors (ConstructRectangularMesh etc).
     *
     * Gets the internal node index between two vertex node indices assuming ordered pairs have
     * been used as keys in the map
     *
     * @param globalIndex1  is the index of one of the vertex nodes
     * @param globalIndex2  is the index of the other vertex node (ordering is unimportant)
     * @param rEdgeMap  the map from ordered pairs of vertex indices to internal node index
     * @return  global node index of the internal node between globalIndex1 and globalIndex2
     */
    unsigned LookupInternalNode(unsigned globalIndex1, unsigned globalIndex2, std::map<std::pair<unsigned, unsigned>, unsigned>& rEdgeMap);


public:

    /**
     * Constructor
     *
     */
    QuadraticMesh()
    {
        this->mMeshIsLinear=false;
    }


    /**
     * Create a quadratic mesh on a slab (on a line in 1D, rectangle in 2d, cuboid in 3D),
     * with the given widths and given space step. In 1D height and depth need to be passed in
     * as 0 (the default value); in 2D depth must be zero
     *
     * @param spaceStep The spatial stepsize
     * @param width the width of the cuboid
     * @param height the height of the cuboid
     * @param depth the depth of the cuboid
     */
    QuadraticMesh(double spaceStep, double width, double height=0, double depth=0);


    /**
     * Load a quadratic mesh from a file.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader);

    /**
     * Load a quadratic mesh from a linear mesh file.
     *
     * Constructs as linear mesh, then exports to triangle/tetgen
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromLinearMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader);

    /**
     *  @return the number of vertices, ie non-internal (non-quadratic), nodes.
     */
    unsigned GetNumVertices() const;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(QuadraticMesh)


#endif /*QUADRATICMESH_HPP_*/
