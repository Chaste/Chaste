/*

Copyright (c) 2005-2012, University of Oxford.
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

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"

#include <vector>

/**
 * A concrete quadratic mesh class that inherits from TetrahedralMesh.
 */
template<unsigned DIM>
class QuadraticMesh : public TetrahedralMesh<DIM, DIM>
{
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
     * with that number of elements in each direction. This writes a temporary node file and uses
     * triangle to mesh this nodefile. The method overloads the equivalent method in
     * AbstractTetrahedralMesh. This is private, users should call ConstructRegularSlabMesh();
     *
     * @param numElemX Number of elements in x-direction (also, the width of the final mesh)
     * @param numElemY Number of elements in y-direction (also, the height of the final mesh)
     * @param unused (defaults to true; must always be true) is for compatibility of the
     *   interface of this method with same name in AbstractTetrahedralMesh.
     */
    void ConstructRectangularMesh(unsigned numElemX, unsigned numElemY, bool unused=true);

public:
    /**
     * Create a quadratic mesh on a rectangle from (0,0) to (numElemX,numElemY)
     * with that number of elements in each direction. 
     * @param numElemX Number of elements in x-direction (also, the width of the final mesh)
     * @param numElemY Number of elements in y-direction (also, the height of the final mesh)
     * @param unused (defaults to true; must always be true) is for compatibility of the
     *   interface of this method with same name in AbstractTetrahedralMesh.
     */
    void ConstructRectangularMeshNewImp(unsigned numElemX, unsigned numElemY, bool unused=true);
protected:
    /**
     * Create a quadratic mesh on a cuobid from (0,0,0) to (numElemX,numElemY,numElemZ)
     * with that number of elements in each direction. This writes a temporary node file and uses
     * tetgen to mesh this nodefile. The method overloads the equivalent method in
     * AbstractTetrahedralMesh. This is private, users should call ConstructRegularSlabMesh();
     *
     * @param numElemX Number of elements in x-direction (also, the width of the final mesh)
     * @param numElemY Number of elements in y-direction (also, the height of the final mesh)
     * @param numElemZ Number of elements in y-direction (also, the depth of the final mesh)
     */
    void ConstructCuboid(unsigned numElemX, unsigned numElemY, unsigned numElemZ);

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
     * with the given widths and given spacestep. In 1D height and depth need to passed in
     * as 0 (the default value), in 2D depth must be zero
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
     *  Get the number of vertices, ie non-internal (non-quadratic), nodes.
     */
    unsigned GetNumVertices() const;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(QuadraticMesh)


#endif /*QUADRATICMESH_HPP_*/
