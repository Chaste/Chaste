/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
     * Check that internals appear after vertices in the list
     */
    void CountAndCheckVertices();

    /**
     * Top level method for making 2D edges have 3 nodes not 2 and making 3D faces have 6 nodes not 3  (ie linear to quadratic).
     * @param pMeshReader Pointer to the reader. Only used if boundaryElemFileHasContainElementInfo==true (can be null if not).
     */
    void AddNodesToBoundaryElements(TrianglesMeshReader<DIM,DIM>* pMeshReader);

    /**
     * This method adds the given node (defined by an element and a node index)
     * to the given boundary element, and also sets the node as a boundary
     * element and adds it to the std::vector of boundary elements.
     *
     * @param pBoundaryElement  pointer to a boundary element in the mesh
     * @param pElement  pointer to an element in the mesh
     * @param internalNode  index of a node in the mesh
     */
    void AddNodeToBoundaryElement(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                                  Element<DIM,DIM>* pElement,
                                  unsigned internalNode);

    /**
     * Given a face in an element (defined by giving an element and the opposite
     * node number to the face) that corresponds to a given boundary element,
     * this method adds in the face's internal nodes to the boundary element
     * (in the correct order).
     *
     * @param pBoundaryElement  pointer to a boundary element in the mesh
     * @param pElement  pointer to an element in the mesh
     * @param nodeIndexOppositeToFace  index of a node in the mesh
     */
    void AddExtraBoundaryNodes(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                               Element<DIM,DIM>* pElement,
                               unsigned nodeIndexOppositeToFace);

    /**
     * Nasty helper method for AddNodeToBoundaryElement() in 3D.
     *
     * This method takes in the three vertices of a face which match the given boundary
     * element, and figure out if the order of the nodes in the face is reversed in
     * the boundary element (returned in the bool 'rReverse'). Also, the offset between
     * the first node in the face (as given to this method) and the first node in
     * the boundary element is computed (returned in the variable 'rOffset'). Offset
     * should then be applied before reverse to match the face nodes to the boundary
     * element nodes.
     *
     * \todo document these parameters
     *
     * @param boundaryElemNode0
     * @param boundaryElemNode1
     * @param pElement
     * @param node0
     * @param node1
     * @param node2
     * @param rOffset
     * @param rReverse
     */
    void HelperMethod1(unsigned boundaryElemNode0, unsigned boundaryElemNode1,
                       Element<DIM,DIM>* pElement,
                       unsigned node0, unsigned node1, unsigned node2,
                       unsigned& rOffset,
                       bool& rReverse);

    /**
     * Nasty helper method for AddNodeToBoundaryElement() in 3D.
     *
     * This method takes the three internal nodes for some face in some element,
     * applies the given offset and reverse (see HelperMethod1) to them, to get
     * the ordered internal nodes which should given to the boundary element.
     * It then calls AddNodeToBoundaryElement with each of the three internal nodes.
     *
     * \todo document these parameters
     *
     * @param pBoundaryElement
     * @param pElement
     * @param internalNode0
     * @param internalNode1
     * @param internalNode2
     * @param offset
     * @param reverse
     */
    void HelperMethod2(BoundaryElement<DIM-1,DIM>* pBoundaryElement,
                       Element<DIM,DIM>* pElement,
                       unsigned internalNode0, unsigned internalNode1, unsigned internalNode2,
                       unsigned offset,
                       bool reverse);

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
     *  Write the boundary elements to file (in case the boundary elements were linear when read and the
     *  quadratic versions have been computed.
     *
     *  @param directory Directory relative to CHASTE_TEST_OUTPUT. Not cleaned
     *  @param fileName Boundary element file name.
     */
    void WriteBoundaryElementFile(std::string directory, std::string fileName);

    /**
     *  Get the number of vertices, ie non-internal (non-quadratic), nodes.
     */
    unsigned GetNumVertices();


//// These methods are required for the adaptive subclass of QuadraticMesh (projects work) to run.
//// In that class the vertices are not assumed to be the first num_vertices nodes in this->mNodes.
//// Therefore need a map from 'vertex index' to 'node index'. Used in the solvers
//// Probably want to go to this design in the future.
//    virtual unsigned GetVertexIndexOfNode(unsigned nodeIndex)
//    {
//        assert(nodeIndex < GetNumVertices());
//        return nodeIndex;
//    }
//
//    virtual unsigned GetNodeIndexOfVertex(unsigned vertexIndex)
//    {
//        assert(vertexIndex < GetNumVertices());
//        return vertexIndex;
//    }

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(QuadraticMesh)


#endif /*QUADRATICMESH_HPP_*/
