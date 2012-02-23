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
#ifndef VERTEXMESH_HPP_
#define VERTEXMESH_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMeshWriter;

#include <iostream>
#include <map>
#include <algorithm>

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "AbstractMesh.hpp"
#include "ArchiveLocationInfo.hpp"
#include "VertexMeshReader.hpp"
#include "VertexMeshWriter.hpp"
#include "VertexElement.hpp"
#include "VertexElementMap.hpp"
#include "TetrahedralMesh.hpp"

/**
 * A vertex-based mesh class, for use in vertex-based simulations.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh : public AbstractMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestVertexMesh;

protected:

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> mElements;

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<ELEMENT_DIM-1, SPACE_DIM>*> mFaces;

    /**
     * Map that is used only when the vertex mesh is used to represent
     * a Voronoi tessellation, the dual to a Delaunay tetrahedral mesh.
     * The map consists of pairs (index1, index2), where index1 denotes
     * the global index of a node in the Delaunay mesh and index2 denotes
     * the global index of the corresponding element in the Voronoi mesh.
     */
    std::map<unsigned, unsigned> mVoronoiElementIndexMap;

    /**
     * Delaunay tetrahedral mesh that is used only when the vertex mesh
     * is used to represent a Voronoi tessellation. A pointer to the
     * Delaunay mesh is required in this case because the Delaunay mesh
     * may be a subclass of TetrahedralMesh, which overrides methods such as
     * GetVectorFromAtoB().
     */
    TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* mpDelaunayMesh;

    /**
     * Solve node mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the node
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Solve element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the element
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Solve boundary element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the boundary element
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /**
     * Populate mNodes with locations corresponding to the element
     * circumcentres of a given TetrahedralMesh. Used by 'Voronoi'
     * constructors.
     *
     * @param rMesh a tetrahedral mesh
     */
    void GenerateVerticesFromElementCircumcentres(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    //////////////////////////////////////////////////////////////////////
    //                        2D-specific methods                       //
    //////////////////////////////////////////////////////////////////////

    /**
     * Test whether a given point lies inside a given element.
     *
     * We use a ray-casting algorithm, which relies on the following result:
     * if the point in question is not on the boundary of the element, then
     * the number of intersections is an even number if the point is outside,
     * and it is odd if inside.
     *
     * Currently the method is coded 'strictly', such that points lying on
     * an edge or at a vertex are considered to lie outside the element.
     *
     * @param rTestPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return if the point is included in the element.
     */
    bool ElementIncludesPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex);

    /**
     * Get the local index of a given element which is the start vertex of the edge
     * of the element that the overlapping point rTestPoint is closest to.
     *
     * @param rTestPoint the point to test
     * @param elementIndex global index of the element in the mesh
     *
     * @return the local index
     */
    unsigned GetLocalIndexForElementEdgeClosestToPoint(const c_vector<double, SPACE_DIM>& rTestPoint, unsigned elementIndex);

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the VertexMesh and its member variables. Note that this will
     * write out a VertexMeshWriter file to wherever ArchiveLocationInfo has specified.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);

        // Create a mesh writer pointing to the correct file and directory
        VertexMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                                             ArchiveLocationInfo::GetMeshFilename(),
                                                             false);
        mesh_writer.WriteFilesUsingMesh(*(const_cast<VertexMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
    }

    /**
     * Loads a mesh by using VertexMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);

        VertexMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
        this->ConstructFromMeshReader(mesh_reader);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Forward declaration */
    class VertexElementIterator;

    /**
     * Get an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline VertexElementIterator GetElementIteratorBegin(bool skipDeletedElements=true);

    /**
     * Get an iterator to one past the last element in the mesh.
     */
    inline VertexElementIterator GetElementIteratorEnd();

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param vertexElements vector of pointers to VertexElements
     */
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements);

    /**
     * Constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param faces vector of pointer to VertexElements
     * @param vertexElements vector of pointers to VertexElement<3,3>s
     */
    VertexMesh(std::vector<Node<SPACE_DIM>*> nodes,
               std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*> faces,
               std::vector<VertexElement<ELEMENT_DIM,SPACE_DIM>*> vertexElements);

    /**
     * Alternative 2D 'Voronoi' constructor. Creates a Voronoi tessellation of a given tetrahedral mesh,
     * which must be Delaunay (see TetrahedralMesh::CheckIsVoronoi).
     *
     * \todo Merge with 3D Voronoi constructor? (#1075)
     *
     * @param rMesh a tetrahedral mesh
     * @param isPeriodic a boolean that indicates whether the mesh is periodic or not
     */
    VertexMesh(TetrahedralMesh<2,2>& rMesh, bool isPeriodic=false);

    /**
     * Alternative 3D 'Voronoi' constructor. Creates a Voronoi tessellation of a given tetrahedral mesh,
     * which must be Delaunay (see TetrahedralMesh::CheckIsVoronoi).
     *
     * \todo Merge with 2D Voronoi constructor? (#1075)
     *
     * @param rMesh a tetrahedral mesh
     */
    VertexMesh(TetrahedralMesh<3,3>& rMesh);

    /**
     * Default constructor for use by serializer.
     */
    VertexMesh();

    /**
     * Destructor.
     */
    virtual ~VertexMesh();

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of VertexElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of VertexElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements() const;

    /**
     * @return the number of Faces in the mesh.
     */
    virtual unsigned GetNumFaces() const;

    /**
     * @param index  the global index of a specified vertex element.
     *
     * @return a pointer to the vertex element
     */
    VertexElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * @param index  the global index of a specified face.
     *
     * @return a pointer to the face
     */
    VertexElement<ELEMENT_DIM-1, SPACE_DIM>* GetFace(unsigned index) const;

    /**
     * Compute the centroid of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x,centroid_y).
     */
    virtual c_vector<double, SPACE_DIM> GetCentroidOfElement(unsigned index);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM,SPACE_DIM>& rMeshReader);

    /**
     * Delete mNodes, mFaces and mElements.
     */
    virtual void Clear();

    /**
     * Given the global index of an element in the Voronoi mesh, returns the
     * global index of the corresponding element in the Delaunay mesh.
     *
     * @param elementIndex global index of an element in the Voronoi mesh
     */
    unsigned GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(unsigned elementIndex);

    /**
     * Given the global index of a node in the Delaunay mesh, returns the
     * global index of the corresponding element in the Voronoi mesh or
     * throws an exception if this does not exist.
     *
     * @param nodeIndex global index of a node in the Delaunay mesh
     */
    unsigned GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(unsigned nodeIndex);

    /**
     * Overridden GetVectorFromAtoB() method. Returns a vector between two points in space.
     *
     * If the mesh is being used to represent a Voronoi tessellation, and mpDelaunayMesh
     * is not NULL, then use that to compute GetVectorFromAtoB.
     *
     * @param rLocationA a c_vector of coordinates
     * @param rLocationB a c_vector of coordinates
     *
     * @return c_vector from location A to location B.
     */
    virtual c_vector<double, SPACE_DIM> GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocationA,
                                                          const c_vector<double, SPACE_DIM>& rLocationB);

    /**
     * Get the volume (or area in 2D, or length in 1D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the volume of the element
     */
    virtual double GetVolumeOfElement(unsigned index);

    /**
     * Compute the surface area (or perimeter in 2D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the surfacearea of the element
     */
    virtual double GetSurfaceAreaOfElement(unsigned index);

    //////////////////////////////////////////////////////////////////////
    //                        2D-specific methods                       //
    //////////////////////////////////////////////////////////////////////

    /**
     * Compute the area gradient of a 2D element at one of its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the area of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetAreaGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the edge of a 2D element ending at its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the edge of the element that ends at this node.
     */
    c_vector<double, SPACE_DIM> GetPreviousEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the edge of a 2D element starting at its nodes.
     *
     * N.B. This calls GetVectorFromAtoB(), which can be overridden
     * in daughter classes for non-Euclidean metrics.
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the edge of the element that starts at this node.
     */
    c_vector<double, SPACE_DIM> GetNextEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the perimeter of a 2D element at its nodes.
     * This returns the sum of GetPreviousEdgeGradientAtNode() and GetNextEdgeGradientAtNode().
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the perimeter of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetPerimeterGradientOfElementAtNode(VertexElement<ELEMENT_DIM,SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the second moments of area of a given 2D element.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (Ixx,Iyy,Ixy).
     */
    virtual c_vector<double, 3> CalculateMomentsOfElement(unsigned index);

    /**
     * Get the length of the edge separating two given elements in 2D.
     *
     * @param elementIndex1 index of an element in the mesh
     * @param elementIndex2 index of an element in the mesh
     */
    double GetEdgeLength(unsigned elementIndex1, unsigned elementIndex2);

    //////////////////////////////////////////////////////////////////////
    //                        3D-specific methods                       //
    //////////////////////////////////////////////////////////////////////

    /**
     * Compute the unit normal vector to a given face in 3D. This is achieved from a triangle
     * of vertices of the face. Note: this may return the outward or inward normal, depending
     * on your point of view.
     *
     * @param pFace a face in the mesh
     *
     * @return the unit normal
     */
    c_vector<double, SPACE_DIM> GetUnitNormalToFace(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace);

    /**
     * Get the area of a given face in 3D. This is achieved by projecting the face onto a 2D plane.
     * To avoid degeneracy and optimize robustness, we choose to ignore the dimension of the component
     * of the unit normal to the plane with the greatest absolute value.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param pFace a face in the mesh
     *
     * @return the area
     */
    virtual double GetAreaOfFace(VertexElement<ELEMENT_DIM-1, SPACE_DIM>* pFace);

    /**
     * Calculate the vector of the shortest axis of a given 2D element.
     * This is the eigenvector associated with the largest eigenvalue
     * of the inertial tensor. If the polygon is regular then the
     * eigenvalues are the same, so we return a random unit vector.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (short_axis_x, short_axis_y).
     */
    c_vector<double, SPACE_DIM> GetShortAxisOfElement(unsigned index);

    /**
     * Given a node, find a set containing the indices of its neighbouring nodes.
     *
     * @param nodeIndex global index of the node
     * @return its neighbouring node indices
     */
    std::set<unsigned> GetNeighbouringNodeIndices(unsigned nodeIndex);

    /**
     * Given a node and one of its containing elements, find a set containing
     * the indices of those neighbouring node(s) that are NOT also in the element.
     *
     * Note that we allow for more than one such index, since there is no reason
     * a priori to assume that each node is contained by exactly three elements.
     *
     * @param nodeIndex global index of the node
     * @param elemIndex global index of the element
     *
     * @return its neighbouring nodes that are not in the element
     */
    std::set<unsigned> GetNeighbouringNodeNotAlsoInElement(unsigned nodeIndex, unsigned elemIndex);

    /**
     * Given an element, find a set containing the indices of its neighbouring elements.
     *
     * @param elementIndex global index of the element
     * @return its neighbouring element indices
     */
    std::set<unsigned> GetNeighbouringElementIndices(unsigned elementIndex);

    //////////////////////////////////////////////////////////////////////
    //                         Nested classes                           //
    //////////////////////////////////////////////////////////////////////

    /**
     * A smart iterator over the elements in the mesh.
     *
     * \todo This is the same as in AbstractTetrahedralMesh and PottsMesh - merge? (#1379)
     */
    class VertexElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         *
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline VertexElement<ELEMENT_DIM, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         */
        inline VertexElement<ELEMENT_DIM, SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& rOther);

        /**
         * Prefix increment operator.
         */
        inline VertexElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * VertexMesh::GetElementIteratorBegin and VertexMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        VertexElementIterator(VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                        typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
                        bool skipDeletedElements=true);

    private:
        /** The mesh we're iterating over. */
        VertexMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM> *>::iterator mElementIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedElements;

        /**
         * Helper method to say when we're at the end.
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         */
        inline bool IsAllowedElement();
    };
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexMesh)


//////////////////////////////////////////////////////////////////////////////
// VertexElementIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorBegin(
        bool skipDeletedElements)
{
    return VertexElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorEnd()
{
    return VertexElementIterator(*this, mElements.end());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>& VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator!=(const VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    }
    while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::VertexElementIterator(
        VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
        bool skipDeletedElements)
    : mrMesh(rMesh),
      mElementIter(elementIter),
      mSkipDeletedElements(skipDeletedElements)
{
    if (mrMesh.mElements.empty())
    {
        // Cope with empty meshes
        mElementIter = mrMesh.mElements.end();
    }
    else
    {
        // Make sure we start at an allowed element
        if (mElementIter == mrMesh.mElements.begin() && !IsAllowedElement())
        {
            ++(*this);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}


#endif /*VERTEXMESH_HPP_*/
