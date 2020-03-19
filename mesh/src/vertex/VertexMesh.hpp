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
#ifndef VERTEXMESH_HPP_
#define VERTEXMESH_HPP_

// Forward declaration prevents circular include chain
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMeshWriter;

#include <algorithm>
#include <iostream>
#include <map>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractMesh.hpp"
#include "ArchiveLocationInfo.hpp"
#include "TetrahedralMesh.hpp"
#include "VertexElement.hpp"
#include "VertexElementMap.hpp"
#include "VertexMeshReader.hpp"
#include "VertexMeshWriter.hpp"


/**
 * A vertex-based mesh class, in which elements may contain different numbers of nodes.
 * This is facilitated by the VertexElement class.
 *
 * This class has two applications in the cell_based code.
 *
 * First, VertexMesh is used as a member of the MeshBasedCellPopulation class to represent
 * a Voronoi tessellation, the dual to a Delaunay mesh, which allows the shapes of cells
 * to be visualised in simulations of a class of off-lattice cell centre-based models.
 *
 * Second, VertexMesh serves as a parent class for MutableVertexMesh, which is used as a
 * member of the VertexBasedCellPopulation class to represent the junctional network of
 * cells that forms the basis of simulations of off-lattice vertex-based models.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class VertexMesh : public AbstractMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestVertexMesh;


protected:
    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> mElements;

    /** Vector of pointers to VertexElements. */
    std::vector<VertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> mFaces;

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
     * @return local index
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Solve element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the element
     * @return local index
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Solve boundary element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the boundary element
     * @return local index
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /**
     * Build edges from elements. Populates edges in EdgeHelper class
     * @param elements from which edges are built
     */
    void GenerateEdgesFromElements(std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> &elements);

    /**
     * Populate mNodes with locations corresponding to the element
     * circumcentres of a given TetrahedralMesh. Used by 'Voronoi'
     * constructors.
     *
     * @param rMesh a tetrahedral mesh
     */
    void GenerateVerticesFromElementCircumcentres(TetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh);

    /**
     * Test whether a given point lies inside a given element.
     *
     * We use a winding number test, which counts the number of times the
     * polygon associated with the element winds around the given point.
     * The point is outside only when this "winding number" vanishes;
     * otherwise, the point is inside.
     *
     * One must decide whether a point on the polygon's boundary is inside
     * or outside: we adopt the standard convention that a point on a left
     * or bottom edge is inside, and a point on a right or top edge is outside.
     * This way, if two distinct polygons share a common boundary segment,
     * then a point on that segment will be in one polygon or the other, but
     * not both at the same time.
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
    template <class Archive>
    void save(Archive& archive, const unsigned int version) const
    {
        archive& boost::serialization::base_object<AbstractMesh<ELEMENT_DIM, SPACE_DIM> >(*this);
        // Create a mesh writer pointing to the correct file and directory
        VertexMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                                             ArchiveLocationInfo::GetMeshFilename(),
                                                             false);
        mesh_writer.WriteFilesUsingMesh(*(const_cast<VertexMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
    }

    /**
     * Load a mesh by using VertexMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void load(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractMesh<ELEMENT_DIM, SPACE_DIM> >(*this);

        VertexMeshReader<ELEMENT_DIM, SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
        this->ConstructFromMeshReader(mesh_reader);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
    /** Forward declaration of element iterator. */
    class VertexElementIterator;

    /**
     * @return an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline VertexElementIterator GetElementIteratorBegin(bool skipDeletedElements = true);

    /**
     * @return an iterator to one past the last element in the mesh.
     */
    inline VertexElementIterator GetElementIteratorEnd();

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
               std::vector<VertexElement<ELEMENT_DIM - 1, SPACE_DIM>*> faces,
               std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*> vertexElements);

    /**
     * @brief  Alternative 2D 'Voronoi' constructor.
     *
     * This VertexMesh constructor is currently only defined for 2D meshes.
     *
     * Creates a Voronoi tessellation of a given tetrahedral mesh,
     * which must be Delaunay (see TetrahedralMesh::CheckIsVoronoi).
     *
     * \todo Merge with 3D Voronoi constructor? (see #1075)
     *
     * @param rMesh a tetrahedral mesh
     * @param isPeriodic a boolean that indicates whether the mesh is periodic or not
     */
    VertexMesh(TetrahedralMesh<2, 2>& rMesh, bool isPeriodic = false);

    /**
     * Alternative 3D 'Voronoi' constructor. Creates a Voronoi tessellation of a given tetrahedral mesh,
     * which must be Delaunay (see TetrahedralMesh::CheckIsVoronoi).
     *
     * \todo Merge with 2D Voronoi constructor? (see #1075)
     *
     * @param rMesh a tetrahedral mesh
     */
    VertexMesh(TetrahedralMesh<3, 3>& rMesh);

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
    VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* GetFace(unsigned index) const;

    /**
     * Compute the centroid of an element.
     *
     * A formula for the centroid of a plane polygon may be found e.g. in the following reference:
     *
     * Mechanics of Materials
     * James M. Gere (Author), Barry J. Goodno.
     * Cengage Learning; 8th edition (January 1, 2012)
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x, centroid_y).
     */
    virtual c_vector<double, SPACE_DIM> GetCentroidOfElement(unsigned index);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader);

    /**
     * Delete mNodes, mFaces and mElements.
     */
    virtual void Clear();

    /**
     * @return the global index of the corresponding element in the Delaunay mesh,
     * given the global index of an element in the Voronoi mesh.
     * @param elementIndex global index of an element in the Voronoi mesh
     */
    unsigned GetDelaunayNodeIndexCorrespondingToVoronoiElementIndex(unsigned elementIndex);

    /**
     * @return the global index of the corresponding element in the Voronoi
     * mesh given the global index of a node in the Delaunay mesh,  or
     * throws an exception if this does not exist.
     *
     * @param nodeIndex global index of a node in the Delaunay mesh
     */
    unsigned GetVoronoiElementIndexCorrespondingToDelaunayNodeIndex(unsigned nodeIndex);

    /**
     * Get the "rosette rank" of an element.
     *
     * This is defined as the maximum number of elements shared by any node in the specified element.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return the rosette rank of the element
     */
    unsigned GetRosetteRankOfElement(unsigned index);

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
    c_vector<double, SPACE_DIM> GetAreaGradientOfElementAtNode(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, unsigned localIndex);

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
    c_vector<double, SPACE_DIM> GetPreviousEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, unsigned localIndex);

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
    c_vector<double, SPACE_DIM> GetNextEdgeGradientOfElementAtNode(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the gradient of the perimeter of a 2D element at its nodes.
     * This returns the sum of GetPreviousEdgeGradientAtNode() and GetNextEdgeGradientAtNode().
     *
     * @param pElement  pointer to a specified vertex element
     * @param localIndex  local index of a node in this element
     *
     * @return the gradient of the perimeter of the element, evaluated at this node.
     */
    c_vector<double, SPACE_DIM> GetPerimeterGradientOfElementAtNode(VertexElement<ELEMENT_DIM, SPACE_DIM>* pElement, unsigned localIndex);

    /**
     * Compute the second moments and product moment of area for a given 2D element
     * about its centroid. These are:
     *
     * I_xx, the second moment of area about an axis through the centroid of the
     * element parallel to the x-axis;
     *
     * I_yy, the second moment of area about an axis through the centroid of the
     * element parallel to the y-axis;
     *
     * and I_xy, product moment of area through the centroid of the element.
     *
     * Formulae for these quantities may be found e.g. in the following reference:
     *
     * Mechanics of Materials
     * James M. Gere (Author), Barry J. Goodno.
     * Cengage Learning; 8th edition (January 1, 2012)
     *
     * This method is used within GetShortAxisOfElement() to compute the direction
     * of the shortest principal axis passing through the centroid, or 'short axis',
     * of the element.
     *
     * Note that by definition, the second moments of area must be non-negative,
     * while the product moment of area may not be.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (Ixx,Iyy,Ixy).
     */
    virtual c_vector<double, 3> CalculateMomentsOfElement(unsigned index);

    /**
     * @return the length of the edge separating two given elements in 2D.
     *
     * @param elementIndex1 index of an element in the mesh
     * @param elementIndex2 index of an element in the mesh
     */
    double GetEdgeLength(unsigned elementIndex1, unsigned elementIndex2);

    /**
     * Get the elongation shape factor of a given element.
     * This is defined as the square root of the ratio of
     * the two second moments of the element around its
     * principal axes.
     *
     * @param elementIndex index of an element in the mesh
     *
     * @return the elongation shape factor of the element.
     */
    double GetElongationShapeFactorOfElement(unsigned elementIndex);

    /**
     * Compute the unit normal vector to a given face in 3D. This is achieved by calculating scaled normal,
     * which is the effective sum of signed areas of triangle forming the face.
     * Note: this may return the outward or inward normal, depending
     * on the face chirality.
     *
     * @param pFace a face in the mesh
     * @param rNormal vector in which to return the unit normal
     *
     * @return the area
     */
    double CalculateUnitNormalToFaceWithArea(VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace, c_vector<double, SPACE_DIM>& rNormal);

    /**
     * Get the area of a given face in 3D.  Uses CalculateUnitNormalToFaceWithArea
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param pFace a face in the mesh
     *
     * @return the area
     */
    virtual double CalculateAreaOfFace(VertexElement<ELEMENT_DIM - 1, SPACE_DIM>* pFace);

    /**
     * Compute the direction of the shortest principal axis passing through the centroid,
     * or 'short axis', of a given element. This is the eigenvector associated with the
     * eigenvalue of largest magnitude of the inertia matrix
     *
     * J = (  I_xx  -I_xy )
     *     ( -I_xy   I_yy )
     *
     * whose entries are computed by calling the method CalculateMomentsOfElement().
     *
     * Note that if the nodes owned by the element are supplied in clockwise rather than
     * anticlockwise manner, or if this arises when any periodicity is enforced, then the
     * sign of each moment may be incorrect change. This means that we need to consider the eigenvalue
     * of largest magnitude rather than largest value when computing the short axis of the
     * element.
     *
     * If the element is a regular polygon then the eigenvalues of the inertia tensor are
     * equal: in this case we return a random unit vector.
     *
     * This method is only implemented in 2D at present.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return a unit vector giving the direction of the short axis
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

    /**
     * Return a pointer to the vertex mesh
     *
     * This method may be overridden in daughter classes for non-Euclidean metrics.
     * This can then be used when writing to VTK.
     *
     * @return a pointer to the vertex mesh
     */
    virtual VertexMesh<ELEMENT_DIM, SPACE_DIM>* GetMeshForVtk();

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
         * @return reference
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline VertexElement<ELEMENT_DIM, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline VertexElement<ELEMENT_DIM, SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * @return true if not equal
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
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
                              typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*>::iterator elementIter,
                              bool skipDeletedElements = true);

    private:
        /** The mesh we're iterating over. */
        VertexMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*>::iterator mElementIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedElements;

        /**
         * Helper method to say when we're at the end.
         * @return true if at end
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         * @return true if allowed
         */
        inline bool IsAllowedElement();
    };
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VertexMesh)

//////////////////////////////////////////////////////////////////////////////
// VertexElementIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorBegin(
    bool skipDeletedElements)
{
    return VertexElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator VertexMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorEnd()
{
    return VertexElementIterator(*this, mElements.end());
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>& VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>* VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator!=(const typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator& VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    } while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::VertexElementIterator(
    VertexMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
    typename std::vector<VertexElement<ELEMENT_DIM, SPACE_DIM>*>::iterator elementIter,
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

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexMesh<ELEMENT_DIM, SPACE_DIM>::VertexElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}

#endif /*VERTEXMESH_HPP_*/
