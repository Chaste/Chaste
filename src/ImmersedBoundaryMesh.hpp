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
#ifndef IMMERSEDBOUNDARYMESH_HPP_
#define IMMERSEDBOUNDARYMESH_HPP_

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ImmersedBoundaryMeshWriter;

#include <iostream>
#include <map>
#include <algorithm>

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "AbstractMesh.hpp"
#include "ArchiveLocationInfo.hpp"
#include "ImmersedBoundaryMeshReader.hpp"
#include "ImmersedBoundaryMeshWriter.hpp"
#include "ImmersedBoundaryElement.hpp"
#include "ImmersedBoundaryArray.hpp"
#include "FluidSource.hpp"

/**
 * An immersed boundary mesh class, in which elements may contain different numbers of nodes.
 * This is facilitated by the ImmersedBoundaryElement class.
 *
 * This class allows immersed boundary simulations.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ImmersedBoundaryMesh : public AbstractMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestImmersedBoundaryMesh;

protected:

    /** Number of grid points in x direction */
    unsigned mNumGridPtsX;

    /** Number of grid points in y direction */
    unsigned mNumGridPtsY;

    /** Whether there is a membrane */
    bool mMeshHasMembrane;

    /** A pointer to the immersed boundary membrane */
    unsigned mMembraneIndex;

    /** Characteristic node spacing */
    double mCharacteristicNodeSpacing;

    /** 2D grid for fluid x velocity */
    multi_array<double, 3> m2dVelocityGrids;

    /** 3D grid for fluid x velocity */
    multi_array<double, 4> m3dVelocityGrids;

    /** Vector of pointers to ImmersedBoundaryElements. */
    std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>*> mElements;

    /** Vector of fluid sources related to elements. */
    std::vector<FluidSource<SPACE_DIM>*> mElementFluidSources;

    /** Vector of fluid sources used to balance those of the elements. */
    std::vector<FluidSource<SPACE_DIM>*> mBalancingFluidSources;

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

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the ImmersedBoundaryMesh and its member variables. Note that this will
     * write out an ImmersedBoundaryMeshWriter file to wherever ArchiveLocationInfo has specified.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);

        // Create a mesh writer pointing to the correct file and directory
        ImmersedBoundaryMeshWriter<ELEMENT_DIM, SPACE_DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                                             ArchiveLocationInfo::GetMeshFilename(),
                                                             false);
        mesh_writer.WriteFilesUsingMesh(*(const_cast<ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
    }

    /**
     * Load a mesh by using ImmersedBoundaryMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);

        ImmersedBoundaryMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
        this->ConstructFromMeshReader(mesh_reader);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:

    /** Forward declaration of element iterator. */
    class ImmersedBoundaryElementIterator;

    /**
     * @return an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline ImmersedBoundaryElementIterator GetElementIteratorBegin(bool skipDeletedElements=true);

    /**
     * @return an iterator to one past the last element in the mesh.
     */
    inline ImmersedBoundaryElementIterator GetElementIteratorEnd();

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param elements vector of pointers to ImmersedBoundaryElements
     * @param numGridPtsX the number of grid points in the x direction
     * @param numGridPtsY the number of grid points in the y direction
     * @param the index of the basement membrane element
     */
    ImmersedBoundaryMesh(std::vector<Node<SPACE_DIM>*> nodes,
                         std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>*> elements,
                         unsigned numGridPtsX=128,
                         unsigned numGridPtsY=128,
                         unsigned membraneIndex=UINT_MAX);

    /**
     * Default constructor for use by serializer.
     */
    ImmersedBoundaryMesh();

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundaryMesh();

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of ImmersedBoundaryElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of ImmersedBoundaryElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements() const;

    /**
     * @return the number of fluid mesh points in the x direction.
     */
    unsigned GetNumGridPtsX() const;

    /**
     * @return the number of fluid mesh points in the y direction.
     */
    unsigned GetNumGridPtsY() const;

    /**
     * @return the characteristic node spacing.
     */
    double GetCharacteristicNodeSpacing() const;

    /**
     * @return the spacing ratio: Characteristic Node Spacing / Fluid Grid Spacing
     */
    double GetSpacingRatio() const;

    /**
     * Overridden GetVectorFromAtoB() method.
     *
     * Evaluates the (surface) distance between two points in a periodic
     * geometry.
     *
     * @param rLocation1 the x and y co-ordinates of point 1
     * @param rLocation2 the x and y co-ordinates of point 2
     * @return the vector from location1 to location2
     */
    c_vector<double, SPACE_DIM> GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocation1, const c_vector<double, SPACE_DIM>& rLocation2);

    /**
     * Move the node with a particular index to a new point in space.
     *
     * @param nodeIndex the index of the node to be moved
     * @param point the new target location of the node
     */
    void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point);

    /**
     * @return reference to non-modifiable 2d fluid velocity grids.
     */
    const multi_array<double, 3>& rGet2dVelocityGrids() const;

    /**
     * @return reference to non-modifiable 3d fluid velocity grids.
     */
    const multi_array<double, 4>& rGet3dVelocityGrids() const;

    /**
     * @return reference to modifiable 2d fluid velocity grids.
     */
    multi_array<double, 3>& rGetModifiable2dVelocityGrids();

    /**
     * @return reference to modifiable 3d fluid velocity grids.
     */
    multi_array<double, 4>& rGetModifiable3dVelocityGrids();

    /**
     * @return reference to the vector of nodes
     */
    std::vector<Node<SPACE_DIM>*>& rGetNodes();

    /**
     * @param the new number of fluid mesh points in the x direction.
     */
    void SetNumGridPtsX(unsigned mesh_points_x);

    /**
     * @param the new number of fluid mesh points in the x direction.
     */
    void SetNumGridPtsY(unsigned mesh_points_y);

    /**
     * @param the new number of fluid mesh points in both directions.
     */
    void SetNumGridPtsXAndY(unsigned numGridPts);

    /**
     * @param the new characteristic node spacing.
     */
    void SetCharacteristicNodeSpacing(double node_spacing);

    /**
     * @param unsigned index of the membrane element
     */
    void SetMembraneIndex(unsigned membrane_index);

    /**
     * @return a pointer to the membrane element, NULL if no membrane
     */
    ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* GetMembraneElement();

    /**
     * @return the global index of the membrane element (UINT_MAX if no membrane element)
     */
    unsigned GetMembraneIndex();

    /**
     * @return reference to vector of element-associated fluid sources
     */
    std::vector<FluidSource<SPACE_DIM>*>& rGetElementFluidSources();

    /**
     * @return reference to vector of balancing fluid sources
     */
    std::vector<FluidSource<SPACE_DIM>*>& rGetBalancingFluidSources();

    /**
     * @param index  the global index of a specified immersed boundary element.
     *
     * @return a pointer to the immersed boundary element
     */
    ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

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
     * @param index  the global index of a specified immersed boundary element
     *
     * @return (centroid_x, centroid_y).
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
     * Get the volume (or area in 2D, or length in 1D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified immersed boundary element
     *
     * @return the volume of the element
     */
    virtual double GetVolumeOfElement(unsigned index);

    /**
     * Compute the surface area (or perimeter in 2D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified immersed boundary element
     *
     * @return the surface area of the element
     */
    virtual double GetSurfaceAreaOfElement(unsigned index);

    /**
     * Compute the average node spacing of an element.
     *
     * @param index  the global index of a specified immersed boundary element
     * @return the surface area of the element
     */
    double GetAverageNodeSpacingOfElement(unsigned index);

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
     * @param index  the global index of a specified immersed boundary element
     *
     * @return (Ixx,Iyy,Ixy).
     */
    virtual c_vector<double, 3> CalculateMomentsOfElement(unsigned index);

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
     * @return the tortuosity of the mesh
     */
    double GetTortuosityOfMesh();

    /**
     * Helper function for GetSkewnessOfElementMassDistributionAboutAxis.
     *
     * @param pNodeA pointer to one node in the comparison
     * @param pNodeB pointer to one node in the comparison
     * @return if the x-location of node A is < the x-location of node B
     */
    bool CompareNodesAlongX(Node<SPACE_DIM>* pNodeA, Node<SPACE_DIM>* pNodeB);

    /**
     * Calculate the skewness of the mass-distribution of a given element along a line perpendicular to a given axis,
     * as a measure of element asymmetry.
     *
     * This method is only valid in 2D.
     *
     * @param elemIndex the index of the element in this mesh
     * @param axis the axis perpendicular to which the distribution is found
     * @return the skewness of the mass distribution for the given element about the given axis
     */
    double GetSkewnessOfElementMassDistributionAboutAxis(unsigned elemIndex, c_vector<double, SPACE_DIM> axis);

    /**
     * Calculate the bounding box of an element specified by its index.
     *
     * @param index the index of the element
     * @return a Chaste cuboid representing the element bounding box
     */
    ChasteCuboid<SPACE_DIM> CalculateBoundingBoxOfElement(unsigned index);


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
     * @param index  the global index of a specified immersed boundary element
     *
     * @return a unit vector giving the direction of the short axis
     */
    c_vector<double, SPACE_DIM> GetShortAxisOfElement(unsigned index);


    /**
     * A smart iterator over the elements in the mesh.
     *
     * \todo This is the same as in AbstractTetrahedralMesh and PottsMesh - merge? (#1379)
     */
    class ImmersedBoundaryElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         * @return reference
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * @return true if not equal
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        inline ImmersedBoundaryElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * ImmersedBoundaryMesh::GetElementIteratorBegin and ImmersedBoundaryMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        ImmersedBoundaryElementIterator(ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                        typename std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
                        bool skipDeletedElements=true);

    private:
        /** The mesh we're iterating over. */
        ImmersedBoundaryMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM> *>::iterator mElementIter;

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
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ImmersedBoundaryMesh)


//////////////////////////////////////////////////////////////////////////////
// ImmersedBoundaryElementIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorBegin(
        bool skipDeletedElements)
{
    return ImmersedBoundaryElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorEnd()
{
    return ImmersedBoundaryElementIterator(*this, mElements.end());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM>* ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::operator!=(const typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator& ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    }
    while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::ImmersedBoundaryElementIterator(
        ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        typename std::vector<ImmersedBoundaryElement<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
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
bool ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ImmersedBoundaryMesh<ELEMENT_DIM, SPACE_DIM>::ImmersedBoundaryElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}

#endif /*IMMERSEDBOUNDARYMESH_HPP_*/
