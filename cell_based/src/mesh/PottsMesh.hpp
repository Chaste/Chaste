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

#ifndef POTTSMESH_HPP_
#define POTTSMESH_HPP_

// Forward declaration prevents circular include chain
template<unsigned DIM>
class PottsMeshWriter;

#include <iostream>
#include <map>
#include <algorithm>

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "AbstractMesh.hpp"
#include "ArchiveLocationInfo.hpp"
#include "PottsMeshReader.hpp"
#include "PottsMeshWriter.hpp"
#include "PottsElement.hpp"

/**
 * A Potts-based mesh class, for use in Cellular Potts model simulations.
 */
template<unsigned DIM>
class PottsMesh : public AbstractMesh<DIM, DIM>
{
    friend class TestPottsMesh;

protected:
    /** Vector of pointers to PottsElements. */
    std::vector<PottsElement<DIM>*> mElements;

    /**
     * Indices of elements that have been marked as deleted.
     * These indices can be reused when adding new elements.
     */
    std::vector<unsigned> mDeletedElementIndices;

    /** Vector of set of Von Neumann neighbours for each node. */
    std::vector< std::set<unsigned> > mVonNeumannNeighbouringNodeIndices;

    /** Vector of set of Moore neighbours for each node. */
    std::vector< std::set<unsigned> > mMooreNeighbouringNodeIndices;

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
     * Archive the PottsMesh and its member variables. Note that this will
     * write out a PottsMeshWriter file to wherever ArchiveLocationInfo has specified.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        // NOTE - Subclasses must archive their member variables BEFORE calling this method.
        archive & mDeletedElementIndices;
        archive & mVonNeumannNeighbouringNodeIndices;
        archive & mMooreNeighbouringNodeIndices;
        archive & boost::serialization::base_object<AbstractMesh<DIM, DIM> >(*this);

        // Create a mesh writer pointing to the correct file and directory
        PottsMeshWriter<DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                         ArchiveLocationInfo::GetMeshFilename(),
                                         false);
        mesh_writer.WriteFilesUsingMesh(*(const_cast<PottsMesh<DIM>*>(this)));
    }

    /**
     * Loads a mesh by using PottsMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        // NOTE - Subclasses must archive their member variables BEFORE calling this method.
        archive & mDeletedElementIndices;
        archive & mVonNeumannNeighbouringNodeIndices;
        archive & mMooreNeighbouringNodeIndices;
        archive & boost::serialization::base_object<AbstractMesh<DIM, DIM> >(*this);

        PottsMeshReader<DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
        this->ConstructFromMeshReader(mesh_reader);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Forward declaration */
    class PottsElementIterator;

    /**
     * @return an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline PottsElementIterator GetElementIteratorBegin(bool skipDeletedElements=true);

    /**
     * @return an iterator to one past the last element in the mesh.
     */
    inline PottsElementIterator GetElementIteratorEnd();

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param pottsElements vector of pointers to PottsElements
     * @param vonNeumannNeighbouringNodeIndices vector of set of Moore neighbours for each node
     * @param mooreNeighbouringNodeIndices vector of set of Von Neumann neighbours for each node
     */
    PottsMesh(std::vector<Node<DIM>*> nodes,
              std::vector<PottsElement<DIM>*> pottsElements,
              std::vector< std::set<unsigned> > vonNeumannNeighbouringNodeIndices,
              std::vector< std::set<unsigned> > mooreNeighbouringNodeIndices);

    /**
     * Default constructor for use by serializer.
     */
    PottsMesh();

    /**
     * Destructor.
     */
    virtual ~PottsMesh();

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of PottsElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of PottsElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements() const;

    /**
     * @param index  the global index of a specified PottsElement.
     *
     * @return a pointer to the PottsElement
     */
    PottsElement<DIM>* GetElement(unsigned index) const;

    /**
     * Compute the centroid of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x,centroid_y).
     */
    virtual c_vector<double, DIM> GetCentroidOfElement(unsigned index);

    /**
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    void ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader);

    /**
     * Delete mNodes and mElements.
     */
    virtual void Clear();

    /**
     * Get the volume (or area in 2D, or length in 1D) of a PottsElement.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified PottsElement element
     *
     * @return the volume of the element
     */
    virtual double GetVolumeOfElement(unsigned index);

    /**
     * Compute the surface area (or perimeter in 2D) of a PottsElement.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified PottsElement
     *
     * @return the surface area of the element
     */
    virtual double GetSurfaceAreaOfElement(unsigned index);

    /**
     * Given a node, return a set containing the indices of its Moore neighbouring nodes.
     *
     * @param nodeIndex global index of the node
     * @return neighbouring node indices in Moore neighbourhood
     */
    std::set<unsigned> GetMooreNeighbouringNodeIndices(unsigned nodeIndex);

    /**
     * Given a node, return a set containing the indices of its Von Neumann neighbouring nodes.
     *
     * @param nodeIndex global index of the node
     * @return neighbouring node indices in Von Neumann neighbourhood
     */
    std::set<unsigned> GetVonNeumannNeighbouringNodeIndices(unsigned nodeIndex);

    /**
     * Mark a node as deleted. Note that in a Potts mesh this requires the elements and connectivity to be updated accordingley.
     *
     * @param index  the global index of a specified node
     */
    void DeleteNode(unsigned index);

    /**
     * Mark an element as deleted. Note that in a Potts mesh this does not
     * delete the nodes so no remeshing is required.
     *
     * @param index  the global index of a specified Potts element
     */
    void DeleteElement(unsigned index);

    /**
     * Remove deleted elements and reorder them appropriately
     *
     */
    void RemoveDeletedElements();

    /**
     * Divide an element by assigning half the nodes to each new element in numerical order.
     * If an odd number of nodes then the existing element has one more node than the new element.
     *
     *
     * @param pElement the element to divide
     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
     *
     * @return the index of the new element
     */
    unsigned DivideElement(PottsElement<DIM>* pElement,
                           bool placeOriginalElementBelow=false);

    /**
     * Add an element to the mesh.
     *
     * @param pNewElement the new element
     *
     * @return the index of the new element in the mesh
     */
    unsigned AddElement(PottsElement<DIM>* pNewElement);

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
     * \todo This is the same as in AbstractTetrahedralMesh and VertexMesh- merge? (#1379)
     */
    class PottsElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         *
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         * @return reference
         */
        inline PottsElement<DIM>& operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline PottsElement<DIM>* operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         * @return not-equal
         */
        inline bool operator!=(const typename PottsMesh<DIM>::PottsElementIterator& rOther);

        /**
         * Prefix increment operator.
         * @return incremented object
         */
        inline PottsElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * PottsMesh::GetElementIteratorBegin and PottsMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements (defaults to true)
         */
        PottsElementIterator(PottsMesh<DIM>& rMesh,
                             typename std::vector<PottsElement<DIM> *>::iterator elementIter,
                             bool skipDeletedElements=true);

    private:
        /** The mesh we're iterating over. */
        PottsMesh<DIM>& mrMesh;

        /** The actual element iterator. */
        typename std::vector<PottsElement<DIM> *>::iterator mElementIter;

        /** @return whether to skip deleted elements. */
        bool mSkipDeletedElements;

        /**
         *  @return whether we're at the end.
         */
        inline bool IsAtEnd();

        /**
         *  @return whether we're allowed to point at this element.
         */
        inline bool IsAllowedElement();
    };
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PottsMesh)

//////////////////////////////////////////////////////////////////////////////
// PottsElementIterator class implementation - most methods are inlined     //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
typename PottsMesh<DIM>::PottsElementIterator PottsMesh<DIM>::GetElementIteratorBegin(
        bool skipDeletedElements)
{
    return PottsElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template<unsigned DIM>
typename PottsMesh<DIM>::PottsElementIterator PottsMesh<DIM>::GetElementIteratorEnd()
{
    return PottsElementIterator(*this, mElements.end());
}

template<unsigned DIM>
PottsElement<DIM>& PottsMesh<DIM>::PottsElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template<unsigned DIM>
PottsElement<DIM>* PottsMesh<DIM>::PottsElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template<unsigned DIM>
bool PottsMesh<DIM>::PottsElementIterator::operator!=(const typename PottsMesh<DIM>::PottsElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template<unsigned DIM>
typename PottsMesh<DIM>::PottsElementIterator& PottsMesh<DIM>::PottsElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    }
    while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template<unsigned DIM>
PottsMesh<DIM>::PottsElementIterator::PottsElementIterator(
        PottsMesh<DIM>& rMesh,
        typename std::vector<PottsElement<DIM>*>::iterator elementIter,
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

template<unsigned DIM>
bool PottsMesh<DIM>::PottsElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template<unsigned DIM>
bool PottsMesh<DIM>::PottsElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}

#endif /*POTTSMESH_HPP_*/
