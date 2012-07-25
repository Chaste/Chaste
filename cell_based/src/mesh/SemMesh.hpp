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

#ifndef SEMMESH_HPP_
#define SEMMESH_HPP_

//// Forward declaration prevents circular include chain
//template<unsigned DIM>
//class SemMeshWriter;

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
//#include "SemMeshReader.hpp"
//#include "SemMeshWriter.hpp"
#include "PottsElement.hpp"

/**
 * A Subcellular Element mesh class, for use in subcellular element model simulations.
 */
template<unsigned DIM>
class SemMesh : public AbstractMesh<DIM, DIM>
{
    friend class TestSemMesh;

protected:

    /** Vector of pointers to PottsElements. */
    std::vector<PottsElement<DIM>*> mElements;

    /**
     * Indices of elements that have been marked as deleted.
     * These indices can be reused when adding new elements.
     */
    std::vector<unsigned> mDeletedElementIndices;

    /**
     * Indices of elements that have been marked as deleted.
     * These indices can be reused when adding new elements.
     */
    std::vector<unsigned> mDeletedNodeIndices;

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

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the SemMesh and its member variables. Note that this will
     * write out a SemMeshWriter file to wherever ArchiveLocationInfo has specified.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        // NOTE - Subclasses must archive their member variables BEFORE calling this method.
        archive & mDeletedElementIndices;
        archive & mDeletedNodeIndices;
        archive & boost::serialization::base_object<AbstractMesh<DIM, DIM> >(*this);

//        // Create a mesh writer pointing to the correct file and directory
//        SemMeshWriter<DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
//                                         ArchiveLocationInfo::GetMeshFilename(),
//                                         false);
//        mesh_writer.WriteFilesUsingMesh(*(const_cast<SemMesh<DIM>*>(this)));
    }

    /**
     * Loads a mesh by using SemMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        // NOTE - Subclasses must archive their member variables BEFORE calling this method.
        archive & mDeletedElementIndices;
        archive & mDeletedNodeIndices;
        archive & boost::serialization::base_object<AbstractMesh<DIM, DIM> >(*this);

//        SemMeshReader<DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
//        this->ConstructFromMeshReader(mesh_reader);
    }
//    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Forward declaration */
    class SemElementIterator;

    /**
     * Get an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline SemElementIterator GetElementIteratorBegin(bool skipDeletedElements=true);

    /**
     * Get an iterator to one past the last element in the mesh.
     */
    inline SemElementIterator GetElementIteratorEnd();

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Default constructor.
     *
     * @param pottsElements vector of pointers to PottsElements
     */
    SemMesh( std::vector<Node<DIM>*> nodes,
            std::vector<PottsElement<DIM>*> pottsElements);

    /**
     * Default constructor for use by serializer.
     */
    SemMesh();

    /**
     * Destructor.
     */
    virtual ~SemMesh();

    /**
     * @return the number of undeleted Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;

    /**
     * @return the number of Nodes in the mesh including deleted nodes.
     */
    virtual unsigned GetNumAllNodes() const;

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
     * Construct the mesh using a MeshReader.
     *
     * @param rMeshReader the mesh reader
     */
    //void ConstructFromMeshReader(AbstractMeshReader<DIM, DIM>& rMeshReader);

    /**
     * Delete mNodes and mElements.
     */
    virtual void Clear();

    /**
     * Remesh the SemMesh. This cleans out elements and nodes that are marked as
     * deleted and tidys up the vectors mElements and mNodes.
     */
    void ReMesh();

    /**
     * Mark an element as deleted.
     *
     * @param index  the global index of a specified Potts element
     */
    void DeleteElement(unsigned index);
//
//    /**
//     * Divide an element by assigning half the nodes to each new element in numerical order.
//     * If an odd number of nodes then the existing element has one more node than the new element.
//     *
//     *
//     * @param pElement the element to divide
//     * @param placeOriginalElementBelow whether to place the original element below (in the y direction) the new element (defaults to false)
//     *
//     * @return the index of the new element
//     */
//    unsigned DivideElement(PottsElement<DIM>* pElement,
//                           bool placeOriginalElementBelow=false);

    /**
     * Add an element to the mesh.
     *
     * @param pNewElement the new element
     * @param pNewNodes the new nodes associated with the element.
     * @return the index of the new element in the mesh
     */
    unsigned AddElement(PottsElement<DIM>* pNewElement, std::vector<Node<DIM>* > pNewNodes);

    //////////////////////////////////////////////////////////////////////
    //                         Nested classes                           //
    //////////////////////////////////////////////////////////////////////

    /**
     * A smart iterator over the elements in the mesh.
     *
     * \todo This is the same as in AbstractTetrahedralMesh and VertexMesh- merge? (#1379)
     */
    class SemElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         *
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline PottsElement<DIM>& operator*();

        /**
         * Member access from a pointer.
         */
        inline PottsElement<DIM>* operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const SemMesh<DIM>::SemElementIterator& rOther);

        /**
         * Prefix increment operator.
         */
        inline SemElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * SemMesh::GetElementIteratorBegin and SemMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements (defaults to true)
         */
        SemElementIterator(SemMesh<DIM>& rMesh,
                             typename std::vector<PottsElement<DIM> *>::iterator elementIter,
                             bool skipDeletedElements=true);

    private:
        /** The mesh we're iterating over. */
        SemMesh<DIM>& mrMesh;

        /** The actual element iterator. */
        typename std::vector<PottsElement<DIM> *>::iterator mElementIter;

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

//#include "SerializationExportWrapper.hpp"
//EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemMesh)

//////////////////////////////////////////////////////////////////////////////
// SemElementIterator class implementation - most methods are inlined     //
//////////////////////////////////////////////////////////////////////////////

template<unsigned DIM>
typename SemMesh<DIM>::SemElementIterator SemMesh<DIM>::GetElementIteratorBegin(
        bool skipDeletedElements)
{
    return SemElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template<unsigned DIM>
typename SemMesh<DIM>::SemElementIterator SemMesh<DIM>::GetElementIteratorEnd()
{
    return SemElementIterator(*this, mElements.end());
}
template<unsigned DIM>
PottsElement<DIM>& SemMesh<DIM>::SemElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template<unsigned DIM>
PottsElement<DIM>* SemMesh<DIM>::SemElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template<unsigned DIM>
bool SemMesh<DIM>::SemElementIterator::operator!=(const SemMesh<DIM>::SemElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template<unsigned DIM>
typename SemMesh<DIM>::SemElementIterator& SemMesh<DIM>::SemElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    }
    while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template<unsigned DIM>
SemMesh<DIM>::SemElementIterator::SemElementIterator(
        SemMesh<DIM>& rMesh,
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
bool SemMesh<DIM>::SemElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template<unsigned DIM>
bool SemMesh<DIM>::SemElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}

#endif /*SEMMESH_HPP_*/
