/*

Copyright (c) 2005-2023, University of Oxford.
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

// Forward declaration prevents circular include chain
template<unsigned DIM>
class SemMeshWriter;

#include <algorithm>
#include <iostream>
#include <map>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp>
#include "ChasteSerialization.hpp"

#include "AbstractMesh.hpp"
#include "ArchiveLocationInfo.hpp"
#include "SemElement.hpp"
#include "Element.hpp"
#include "SemMeshReader.hpp"
#include "SemMeshWriter.hpp"
#include "NodeMap.hpp"

/**
 * \todo Document class
 */
template<unsigned DIM>
class SemMesh : public AbstractMesh<DIM, DIM>
{
private:
    /** Vector of pointers to SemElements. */
    std::vector<SemElement<DIM>*> mElements;

    /**
     * Solve node mapping method. This overridden method is required as it is 
     * pure virtual in the base class.
     *
     * @param index the global index of the node
     * 
     * @return local index
     */
    unsigned SolveNodeMapping(unsigned index) const;

    /**
     * Solve element mapping method. This overridden method is required as it is 
     * pure virtual in the base class.
     *
     * @param index the global index of the element
     * 
     * @return local index
     */
    unsigned SolveElementMapping(unsigned index) const;

    /**
     * Solve boundary element mapping method. This overridden method is required
     * as it is pure virtual in the base class.
     *
     * @param index the global index of the boundary element
     * 
     * @return local index
     */
    unsigned SolveBoundaryElementMapping(unsigned index) const;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the SemMesh and its member variables. Note that this will write 
     * out a SemMeshWriter file to wherever ArchiveLocationInfo has specified.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void save(Archive& archive, const unsigned int version) const
    {
        archive& boost::serialization::base_object<AbstractMesh<DIM, DIM> >(*this);

        // Create a mesh writer pointing to the correct file and directory
        SemMeshWriter<DIM> writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                  ArchiveLocationInfo::GetMeshFilename(),
                                  false);
        writer.WriteFilesUsingMesh(*(const_cast<SemMesh<DIM>*>(this)));
    }

    /**
     * Load a mesh using SemMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void load(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractMesh<DIM, DIM> >(*this);

        SemMeshReader<DIM> reader(ArchiveLocationInfo::GetArchiveDirectory() 
                                  + ArchiveLocationInfo::GetMeshFilename());
        this->ConstructFromMeshReader(reader);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
    /** Forward declaration of SemElement iterator. */
    class SemElementIterator;

    /**
     * @return an iterator to the first SemElement in the SemMesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline SemElementIterator GetElementIteratorBegin(bool skipDeletedElements = true);

    /**
     * @return an iterator to one past the last SemElement in the SemMesh.
     */
    inline SemElementIterator GetElementIteratorEnd();

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param semElements vector of pointers to SemElements
     */
    SemMesh(std::vector<Node<DIM>*> nodes,
            std::vector<SemElement<DIM>*> semElements);

    /**
     * Default constructor for use by serializer.
     */
    SemMesh();
    

    /**
     * Destructor.
     */
    virtual ~SemMesh();

    /**
     * @return the number of Nodes in the mesh.
     */
    virtual unsigned GetNumNodes() const;
    
    /**
     * Adds a node to the mesh
     * 
     * @param pNewNode a pointer to the node to add
    */
   virtual unsigned AddNode(Node<DIM>* pNewNode);

    /**
     * @return the number of SemElements in the mesh.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of SemElements in the mesh, including those marked as deleted.
     */
    unsigned GetNumAllElements() const;

    /**
     * @param index the global index of a specified SemElement
     *
     * @return a pointer to the SemElement.
     */
    virtual SemElement<DIM>* GetElement(unsigned index) const;

    /**
     * Add a node to the mesh
     * 
     * @param pNewNode a pointer to the node to add to the mesh
     * 
     * @return the id of the node within the mesh
    */
    virtual unsigned AddElement(SemElement<DIM>* pNewElement);
    
    /**
     * Compute the centroid of an element.
     *
     * This must be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index the global index of a specified SemElement
     *
     * @return centroid of the element as a c_vector
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
     * Get the volume (or area in 2D, or length in 1D) of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified SemElement
     *
     * @return the volume of the element
     */
    virtual double GetVolumeOfElement(unsigned index);

    /**
     * Return a pointer to the SemMesh.
     *
     * This method may be overridden in daughter classes for non-Euclidean metrics.
     * This can then be used when writing to VTK.
     *
     * @return a pointer to the vertex mesh
     */
    virtual SemMesh<DIM>* GetMeshForVtk();
    
    void DeleteNodePriorToReMesh(unsigned int node);
    void ReMesh(NodeMap map);

    /**
     * A smart iterator over the SemElements in the SemMesh.
     */
    class SemElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current 
         * element. Make sure to use a reference for the result to avoid 
         * copying elements unnecessarily.
         * 
         * @return reference to the current SemElement.
         */
        inline SemElement<DIM>& operator*();

        /**
         * Member access from a pointer.
         * 
         * @return pointer to the current SemElement
         */
        inline SemElement<DIM>* operator->();

        /**
         * Comparison not-equal-to.
         * 
         * @param rOther iterator with which comparison is made
         * 
         * @return true if not equal.
         */
        inline bool operator!=(const typename SemMesh<DIM>::SemElementIterator& rOther);

        /**
         * Prefix increment operator.
         * 
         * @return reference to incremented object.
         */
        inline SemElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * SemMesh::GetElementIteratorBegin() and 
         * SemMesh::GetElementIteratorEnd() instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        SemElementIterator(SemMesh<DIM>& rMesh,
                              typename std::vector<SemElement<DIM>*>::iterator elementIter,
                              bool skipDeletedElements = true);

    private:
        /** The mesh we're iterating over. */
        SemMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<SemElement<DIM>*>::iterator mElementIter;

        /** Whether to skip deleted elements. */
        bool mSkipDeletedElements;

        /**
         * Helper method to say when we're at the end.
         * 
         * @return true if at end.
         */
        inline bool IsAtEnd();

        /**
         * Helper method to say if we're allowed to point at this element.
         * 
         * @return true if allowed.
         */
        inline bool IsAllowedElement();
    };
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemMesh)

// SemElementIterator class implementation - most methods are inlined

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
SemElement<DIM>& SemMesh<DIM>::SemElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template<unsigned DIM>
SemElement<DIM>* SemMesh<DIM>::SemElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template<unsigned DIM>
bool SemMesh<DIM>::SemElementIterator::operator!=(const typename SemMesh<DIM>::SemElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template<unsigned DIM>
typename SemMesh<DIM>::SemElementIterator& SemMesh<DIM>::SemElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    } while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template<unsigned DIM>
SemMesh<DIM>::SemElementIterator::SemElementIterator(
    SemMesh<DIM>& rMesh,
    typename std::vector<SemElement<DIM>*>::iterator elementIter,
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