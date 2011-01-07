/*

Copyright (C) University of Oxford, 2005-2011

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
#ifndef POTTSMESH_HPP_
#define POTTSMESH_HPP_

// Forward declaration prevents circular include chain
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
//class PottsMeshWriter;

#include <iostream>
#include <map>
#include <algorithm>

#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include "AbstractMesh.hpp"
#include "PottsElement.hpp"

/**
 * A potts-based mesh class, for use in potts-based simulations.
 */
//template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class PottsMesh : public AbstractMesh<2, 2>
{
    friend class TestPottsMesh;

protected:
    /** Vector of pointers to PottsElements. */
    std::vector<PottsElement *> mElements;

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

//    /**
//     * Archive the PottsMesh and its member variables. Note that this will
//     * write out a PottsMeshWriter file to wherever ArchiveLocationInfo has specified.
//     *
//     * @param archive the archive
//     * @param version the current version of this class
//     */
//    template<class Archive>
//    void save(Archive & archive, const unsigned int version) const
//    {
//        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);
//
//        // Create a mesh writer pointing to the correct file and directory
//        PottsMeshWriter mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
//                                                             ArchiveLocationInfo::GetMeshFilename(),
//                                                             false);
//        mesh_writer.WriteFilesUsingMesh(*(const_cast<PottsMesh*>(this)));
//    }
//
//    /**
//     * Loads a mesh by using PottsMeshReader and the location in ArchiveLocationInfo.
//     *
//     * @param archive the archive
//     * @param version the current version of this class
//     */
//    template<class Archive>
//    void load(Archive & archive, const unsigned int version)
//    {
//        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);
//
//        PottsMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());
//        this->ConstructFromMeshReader(mesh_reader);
//    }
//    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Forward declaration */
    class PottsElementIterator;

    /**
     * Get an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline PottsElementIterator GetElementIteratorBegin(bool skipDeletedElements=true);

    /**
     * Get an iterator to one past the last element in the mesh.
     */
    inline PottsElementIterator GetElementIteratorEnd();

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Default constructor.
     *
     * @param nodes vector of pointers to nodes
     * @param PottsElements vector of pointers to PottsElements
     */
    PottsMesh(std::vector<Node<2>*> nodes,
               std::vector<PottsElement*> PottsElements);

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
     * @param index  the global index of a specified vertex element.
     *
     * @return a pointer to the vertex element
     */
    PottsElement* GetElement(unsigned index) const;

    /**
     * Compute the centroid of an element.
     *
     * This needs to be overridden in daughter classes for non-Euclidean metrics.
     *
     * @param index  the global index of a specified vertex element
     *
     * @return (centroid_x,centroid_y).
     */
    virtual c_vector<double, 2> GetCentroidOfElement(unsigned index);

//    /*
//     * Construct the mesh using a MeshReader.
//     *
//     * @param rMeshReader the mesh reader
//     */
//    void ConstructFromMeshReader(AbstractMeshReader<2,2>& rMeshReader);

    /**
     * Delete mNodes, mFaces and mElements.
     */
    virtual void Clear();

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
    virtual c_vector<double, 2> GetVectorFromAtoB(const c_vector<double, 2>& rLocationA,
                                                          const c_vector<double, 2>& rLocationB);

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
    //                         Nested classes                           //
    //////////////////////////////////////////////////////////////////////

    /**
     * A smart iterator over the elements in the mesh.
     *
     * \todo This is the same as in AbstractTetrahedralMesh - merge? (#1379)
     */
    class PottsElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         *
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline PottsElement& operator*();

        /**
         * Member access from a pointer.
         */
        inline PottsElement* operator->();

        /**
         * Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const PottsMesh::PottsElementIterator& rOther);

        /**
         * Prefix increment operator.
         */
        inline PottsElementIterator& operator++();

        /**
         * Constructor for a new iterator. NOTE this is moved here as this class is currently not templated.
         *
         * This should not be called directly by user code; use the mesh methods
         * PottsMesh::GetElementIteratorBegin and PottsMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        PottsElementIterator(PottsMesh& rMesh,
                             std::vector<PottsElement *>::iterator elementIter,
                             bool skipDeletedElements=true)
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



    private:
        /** The mesh we're iterating over. */
        PottsMesh& mrMesh;

        /** The actual element iterator. */
        std::vector<PottsElement *>::iterator mElementIter;

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
//EXPORT_TEMPLATE_CLASS_ALL_DIMS(PottsMesh)


//////////////////////////////////////////////////////////////////////////////
// PottsElementIterator class implementation - most methods are inlined    //
//////////////////////////////////////////////////////////////////////////////

PottsMesh::PottsElementIterator PottsMesh::GetElementIteratorBegin(
        bool skipDeletedElements)
{
    return PottsElementIterator(*this, mElements.begin(), skipDeletedElements);
}


PottsMesh::PottsElementIterator PottsMesh::GetElementIteratorEnd()
{
    return PottsElementIterator(*this, mElements.end());
}

PottsElement& PottsMesh::PottsElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

PottsElement* PottsMesh::PottsElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

bool PottsMesh::PottsElementIterator::operator!=(const PottsMesh::PottsElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

PottsMesh::PottsElementIterator& PottsMesh::PottsElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    }
    while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

bool PottsMesh::PottsElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

bool PottsMesh::PottsElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}


#endif /*POTTSMESH_HPP_*/
