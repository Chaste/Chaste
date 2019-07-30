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

#ifndef ABSTRACTTETRAHEDRALMESH_HPP_
#define ABSTRACTTETRAHEDRALMESH_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "ChasteSerializationVersion.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/split_member.hpp>

#include <vector>
#include <string>
#include <cassert>
#include <boost/foreach.hpp>

#include "AbstractMesh.hpp"
#include "BoundaryElement.hpp"
#include "Element.hpp"
#include "GenericMeshReader.hpp"
#include "AbstractMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshWriter.hpp"
#include "ArchiveLocationInfo.hpp"
#include "FileFinder.hpp"


/// Forward declaration which is going to be used for friendship
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractConductivityTensors;

/**
 * Abstract base class for all tetrahedral meshes (inherits from AbstractMesh).
 */

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralMesh : public AbstractMesh<ELEMENT_DIM, SPACE_DIM>
{
    friend class AbstractConductivityTensors<ELEMENT_DIM, SPACE_DIM>; //A class which needs a global to local element mapping
    friend class CentroidWriter; //A test class which needs access to mElements in order to check that local/global indices match
protected:
    /**
     * Most tet meshes are linear (set to true).  Set to false in quadratics.
     */
    bool mMeshIsLinear;

private:
    /**
     * Pure virtual solve element mapping method. For an element with a given
     * global index, get the local index used by this process.
     * Overridden in TetrahedralMesh and DistributedTetrahedralMesh classes.
     *
     * @param index the global index of the element
     * @return local index
     *
     */
    virtual unsigned SolveElementMapping(unsigned index) const = 0;

    /**
     * Pure virtual solve boundary element mapping method. For a boundary
     * element with a given global index, get the local index used by this process.
     * Overridden in TetrahedralMesh and DistributedTetrahedralMesh classes.
     *
     * @param index the global index of the boundary element
     * @return local index
     */
    virtual unsigned SolveBoundaryElementMapping(unsigned index) const = 0;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the AbstractTetrahedralMesh. Note that this will write out a TrianglesMeshWriter file
     * to wherever ArchiveLocationInfo has specified.
     *
     * If the mesh is MutableMesh (or a subclass) the file is written by examining the current mesh.
     *
     * If the mesh is not mutable then the file is a copy of the original file the mesh was read from.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mMeshIsLinear;
        // Create a mesh writer pointing to the correct file and directory
        TrianglesMeshWriter<ELEMENT_DIM,SPACE_DIM> mesh_writer(ArchiveLocationInfo::GetArchiveRelativePath(),
                                                               ArchiveLocationInfo::GetMeshFilename(),
                                                               false);
        // Binary meshes have similar content to the original Triangle/Tetgen format, but take up less space on disk
        mesh_writer.SetWriteFilesAsBinary();

        // Archive the mesh permutation, so we can just copy the original mesh files whenever possible
        bool permutation_available = (this->rGetNodePermutation().size() != 0);
        archive & permutation_available;

        if (permutation_available)
        {
            const std::vector<unsigned>& rPermutation = this->rGetNodePermutation();
            archive & rPermutation;
        }

        if (!this->IsMeshOnDisk() || this->mMeshChangesDuringSimulation)
        {
            mesh_writer.WriteFilesUsingMesh(*(const_cast<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
        }
        else
        {
            unsigned order_of_element = (mMeshIsLinear?1:2);
            unsigned& order_of_boundary_element = order_of_element;

            // Mesh in disc, copy it to the archiving folder
            std::string original_file=this->GetMeshFileBaseName();
            std::shared_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > p_original_mesh_reader
                = GenericMeshReader<ELEMENT_DIM, SPACE_DIM>(original_file, order_of_element, order_of_boundary_element);

            if (p_original_mesh_reader->IsFileFormatBinary())
            {
                // Mesh is in binary format, we can just copy the files across ignoring the mesh reader
                if (PetscTools::AmMaster())
                {
                    FileFinder mesh_base(this->GetMeshFileBaseName());
                    FileFinder mesh_folder = mesh_base.GetParent();
                    std::string mesh_leaf_name = mesh_base.GetLeafNameNoExtension();
                    std::vector<FileFinder> mesh_files = mesh_folder.FindMatches(mesh_leaf_name + ".*");
                    FileFinder dest_dir(ArchiveLocationInfo::GetArchiveDirectory());
                    BOOST_FOREACH(const FileFinder& r_mesh_file, mesh_files)
                    {
                        FileFinder dest_file(ArchiveLocationInfo::GetMeshFilename() + r_mesh_file.GetExtension(),
                                             dest_dir);
                        ABORT_IF_THROWS(r_mesh_file.CopyTo(dest_file));
                    }
                }
            }
            else
            {
                // Mesh in text format, use the mesh writer to "binarise" it
                mesh_writer.WriteFilesUsingMeshReaderAndMesh(*p_original_mesh_reader,
                                                             *(const_cast<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>*>(this)));
            }
        }

        // Make sure that the files are written before slave processes proceed
        PetscTools::Barrier("AbstractTetrahedralMesh::save");
    }

    /**
     * Loads a mesh by using TrianglesMeshReader and the location in ArchiveLocationInfo.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractMesh<ELEMENT_DIM,SPACE_DIM> >(*this);
        archive & mMeshIsLinear;

        bool permutation_available=false;
        std::vector<unsigned> permutation;

        if (version > 0)
        {
            archive & permutation_available;

            if (permutation_available)
            {
                archive & permutation;
            }
        }

        // Store the DistributedVectorFactory loaded from the archive
        DistributedVectorFactory* p_factory = this->mpDistributedVectorFactory;
        this->mpDistributedVectorFactory = nullptr;

        // Check whether we're migrating, or if we can use the original partition for the mesh
        DistributedVectorFactory* p_our_factory = nullptr;
        if (p_factory)
        {
            p_our_factory = p_factory->GetOriginalFactory();
        }
        if (p_our_factory && p_our_factory->GetNumProcs() == p_factory->GetNumProcs())
        {
            // Specify the node distribution
            this->SetDistributedVectorFactory(p_our_factory);
        }
        else
        {
            // Migrating; let the mesh re-partition if it likes
            /// \todo #1199  make this work for everything else...
            p_our_factory = nullptr;
        }

        if (mMeshIsLinear)
        {
            // I am a linear mesh
            TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename());

            if (permutation_available)
            {
                mesh_reader.SetNodePermutation(permutation);
            }

            this->ConstructFromMeshReader(mesh_reader);
        }
        else
        {
            // I am a quadratic mesh and need quadratic information from the reader
            TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader(ArchiveLocationInfo::GetArchiveDirectory() + ArchiveLocationInfo::GetMeshFilename(), 2, 2);
            this->ConstructFromMeshReader(mesh_reader);
        }

        // Make sure we're using the correct vector factory
        if (p_factory)
        {
            if (!this->mpDistributedVectorFactory)
            {
                // If we're not using a DistributedTetrahedralMesh, ConstructFromMeshReader won't set
                // this->mpDistributedVectorFactory.
                this->mpDistributedVectorFactory = p_factory;
            }
            else
            {
                // We need to update p_factory to match this->mpDistributedVectorFactory, and then use
                // p_factory, since the rest of the code (e.g. AbstractCardiacPde) will be using p_factory.
                p_factory->SetFromFactory(this->mpDistributedVectorFactory);
                if (p_our_factory != this->mpDistributedVectorFactory)
                {
                    // Avoid memory leak
                    delete this->mpDistributedVectorFactory;
                }
                this->mpDistributedVectorFactory = p_factory;
            }
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

protected:  // Give access of these variables to subclasses

    /** Vector of pointers to elements in the mesh. */
    std::vector<Element<ELEMENT_DIM, SPACE_DIM> *> mElements;

    /** Vector of pointers to boundary elements in the mesh. */
    std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *> mBoundaryElements;

    /**
     * Sets the ownership of each element according to which nodes are owned by the
     * process.
     *
     * Information on node ownership comes from the distributed vector factory and
     * an element is "owned" if one or more of its nodes are owned
     */
    void SetElementOwnerships();

public:

    //////////////////////////////////////////////////////////////////////
    //                            Iterators                             //
    //////////////////////////////////////////////////////////////////////

    /** Definition of boundary element Iterator type. */
    typedef typename std::vector<BoundaryElement<ELEMENT_DIM-1, SPACE_DIM> *>::const_iterator BoundaryElementIterator;

    /** Forward declaration */
    class ElementIterator;

    /**
     * @return an iterator to the first element in the mesh.
     *
     * @param skipDeletedElements whether to include deleted element
     */
    inline ElementIterator GetElementIteratorBegin(bool skipDeletedElements=true);

    /**
     * @return an iterator to one past the last element in the mesh.
     */
    inline ElementIterator GetElementIteratorEnd();

    //////////////////////////////////////////////////////////////////////
    //                             Methods                              //
    //////////////////////////////////////////////////////////////////////

    /**
     * Constructor.
     */
    AbstractTetrahedralMesh();

    /**
     * Virtual destructor, since this class has virtual methods.
     */
    virtual ~AbstractTetrahedralMesh();


    /**
     * @return the number of elements that are actually in use.
     */
    virtual unsigned GetNumElements() const;

    /**
     * @return the number of local elements that are in use on this process (only over-ridden when the mesh is distributed).
     */
    virtual unsigned GetNumLocalElements() const;

    /**
     * @return the number of boundary elements that are actually in use.
     */
    virtual unsigned GetNumBoundaryElements() const;

    /**
     * @return the number of local boundary elements that are in use on this process (only over-ridden when the mesh is distributed).
     */
    virtual unsigned GetNumLocalBoundaryElements() const;

    /**
     * @return the total number of elements (including those marked as deleted).
     */
    unsigned GetNumAllElements() const;

    /**
     * @return the total number of boundary elements (including those marked as deleted).
     */
    unsigned GetNumAllBoundaryElements() const;

    /**
     * @return the number of cable elements that are actually in use.
     *
     * This will always return zero until overridden in the MixedDimensionMesh class
     */
    virtual unsigned GetNumCableElements() const;

    /**
     * @return the number of vertices (nodes which are also corners of elements).  For a linear mesh all nodes are vertices,
     * so this method is a synonym for GetNumNodes.  However, it is over-ridden in quadratic meshes.
     */
    virtual unsigned GetNumVertices() const;

    /**
     * @return the largest index of nodes on this process. Overwritten in subclasses and used for setting
     * up NodeMaps
     *
     * @return the largest node index on this process.
     */
    virtual unsigned GetMaximumNodeIndex();

    /**
     * Get the element with a given index in the mesh.
     *
     * @param index the global index of the element
     * @return a pointer to the element.
     */
    Element<ELEMENT_DIM, SPACE_DIM>* GetElement(unsigned index) const;

    /**
     * Get the boundary element with a given index in the mesh.
     *
     * @param index the global index of the boundary element
     * @return a pointer to the boundary element.
     */
    BoundaryElement<ELEMENT_DIM-1, SPACE_DIM>* GetBoundaryElement(unsigned index) const;


    /**
     * Construct the mesh using a MeshReader.
     * This method must be overridden in concrete classes.
     *
     * @param rMeshReader the mesh reader
     */
    virtual void ConstructFromMeshReader(AbstractMeshReader<ELEMENT_DIM, SPACE_DIM>& rMeshReader)=0;

    /**
     * Construct the mesh using another mesh.
     * This takes a mesh of a given concrete class and produces a deep copy.
     *
     * Use with caution when copying between subclasses.
     * @param rOtherMesh the mesh to copy
     * \todo Can we make this const?
     */
    void ConstructFromMesh(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rOtherMesh);


    /**
     * @return a pointer to the first boundary element in the mesh.
     */
    BoundaryElementIterator GetBoundaryElementIteratorBegin() const;

    /**
     * @return a pointer to *one past* the last boundary element in the mesh
     * (for consistency with STL iterators).
     */
    BoundaryElementIterator GetBoundaryElementIteratorEnd() const;

    /**
     * Compute the inverse Jacobian for a given element in the mesh.
     *
     * @param elementIndex index of an element
     * @param rJacobian  the Jacobian matrix
     * @param rJacobianDeterminant  the determinant of the Jacobian matrix
     * @param rInverseJacobian  the inverse Jacobian matrix
     */
    virtual void GetInverseJacobianForElement(unsigned elementIndex, c_matrix<double, SPACE_DIM, ELEMENT_DIM>& rJacobian,
                                              double& rJacobianDeterminant,
                                              c_matrix<double, ELEMENT_DIM, SPACE_DIM>& rInverseJacobian) const;

    /**
     * Compute the weighted direction for a given boundary element.
     *
     * @param elementIndex index of an element
     * @param rWeightedDirection the weighted direction vector
     * @param rJacobianDeterminant  the determinant of the Jacobian matrix
     *
     */
    virtual void GetWeightedDirectionForBoundaryElement(unsigned elementIndex,
                                                        c_vector<double, SPACE_DIM>& rWeightedDirection,
                                                        double& rJacobianDeterminant) const;

    /**
     * Check whether mesh has outward-facing normals.
     *
     * This throws a suitable exception if an inward facing normal is found.
     *
     * @return nothing
     */
    void CheckOutwardNormals();

    /**
     * Construct a 1D linear grid on [0,width]
     *
     * ELEMENT_DIM must be equal to 1. If SPACE_DIM > 1 then the
     * y & z default to 0.0 for every node.
     *
     * @param width  width of the mesh (in the x-direction)
     *
     * In this method the width is also THE NUMBER OF ELEMENTS IN THE x-direction.
     *
     * Overridden in DistributedTetrahedralMesh
     */
    virtual void ConstructLinearMesh(unsigned width);

    /**
     * Construct a 2D rectangular grid on [0,width]x[0,height].
     *
     * Diagonals can be staggered so that there is no preferred
     * diffusion propagation direction.
     *
     * @param width  width of the mesh (in the x-direction)
     * @param height  height of the mesh (in the y-direction)
     * @param stagger  whether the mesh should 'jumble' up the elements (defaults to true)
     *
     * In this method the width is also THE NUMBER OF ELEMENTS IN THE x-direction,
     * and similarly with the y direction.
     *
     * Overridden in DistributedTetrahedralMesh
     */
    virtual void ConstructRectangularMesh(unsigned width, unsigned height, bool stagger=true);

    /** Construct a 3D cuboid grid on [0,width]x[0,height]x[0,depth].
     *
     * @param width  width of the mesh (in the x-direction)
     * @param height  height of the mesh (in the y-direction)
     * @param depth  depth of the mesh (in the z-direction).
     *
     * In this method the width is also THE NUMBER OF ELEMENTS IN THE x-direction,
     * and similarly with the y and z directions.
     *
     * Overridden in DistributedTetrahedralMesh
     */
    virtual void ConstructCuboid(unsigned width, unsigned height, unsigned depth);

    /**
     *  Create a 1D mesh on [0, width], 2D mesh on [0, width]x[0 height] with staggering or
     *  3D mesh on [0, width]x[0 height]x[0 depth with a given axis-aligned space step.
     *  If SPACE_DIM > ELEMENT_DIM then the y & z default to 0.0 for every node.
     *
     *  @param spaceStep The axis-aligned space step
     *  @param width The width (x-dimension)
     *  @param height The height (y-dimension - ignored if ELEMENT_DIM is 1D)
     *  @param depth The depth (z-dimension -ignored in 1D and 2D)
     */
    void ConstructRegularSlabMesh(double spaceStep, double width, double height=0, double depth=0);

    /**
     *  This is a wrapper method to ConstructRegularSlabMesh() which is useful for parallel distributed meshes.
     *  By default slabs are split across processes in the top dimension (y in 2d, z in 3d)
     *  but it may be more useful to split in x or y.
     *
     *  Create a 1D mesh on [0, width], 2D mesh on [0, width]x[0 height] with staggering or
     *  3D mesh on [0, width]x[0 height]x[0 depth with a given axis-aligned space step.
     *  If SPACE_DIM > ELEMENT_DIM then the y & z default to 0.0 for every node.
     *
     *  @param dimension The dimension/axis to be split over the processes.  When dimension=0 the split is on x, dimension=1
     *  indicates split on y etc.  If dimension >= SPACE_DIM then an exception is thrown.
     *  @param spaceStep The axis-aligned space step
     *  @param width The width (x-dimension)
     *  @param height The height (y-dimension - ignored if ELEMENT_DIM is 1D)
     *  @param depth The depth (z-dimension -ignored in 1D and 2D)
     */
    void ConstructRegularSlabMeshWithDimensionSplit(unsigned dimension, double spaceStep, double width, double height=0, double depth=0);


    /**
     * Determine whether or not the current process owns node 0 of this boundary element (tie breaker to determine which process writes
     * to file for when two or more share ownership of a face).
     *
     * @param faceIndex is the global index of the face
     * @return true if this process is designated owner
     */
    virtual bool CalculateDesignatedOwnershipOfBoundaryElement( unsigned faceIndex );

    /**
     * Determine whether or not the current process owns node 0 of this element (tie breaker to determine which process writes
     * to file for when two or more share ownership of an element).
     *
     * @param elementIndex is the global index of the element
     * @return true if this process is designated owner
     */
    virtual bool CalculateDesignatedOwnershipOfElement( unsigned elementIndex );

    /**
     * @return Iterates through local nodes and finds two representative nodes with a maximum number of
     * containing elements for all locally owned nodes.  The two representative nodes are the most connected
     * nodes at the boundary and the most connected node in the interior.  This is because a well-element-connected
     * boundary node is likely to be connected to more nodes since fewer of them will overlap.)  At these
     * representative node the node connectivity (number of nodes in forward star) is determined.
     *
     * Useful for determining FEM matrix fill.
     */
    unsigned CalculateMaximumNodeConnectivityPerProcess() const;

    /**
     * Utility method to give the functionality of iterating through the halo nodes of a process. Will return an empty
     * std::vector (i.e. no halo nodes) unless overridden by distributed derived classes.
     *
     * @param rHaloIndices  A vector to fill with the global indices of the nodes which are locally halos
     */
    virtual void GetHaloNodeIndices(std::vector<unsigned>& rHaloIndices) const;

    /**
     * Get the nodes which will need to be exchanged between remote processes.
     * If we have an element which node indices outside the local [mLo, mHi) region
     * then we know that those nodes will need to be recieved from a remote process, while
     * those inside the range [mLo, mHi) will need to be sent
     *
     * @param rNodesToSendPerProcess (output) a vector which will be of size GetNumProcs()
     * where each internal vector except i=GetMyRank() contains an ordered list of indices of
     * nodes to send to process i
     *
     * @param rNodesToReceivePerProcess (output) a vector which will be of size GetNumProcs()
     * for information to receive for process i
     */
     void CalculateNodeExchange( std::vector<std::vector<unsigned> >& rNodesToSendPerProcess,
                                 std::vector<std::vector<unsigned> >& rNodesToReceivePerProcess);


     /**
      * Computes the minimum and maximum lengths of the edges in the mesh.
      * Overridden in Distributed case
      * \todo Should be const
      *
      * @return The minimum and maximum edge lengths in the mesh
      *
      */
     virtual c_vector<double, 2> CalculateMinMaxEdgeLengths();

     /**
      * Return the element index for the first element that contains a test point
      *
      * @param rTestPoint reference to the point
      * @param strict  Should the element returned contain the point in the interior and
      *      not on an edge/face/vertex (default = not strict)
      * @param testElements  a set of guesses for the element (a set of element indices), to be checked
      *      first for potential efficiency improvements. (default = empty set)
      * @param onlyTryWithTestElements Do not continue with other elements after trying the with testElements
      *      (for cases where you know the testPoint must be in the set of test elements or maybe outside
      *      the mesh).
      * @return element index
      */
     unsigned GetContainingElementIndex(const ChastePoint<SPACE_DIM>& rTestPoint,
                                        bool strict=false,
                                        std::set<unsigned> testElements=std::set<unsigned>(),
                                        bool onlyTryWithTestElements = false);

     /** As with GetNearestElementIndex() except only searches in the given set of elements.
      * @param rTestPoint reference to the point
      * @param testElements a set of elements (element indices) to look in
      * @return element index
      */
     unsigned GetNearestElementIndexFromTestElements(const ChastePoint<SPACE_DIM>& rTestPoint,
                                                     std::set<unsigned> testElements);

    //////////////////////////////////////////////////////////////////////
    //                         Nested classes                           //
    //////////////////////////////////////////////////////////////////////

    /**
     * A smart iterator over the elements in the mesh.
     */
    class ElementIterator
    {
    public:
        /**
         * Dereference the iterator giving you a *reference* to the current element.
         * @return reference
         * Make sure to use a reference for the result to avoid copying elements unnecessarily.
         */
        inline Element<ELEMENT_DIM, SPACE_DIM>& operator*();

        /**
         * Member access from a pointer.
         * @return pointer
         */
        inline Element<ELEMENT_DIM, SPACE_DIM>* operator->();

        /**
         * @return Comparison not-equal-to.
         *
         * @param rOther iterator with which comparison is made
         */
        inline bool operator!=(const typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator& rOther);

        /**
         * Prefix increment operator.
         * @return reference to incremented object
         */
        inline ElementIterator& operator++();

        /**
         * Constructor for a new iterator.
         *
         * This should not be called directly by user code; use the mesh methods
         * AbstractTetrahedralMesh::GetElementIteratorBegin and AbstractTetrahedralMesh::GetElementIteratorEnd instead.
         *
         * @param rMesh the mesh to iterator over
         * @param elementIter where to start iterating
         * @param skipDeletedElements whether to include deleted elements
         */
        ElementIterator(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                        typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
                        bool skipDeletedElements=true);

    private:
        /** The mesh we're iterating over. */
        AbstractTetrahedralMesh& mrMesh;

        /** The actual element iterator. */
        typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::iterator mElementIter;

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

TEMPLATED_CLASS_IS_ABSTRACT_2_UNSIGNED(AbstractTetrahedralMesh)

namespace boost {
namespace serialization {
/**
 * Specify a version number for archive backwards compatibility.
 *
 * This is how to do BOOST_CLASS_VERSION(AbstractCardiacTissue, 1)
 * with a templated class.
 */
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
struct version<AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> >
{
    ///Macro to set the version number of templated archive in known versions of Boost
    CHASTE_VERSION_CONTENT(1);
};
} // namespace serialization
} // namespace boost


//////////////////////////////////////////////////////////////////////////////
//      ElementIterator class implementation - most methods are inlined     //
//////////////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorBegin(
        bool skipDeletedElements)
{
    return ElementIterator(*this, mElements.begin(), skipDeletedElements);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::GetElementIteratorEnd()
{
    return ElementIterator(*this, mElements.end());
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>& AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::operator*()
{
    assert(!IsAtEnd());
    return **mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Element<ELEMENT_DIM, SPACE_DIM>* AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::operator->()
{
    assert(!IsAtEnd());
    return *mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::operator!=(const typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator& rOther)
{
    return mElementIter != rOther.mElementIter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
typename AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator& AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::operator++()
{
    do
    {
        ++mElementIter;
    }
    while (!IsAtEnd() && !IsAllowedElement());

    return (*this);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::ElementIterator(
        AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
        typename std::vector<Element<ELEMENT_DIM, SPACE_DIM> *>::iterator elementIter,
        bool skipDeletedElements)
    : mrMesh(rMesh),
      mElementIter(elementIter),
      mSkipDeletedElements(skipDeletedElements)
{
    if (mrMesh.mElements.size() == 0)
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
bool AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::IsAtEnd()
{
    return mElementIter == mrMesh.mElements.end();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>::ElementIterator::IsAllowedElement()
{
    return !(mSkipDeletedElements && (*this)->IsDeleted());
}

#endif /*ABSTRACTTETRAHEDRALMESH_HPP_*/
