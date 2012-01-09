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

#ifndef _NODE_HPP_
#define _NODE_HPP_

#include "UblasVectorInclude.hpp"

#include <set>
#include <vector>

#include "ChastePoint.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AbstractTetrahedralMesh;

/**
 * A node in a finite element mesh.
 */
template<unsigned SPACE_DIM>
class Node
{
private:

    /** The index of this node within the mesh. */
    unsigned mIndex;

    /** A region ID. */
    unsigned mRegion;

    /** The location of this node within the mesh. */
    c_vector<double, SPACE_DIM> mLocation;

    /** Whether this node is a boundary node. */
    bool mIsBoundaryNode;

    /** Whether this node is an internal node. (For use with QuadraticMesh)*/
    bool mIsInternal;

    /**
     * Whether this node has been deleted, and hence
     * whether its location in the mesh can be re-used. (For use in MutableMesh)
     */
    bool mIsDeleted;

    /** Set of indices of elements containing this node as a vertex. */
    std::set<unsigned> mElementIndices;

    /** Set of indices of boundary elements containing this node as a vertex. */
    std::set<unsigned> mBoundaryElementIndices;

    /** Set of node attributes*/
    std::vector<double> mNodeAttributes;

    /**
     * Extraction of commonality between the constructors.
     *
     * @param index  the index of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node
     */
    void CommonConstructor(unsigned index, bool isBoundaryNode);

public:

    /**
     * There are many ways of creating a node, depending on how you wish to specify its
     * spatial location.
     */

    /**
     * Constructor which takes the node's location as a ChastePoint.
     *
     * @param index  the index of the node in the mesh
     * @param point  the location of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     */
    Node(unsigned index, ChastePoint<SPACE_DIM> point, bool isBoundaryNode=false);

    /**
     * Constructor which takes the node's location as a std::vector.
     *
     * @param index  the index of the node in the mesh
     * @param coords  the location of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     */
    Node(unsigned index, std::vector<double> coords, bool isBoundaryNode=false);

    /**
     * Constructor which takes the node's location as a c_vector.
     *
     * @param index  the index of the node in the mesh
     * @param location  the location of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     */
    Node(unsigned index, c_vector<double, SPACE_DIM> location, bool isBoundaryNode=false);

    /**
     * Constructor which takes the coordinates of the node's location as separate input arguments.
     *
     * @param index  the index of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     * @param v1 the x-coordinate of the node in the mesh (defaults to 0)
     * @param v2 the y-coordinate of the node in the mesh (defaults to 0)
     * @param v3 the z-coordinate of the node in the mesh (defaults to 0)
     */
    Node(unsigned index, bool isBoundaryNode=false, double v1=0, double v2=0, double v3=0);

    /**
     * Constructor which takes the coordinates of the node's location as array pointer.
     *
     * @param index  the index of the node in the mesh
     * @param isBoundaryNode  whether the node is a boundary node (defaults to false)
     * @param location address of the x-coordinate of the node in the mesh
     * (other coordinates are assumed to be in contiguous memory)
     */
    Node(unsigned index,  double *location, bool isBoundaryNode=false);

    /**
     * Set the node's location.
     *
     * Note: setting the point in space is dangerous.
     * Jacobian and JacobianDeterminant of element need to be updated.
     *
     * @param point the new location of the node
     */
    void SetPoint(ChastePoint<SPACE_DIM> point);

    /**
     * Set the index of this node in the mesh.
     *
     * @param index the index of the node
     */
    void SetIndex(unsigned index);

    /**
     * Add an attribute to the list of node attributes.
     *
     * @param attribute: the value of the attribute to be added
     */
    void AddNodeAttribute(double attribute);

    /**
     * Set whether this node is a boundary node.
     *
     * @param value whether the node is a boundary node
     */
    void SetAsBoundaryNode(bool value=true);

    /**
     * Get the node's location as a ChastePoint.
     */
    ChastePoint<SPACE_DIM> GetPoint() const;

    /**
     * Get the node's location as a c_vector.
     *
     * The returned location may not be modified; if you want that functionality use
     * rGetModifiableLocation instead.
     */
    const c_vector<double, SPACE_DIM>& rGetLocation() const;

    /**
     * Get the node's location as a c_vector.
     *
     * If you modify the returned location,
     * Jacobian and JacobianDeterminant of elements need to be updated.
     *
     * Don't forget to assign the result of this call to a reference!
     */
    c_vector<double, SPACE_DIM>& rGetModifiableLocation();

    /**
     * Get the index of this node in the mesh.
     */
    unsigned GetIndex() const;

    /**
     * Get whether this node is a boundary node.
     */
    bool IsBoundaryNode() const;

    /**
     * Add an element that contains this node.
     *
     * @param index of the element to add.
     */
    void AddElement(unsigned index);

    /**
     * Remove an element that contains this node.
     *
     * @param index of the element to be removed.
     */
    void RemoveElement(unsigned index);

    /**
     * Remove an boundary element that contains this node.
     *
     * @param index of the boundary element to be removed.
     */
    void RemoveBoundaryElement(unsigned index);

    /**
     * Add an boundary element that contains this node.
     *
     * @param index of the element to add.
     */
    void AddBoundaryElement(unsigned index);

    /**
     * Return a set of indices of elements containing this node as a vertex.
     */
    std::set<unsigned>& rGetContainingElementIndices();

    /**
     * Returs a vector containing the ndoe attributes.
     */
    std::vector<double>& rGetNodeAttributes();

    /**
     * Return a set of indices of boundary elements containing this node as a vertex.
     */
    std::set<unsigned>& rGetContainingBoundaryElementIndices();

    /**
     * Get the number of elements in the mesh that contain this node.
     */
    unsigned GetNumContainingElements() const;

    /**
     * Get the number of boundary elements in the mesh that contain this node.
     */
    unsigned GetNumBoundaryElements() const;

    /**
     * Mark the node as having been removed from the mesh.
     */
    void MarkAsDeleted();

    /**
     * Get whether the node is marked as deleted.
     */
    bool IsDeleted() const;

    /**
     * Mark the node as being internal (not vertex) in a quadratic element.
     */
    void MarkAsInternal();

    /**
     * Get whether the node is internal (not vertex) in a quadratic element.
     */
    bool IsInternal() const;

    /**
     * Determine if a node lives within a flagged element.
     *
     * @param rMesh the mesh containing the nodes and elements
     */
    template <unsigned ELEMENT_DIM>
    bool IsFlagged(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
    {
        bool in_flagged_element = false;
        for (ContainingElementIterator it = ContainingElementsBegin();
             it != ContainingElementsEnd();
             ++it)
        {
            if (rMesh.GetElement(*it)->IsFlagged())
            {
                in_flagged_element = true;
                break;
            }
        }
        return in_flagged_element;
    }

    /**
     * Set the node's region ID.
     *
     * @param region the new region ID
     */
    void SetRegion(unsigned region);

    /**
     * Get the node's region ID.
     */
    unsigned GetRegion() const;

    /**
     * An iterator over the indices of elements which contain this node.
     */
    class ContainingElementIterator
    {
    public:
        /**
         * Constructor for a new ContainingElementIterator.
         *
         * @param indexIterator  an index iterator
         */
        ContainingElementIterator(std::set<unsigned>::const_iterator indexIterator)
            : mIndexIterator(indexIterator)
        {}
        /**
         * Prefix dereference operator.
         */
        const unsigned& operator*() const
        {
            return *mIndexIterator;
        }
        /**
         * Comparison not-equal-to.
         *
         * @param rOther ContainingElementIterator with which comparison is made
         */
        bool operator!=(const ContainingElementIterator& rOther) const
        {
            return mIndexIterator != rOther.mIndexIterator;
        }
        /**
         * Comparison equal-to.
         *
         * @param rOther ContainingElementIterator with which comparison is made
         */
        bool operator==(const ContainingElementIterator& rOther) const
        {
            return !operator!=(rOther);
        }
        /**
         * Prefix increment operator.
         */
        ContainingElementIterator& operator++()
        {
            ++mIndexIterator;
            return *this;
        }
    private:
        std::set<unsigned>::const_iterator mIndexIterator;  /**< Element index iterator. */
    };

    /**
     * Get a ContainingElementIterator pointing to the first containing element
     */
    ContainingElementIterator ContainingElementsBegin() const
    {
        return ContainingElementIterator(mElementIndices.begin());
    }

    /**
     * Get a ContainingElementIterator pointing to one past the last containing element
     */
    ContainingElementIterator ContainingElementsEnd() const
    {
        return ContainingElementIterator(mElementIndices.end());
    }

    /**
     * An iterator over the indices of boundary elements which contain this node.
     */
    class ContainingBoundaryElementIterator
    {
    public:
        /**
         * Constructor for a new ContainingBoundaryElementIterator.
         *
         * @param indexIterator  an index iterator
         */
        ContainingBoundaryElementIterator(std::set<unsigned>::const_iterator indexIterator)
            : mIndexIterator(indexIterator)
        {}
        /**
         * Prefix dereference operator.
         */
        const unsigned& operator*() const
        {
            return *mIndexIterator;
        }
        /**
         * Comparison not-equal-to.
         *
         * @param rOther ContainingBoundaryElementIterator with which comparison is made
         */
        bool operator!=(const ContainingBoundaryElementIterator& rOther) const
        {
            return mIndexIterator != rOther.mIndexIterator;
        }
        /**
         * Comparison equal-to.
         *
         * @param rOther ContainingBoundaryElementIterator with which comparison is made
         */
        bool operator==(const ContainingBoundaryElementIterator& rOther) const
        {
            return !operator!=(rOther);
        }
        /**
         * Prefix increment operator.
         */
        ContainingBoundaryElementIterator& operator++()
        {
            ++mIndexIterator;
            return *this;
        }
    private:
        std::set<unsigned>::const_iterator mIndexIterator;  /**< Boundary element index iterator. */
    };

    /**
     * Get a ContainingBoundaryElementIterator pointing to the first containing boundary element
     */
    ContainingBoundaryElementIterator ContainingBoundaryElementsBegin() const
    {
        return ContainingBoundaryElementIterator(mBoundaryElementIndices.begin());
    }

    /**
     * Get a ContainingBoundaryElementIterator pointing to one past the last containing boundary element
     */
    ContainingBoundaryElementIterator ContainingBoundaryElementsEnd() const
    {
        return ContainingBoundaryElementIterator(mBoundaryElementIndices.end());
    }
};


#endif //_NODE_HPP_
