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

#ifndef CYLINDRICAL2DNODESONLYMESH_HPP_
#define CYLINDRICAL2DNODESONLYMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "NodesOnlyMesh.hpp"

/**
 * A subclass of NodesOnlyMesh<2> for a rectangular mesh with
 * periodic left and right boundaries, representing a cylindrical geometry.
 *
 * The class works by overriding calls such as ReMesh() and
 * GetVectorFromAtoB() so that simulation classes can treat this
 * class in exactly the same way as a NodesOnlyMesh<2>.
 */
class Cylindrical2dNodesOnlyMesh: public NodesOnlyMesh<2>
{
private:
    /**
     * The periodic width of the domain
     */
    double mWidth;

    friend class TestCylindrical2dNodesOnlyMesh;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archives the member variables of the object which have to be preserved
     * during its lifetime.
     *
     * Note that we must archive any member variables FIRST so that this
     * method can call a ReMesh (to convert from TrianglesMeshReader input
     * format into our native format).
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<NodesOnlyMesh<2> >(*this);
        archive & mWidth;
    }

public:

    /**
     * Constructor.
     *
     * @param width the width of the mesh (circumference)
     */
    Cylindrical2dNodesOnlyMesh(double width);

    /**
     * Set up the box collection
     *
     * @param cutOffLength the cut off length for node neighbours
     * @param domainSize the size of the domain containing the nodes.
     * @param numLocalRows the number of rows of the collection that this process should own.
     * @param isPeriodic whether the box collection should be periodic. Defaults to true.
     */
    virtual void SetUpBoxCollection(double cutOffLength, c_vector<double, 2*2> domainSize, int numLocalRows = PETSC_DECIDE, bool isPeriodic = true);

    /**
     * Overridden GetVectorFromAtoB() method.
     *
     * Evaluates the (surface) distance between two points in a 2D cylindrical
     * geometry.
     *
     * @param rLocation1 the x and y co-ordinates of point 1
     * @param rLocation2 the x and y co-ordinates of point 2
     * @return the vector from location1 to location2
     */
    c_vector<double, 2> GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2);

    /**
     * Overridden GetWidth() method.
     *
     * Calculate the 'width' of any dimension of the mesh, taking periodicity
     * into account.
     *
     * @param rDimension a dimension (0 or 1)
     * @return The maximum distance between any nodes in this dimension.
     */
    double GetWidth(const unsigned& rDimension) const;

    /**
     * Overridden SetNode() method.
     *
     * If the location should be set outside a cylindrical boundary
     * move it back onto the cylinder.
     *
     * @param nodeIndex is the index of the node to be moved
     * @param point is the new target location of the node
     * @param concreteMove is set to false if we want to skip the signed area tests in the parent Class Note this should always be false here
     */
    void SetNode(unsigned nodeIndex, ChastePoint<2> point, bool concreteMove = false);

    /**
     * Overridden AddNode() method.
     *
     * @param pNewNode  pointer to the new node
     * @return index of new node
     */
    unsigned AddNode(Node<2>* pNewNode);

    /**
     * Overridden RefreshMesh() method.
     *
     * If the location is outside the domain width, move it to
     * within the boundary
     */
    void RefreshMesh();
};

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Cylindrical2dNodesOnlyMesh.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Cylindrical2dNodesOnlyMesh * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const double width = t->GetWidth(0);
    ar & width;
}

/**
 * De-serialize constructor parameters and initialise a Cylindrical2dNodesOnlyMesh.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Cylindrical2dNodesOnlyMesh * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    double width;
    ar & width;

    // Invoke inplace constructor to initialise instance
    ::new(t)Cylindrical2dNodesOnlyMesh(width);
}
}
} // namespace ...

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Cylindrical2dNodesOnlyMesh)

#endif /*CYLINDRICAL2DNODESONLYMESH_HPP_*/
