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

#ifndef PERIODICNODESONLYMESH_HPP_
#define PERIODICNODESONLYMESH_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "NodesOnlyMesh.hpp"

/**
 * A subclass of NodesOnlyMesh<SPACE_DIM> for a rectangular mesh with
 * periodic boundaries.
 *
 * The class works by overriding calls such as ReMesh() and
 * GetVectorFromAtoB() so that simulation classes can treat this
 * class in exactly the same way as a NodesOnlyMesh<SPACE_DIM>.
 */
template<unsigned SPACE_DIM>
class PeriodicNodesOnlyMesh: public NodesOnlyMesh<SPACE_DIM>
{
private:
    /**
     * The periodic widths of the domain
     */
    c_vector<double,SPACE_DIM> mWidth;

    /*
     * Whether a dimension is periodic
     */
    c_vector<bool,SPACE_DIM> mIsDimPeriodic;

    friend class TestPeriodicNodesOnlyMesh;

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
        archive & boost::serialization::base_object<NodesOnlyMesh<SPACE_DIM> >(*this);
        archive & mWidth;
        archive & mIsDimPeriodic;
    }

public:

    /**
     * Constructor.
     *
     * @param width the periodic widths of the mesh. (Note if width=0.0 then not periodic)
     *
     */
    PeriodicNodesOnlyMesh(c_vector<double,SPACE_DIM> width);

    /**
     * Set up the box collection
     *
     * @param cutOffLength the cut off length for node neighbours
     * @param domainSize the size of the domain containing the nodes.
     * @param numLocalRows the number of rows of the collection that this process should own.
     * @param isDimPeriodic whether the box collection should be periodic in the x y and/or z direction.
     */
    virtual void SetUpBoxCollection(double cutOffLength, c_vector<double, 2*SPACE_DIM> domainSize, int numLocalRows = PETSC_DECIDE,  c_vector<bool,SPACE_DIM> isDimPeriodic = zero_vector<bool>(SPACE_DIM));

    /**
     * Overridden GetVectorFromAtoB() method.
     *
     * Evaluates the (surface) distance between two points in ND periodic
     * geometry.
     *
     * @param rLocation1 the x y and z co-ordinates of point 1
     * @param rLocation2 the x y and z co-ordinates of point 2
     * @return the vector from location1 to location2
     */
    c_vector<double, SPACE_DIM> GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocation1, const c_vector<double, SPACE_DIM>& rLocation2);

    /**
     * Overridden GetWidth() method.
     *
     * Calculate the 'width' of any dimension of the mesh, taking periodicity
     * into account.
     *
     * @param rDimension a dimension (0,1 or 2)
     * @return The maximum distance between any nodes in this dimension.
     */
    double GetWidth(const unsigned& rDimension) const;

    /*
     * Method to get the width of the periodic dimensions
     * @return The value of the mWidth vector
    */
    c_vector<double,SPACE_DIM> GetPeriodicWidths() const;

    /*
     * Overridden SetNode() method.
     *
     * If the location should be set outside a cylindrical boundary
     * move it back onto the cylinder.
     *
     * @param nodeIndex is the index of the node to be moved
     * @param point is the new target location of the node
     * @param concreteMove is set to false if we want to skip the signed area tests in the parent Class Note this should always be false here
     */

    void SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point, bool concreteMove = false);

    /**
     * Overridden AddNode() method.
     *
     * @param pNewNode  pointer to the new node
     * @return index of new node
     */
    unsigned AddNode(Node<SPACE_DIM>* pNewNode);

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
 * Serialize information required to construct a PeriodicNodesOnlyMesh.
 */
template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const PeriodicNodesOnlyMesh<SPACE_DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const c_vector<double,SPACE_DIM> width = t->GetPeriodicWidths();
    for ( unsigned i=0; i < SPACE_DIM; i++ )
    {
        ar << width(i);
    }
}

/**
 * De-serialize constructor parameters and initialise a PeriodicNodesOnlyMesh.
 */
template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, PeriodicNodesOnlyMesh<SPACE_DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance of the mesh
    c_vector<double,SPACE_DIM> width;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        double current_width;
        ar & current_width;
        width(i) = current_width;
    }

    // Invoke inplace constructor to initialise instance
    ::new(t)PeriodicNodesOnlyMesh<SPACE_DIM>(width);
}
}
} // namespace ...

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PeriodicNodesOnlyMesh)

#endif /*PERIODICNODESONLYMESH_HPP_*/
