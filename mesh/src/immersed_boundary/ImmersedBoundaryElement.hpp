/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef IMMERSEDBOUNDARYELEMENT_HPP_
#define IMMERSEDBOUNDARYELEMENT_HPP_

#include "MutableElement.hpp"
#include "FluidSource.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * An element class for use in the ImmersedBoundaryMesh class. The main
 * difference between this and the Element class is that a
 * ImmersedBoundaryElement can have a variable number of nodes associated
 * with it.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class ImmersedBoundaryElement : public MutableElement<ELEMENT_DIM, SPACE_DIM>
{
private:

    /** Fluid source associated with this element. */
    std::shared_ptr<FluidSource<SPACE_DIM>> mpFluidSource;

    /** Corner nodes associated with this element. */
    std::vector<Node<SPACE_DIM>*> mCornerNodes;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /** The average node spacing. */
    double mAverageNodeSpacing;

    /** Whether this element is on the boundary. This is a non-local property so must be calculated by the mesh. */
    bool mIsBoundaryElement;

    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This needs to be first so that MeshBasedCellPopulation::Validate() doesn't go mental.
        archive & boost::serialization::base_object<MutableElement<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param index global index of the element
     * @param rNodes vector of Nodes associated with the element
     */
    ImmersedBoundaryElement(unsigned index,
                            const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Destructor.
     */
    virtual ~ImmersedBoundaryElement();

    /**
     * Set mpFluidSource.
     *
     * @param fluidSource  the fluid source associated with this element
     */
    void SetFluidSource(std::shared_ptr<FluidSource<SPACE_DIM>> fluidSource);

    /**
     * Get the fluid source associated with this element
     *
     * @return pointer to the fluid source
     */
    std::shared_ptr<FluidSource<SPACE_DIM>> GetFluidSource();

    /**
     * @return the vector of corner nodes.
     *
     * Don't forget to assign the result of this call to a reference!
     */
    std::vector<Node<SPACE_DIM>*>& rGetCornerNodes();

    /**
     * @return the average node spacing.
     */
    double GetAverageNodeSpacing();

    /**
     * @param averageNodeSpacing the new average node spacing.
     */
    void SetAverageNodeSpacing(double averageNodeSpacing);

    /**
     * Overridden method to get whether or not the element is on the boundary.
     *
     * @return whether or not the element is on the boundary.
     */
    virtual bool IsElementOnBoundary() const;

    /** @param isBoundaryElement whether the element is on the boundary */
    void SetIsBoundaryElement(bool isBoundaryElement);

    /**
     * @param node the node to add as a corner node
     */
    void AddCornerNode(Node<SPACE_DIM>* node);
};

//////////////////////////////////////////////////////////////////////
//                  Specialization for 1d elements                  //
//                                                                  //
//                 1d elements are just edges (lines)               //
//////////////////////////////////////////////////////////////////////

/**
 * Specialization for 1d elements so we don't get errors from Boost on some
 * compilers.
 */
template<unsigned SPACE_DIM>
class ImmersedBoundaryElement<1, SPACE_DIM> : public MutableElement<1,SPACE_DIM>
{
private:

    /** Fluid source associated with this element. */
    std::shared_ptr<FluidSource<SPACE_DIM>> mpFluidSource;

    /** Corner nodes associated with this element. */
    std::vector<Node<SPACE_DIM>*> mCornerNodes;

    /** The average node spacing. */
    double mAverageNodeSpacing;

    /** Whether this element is on the boundary. This is a non-local property so must be calculated by the mesh. */
    bool mIsBoundaryElement;

public:

    /**
     * Constructor which takes in a vector of nodes.
     *
     * @param index  the index of the element in the mesh
     * @param rNodes the nodes owned by the element
     */
    ImmersedBoundaryElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Set mpFluidSource
     *
     * @param fluidSource a shared pointer to the fluid source associated with this element
     */
    void SetFluidSource(std::shared_ptr<FluidSource<SPACE_DIM>> fluidSource);

    /**
     * Get the fluid source associated with this element
     *
     * @return pointer to the fluiod source
     */
    std::shared_ptr<FluidSource<SPACE_DIM>> GetFluidSource(void);

    /**
     * @return the vector of corner nodes.
     *
     * Don't forget to assign the result of this call to a reference!
     */
    std::vector<Node<SPACE_DIM>*>& rGetCornerNodes();

    /**
     * @return the average node spacing.
     */
    double GetAverageNodeSpacing();

    /**
     * @param averageNodeSpacing the new average node spacing.
     */
    void SetAverageNodeSpacing(double averageNodeSpacing);

    /**
     * Overridden method to get whether or not the element is on the boundary.
     *
     * @return false.
     */
    virtual bool IsElementOnBoundary() const;

    /** @param isBoundaryElement whether the element is on the boundary */
    void SetIsBoundaryElement(bool isBoundaryElement);

};

#endif /*IMMERSEDBOUNDARYELEMENT_HPP_*/
