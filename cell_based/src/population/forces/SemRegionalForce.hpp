/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef SEMREGIONALFORCE_HPP_
#define SEMREGIONALFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A linear spring force for Subcellular Element (SEM) simulations in which the spring constant
 * and rest length depend on the region label of the interacting nodes.
 *
 * Node regions are set by the mesh generator using the SemNodeRegion enum:
 *   SEM_INTERIOR_REGION (= 0): subcellular nodes in the cell interior
 *   SEM_BOUNDARY_REGION (= 1): subcellular nodes on the cell surface / cortex
 *
 * The effective spring constant and rest length for a pair of interacting nodes is the arithmetic
 * mean of each node's per-region value, allowing smooth transitions at interior-boundary contacts.
 *
 * Forces are applied to node pairs whose separation is less than 0.5 (in mesh units). The force
 * is zero beyond that distance.
 */
template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class SemRegionalForce : public AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM>
{
    friend class TestForces;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mSpringConstants;
        archive & mRestLengths;
    }

protected:

    /** Spring constants indexed by SemNodeRegion: [SEM_INTERIOR_REGION, SEM_BOUNDARY_REGION, ...]. */
    std::vector<double> mSpringConstants = {1.0, 2.0, 3.0};

    /** Rest lengths indexed by SemNodeRegion: [SEM_INTERIOR_REGION, SEM_BOUNDARY_REGION, ...]. */
    std::vector<double> mRestLengths = {0.2, 0.15, 0.1};


public:

    /**
     * Constructor.
     */
    SemRegionalForce();

    /**
     * Destructor.
     */
    virtual ~SemRegionalForce() = default;

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Calculates the force between two nodes using region-dependent spring
     * constants and rest lengths.
     *
     * Note that this assumes the nodes are neighbours and is called by
     * AddForceContribution().
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     *
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                     unsigned nodeBGlobalIndex,
                                                     AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile) override;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SemRegionalForce)

#endif /*SEMREGIONALFORCE_HPP_*/
