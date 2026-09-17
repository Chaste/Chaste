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

#ifndef SIMPLELOGARITHMICREPULSIONFORCE_HPP_
#define SIMPLELOGARITHMICREPULSIONFORCE_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"
#include "NodeBasedCellPopulation.hpp"

/**
 * A class for a simple two-body repulsion force law. Designed
 * for use in node-based simulations.
 *
 * A linear repulsive force is applied between cells that overlap
 * (i.e. whose separation is less than the sum of their radii).
 * No attractive force is applied, and no cell-age or cell-cycle
 * effects are included.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class SimpleLogarithmicRepulsionForce : public AbstractTwoBodyInteractionForce<ELEMENT_DIM, SPACE_DIM>
{
private :

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
        archive & mRepulsionParameter;
    }

protected:

    /** Repulsion parameter. Defaults to 15.0. */
    double mRepulsionParameter;

public :

    /**
     * Constructor.
     */
    SimpleLogarithmicRepulsionForce();

    /**
     * @return mRepulsionParameter
     */
    double GetRepulsionParameter();

    /**
     * Set mRepulsionParameter.
     *
     * @param repulsionParameter the new value of mRepulsionParameter
     */
    void SetRepulsionParameter(double repulsionParameter);

    /**
     * Overridden CalculateForceBetweenNodes() method.
     *
     * Returns a linear repulsive force when nodes overlap, and zero otherwise.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @return The force exerted on Node A by Node B.
     */
    c_vector<double, SPACE_DIM> CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                           unsigned nodeBGlobalIndex,
                                                           AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SimpleLogarithmicRepulsionForce)

#endif /*SIMPLELOGARITHMICREPULSIONFORCE_HPP_*/
