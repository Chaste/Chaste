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

#ifndef ABSTRACTCASWITCHINGUPDATERULE_HPP_
#define ABSTRACTCASWITCHINGUPDATERULE_HPP_

#include "AbstractUpdateRule.hpp"
#include "CaBasedCellPopulation.hpp"

template<unsigned DIM>
class PottsBasedCellPopulation; // Circular definition

/**
 * An abstract CA 'switching' update rule class, for use in CA simulations.
 */
template<unsigned DIM>
class AbstractCaSwitchingUpdateRule : public AbstractUpdateRule<DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractUpdateRule<DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    AbstractCaSwitchingUpdateRule();

    /**
     * Destructor.
     */
    virtual ~AbstractCaSwitchingUpdateRule();

   /**
     * Calculate the probability of a given switch.
     *
     * @param currentNodeIndex The index of the current node/lattice site
     * @param neighbourNodeIndex The index of the neighbour node/lattice site
     * @param rCellPopulation The cell population
     * @param dt is the time interval
     * @param deltaX defines the size of the lattice site
     * @return The probability of the cells associated to the current node and the target node switching.
     */
    double virtual EvaluateSwitchingProbability(unsigned currentNodeIndex,
                                               unsigned neighbourNodeIndex,
                                               CaBasedCellPopulation<DIM>& rCellPopulation,
                                               double dt,
                                               double deltaX)=0;

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile a file stream
     */
    virtual void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#endif /*ABSTRACTCASWITCHINGUPDATERULE_HPP_*/
