/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef DIFFUSIONCAUPDATERULE_HPP_
#define DIFFUSIONCAUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCaUpdateRule.hpp"
#include "CaBasedCellPopulation.hpp"

/**
 * An update rule class to model diffusion on the lattice.
 *
 * A diffusion coefficient (defaulting to 1.0) is passed into the constructor.
 * At each time step, a uniform random number r_i is generated for each cell i;
 * if r_i > D*dt, where D denotes the diffusion coefficient and dt denotes the
 * time step, then the cell is moved to a random free neighbouring site (if there
 * are no free neighbouring sites, then the cell is not moved).
 */
template<unsigned DIM>
class DiffusionCaUpdateRule : public AbstractCaUpdateRule<DIM>
{
private:

    /** Diffusion constant. */
    double mDiffusionConstant;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCaUpdateRule<DIM> >(*this);
        archive & mDiffusionConstant;
    }

public:

    /**
     * Constructor.
     *
     * @param diffusionConstant the diffusion constant (defaults to 1.0)
     */
    DiffusionCaUpdateRule(double diffusionConstant=1.0);

    /**
     * Destructor.
     */
    ~DiffusionCaUpdateRule();

    /**
     * Overridden GetNewLocationOfCell() method.
     *
     * This randomly moves the cell in one of the possible directions defined by the cell population method
     * GetNeighbouringNodeIndices().
     *
     * @param currentLocationIndex the current location index corresponding to a cell
     * @param rCellPopulation reference to the cell population
     * @param dt simulation time step, used to calculate the probability of movement
     *
     * @return the new location index of the cell
     */
    unsigned GetNewLocationOfCell(unsigned currentLocationIndex,
                                  CaBasedCellPopulation<DIM>& rCellPopulation,
                                  double dt);

    /** @return mDiffusionConstant */
    double GetDiffusionConstant();

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DiffusionCaUpdateRule)

#endif /*DIFFUSIONCAUPDATERULE_HPP_*/
