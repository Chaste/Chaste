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

#ifndef ADVECTIONCAUPDATERULE_HPP_
#define ADVECTIONCAUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCaUpdateRule.hpp"
#include "CaBasedCellPopulation.hpp"

/**
 * An update rule class to model uniform steady advection. This class is
 * currently implemented in 2D only.
 *
 * The constructor must be passed an unsigned that corresponds to one of
 * the eight possible flow directions in 2D (N, NW, W,SW, S, SE, E, NE)
 * and a double that gives the speed of the flow. At each time step, a
 * uniform random number r_i is generated for each cell i that has a free
 * neighbour in the direction of the imposed flow; if r_i < s*dt, where
 * s denotes the flow speed and dt denotes the time step, then the cell is
 * moved. This ensures that the mean speed of an advected cell is equal to
 * the flow speed, independent of the time step used.
 */
template<unsigned DIM>
class AdvectionCaUpdateRule : public AbstractCaUpdateRule<DIM>
{
private:

    /**
     * Unsigned that describes the direction of the flow.
     * We use the following convention to encode direction:
     * 0=N, 1=NW, 2=W, 3=SW, 4=S, 5=SE, 6=E, 7=NE.
     */
    unsigned mAdvectionDirection;

    /**
     * The speed of the flow.
     */
    double mAdvectionSpeed;

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
        archive & mAdvectionDirection;
        archive & mAdvectionSpeed;
    }

public:

    /**
     * Constructor.
     *
     * @param advectionDirection the direction of the flow (must take one of the values 0 to 7 inclusive)
     * @param advectionSpeed the speed of the flow
     */
    AdvectionCaUpdateRule(unsigned advectionDirection, double advectionSpeed);

    /**
     * Alternative constructor, for use in archiving.
     */
    AdvectionCaUpdateRule();

    /**
     * Destructor.
     */
    ~AdvectionCaUpdateRule();

    /**
     * Overridden GetNewLocationOfCell() method.
     *
     * This moves the cell in the prescribed direction mAdvectionDirection to its nearest neighbour
     * if this is free, with a probability that scales to ensure that the average speed is mAdvectionSpeed.
     *
     * @param currentLocationIndex the current location index of a cell
     * @param rCellPopulation reference to the cell population
     * @param dt timestep of the simulation to calculate probability of movement in current timestep
     */
    unsigned GetNewLocationOfCell(unsigned currentLocationIndex,
                                  CaBasedCellPopulation<DIM>& rCellPopulation,
                                  double dt);

    /**
     * @return mAdvectionDirection.
     */
    unsigned GetAdvectionDirection();

    /**
     * @return mAdvectionSpeed.
     */
    double GetAdvectionSpeed();

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AdvectionCaUpdateRule)

#endif /*ADVECTIONCAUPDATERULE_HPP_*/
