/*

Copyright (C) University of Oxford, 2005-2011

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
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
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
     * @param currentLocationIndex reference to vector of forces on nodes
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
