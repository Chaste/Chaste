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


#ifndef ABSTRACTALGEBRAICCONTRACTIONMODEL_HPP_
#define ABSTRACTALGEBRAICCONTRACTIONMODEL_HPP_


#include "AbstractContractionModel.hpp"

/**
 *  Contraction models that give the active tension as an algebraic function
 *  of [Ca] or voltage, stretch, stretch rate and time; with no ODEs to
 *  integrate.
 */
class AbstractAlgebraicContractionModel : public AbstractContractionModel
{
protected:
    /** The time at the next timestep. Set in RunDoNotUpdate() */
    double mTime;

public:
    /** Constructor does nothing */
    AbstractAlgebraicContractionModel()
     : AbstractContractionModel()
    {
        mTime = 0.0;
    }

    /** No ODE to run, so this does nothing except save the time (using the
     *  time at the next timestep)
     *  @param startTime start time
     *  @param endTime end time
     *  @param timestep timestep for integrating ODEs if there were any
     */
    void RunDoNotUpdate(double startTime, double endTime, double timestep)
    {
        mTime = endTime;
    }

    /** No ODE to run, so this does nothing except save the time (using the
     *  time at the next timestep)
     *  @param startTime start time
     *  @param endTime end time
     *  @param timestep timestep for integrating ODEs if there were any
     */
    void RunAndUpdate(double startTime, double endTime, double timestep)
    {
        mTime = endTime;
    }

    /**
     *  Same as GetActiveTension() for algebraic models (uses which stretch and
     *  and stretch rate has been passed in).
     */
    double GetNextActiveTension()
    {
        return GetActiveTension();
    }

    /** No ODE so does nothing. */
    void UpdateStateVariables()
    {
    }
};

#endif /*ABSTRACTALGEBRAICCONTRACTIONMODEL_HPP_*/
