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

#ifndef ABSTRACTODEBASEDCONTRACTIONMODEL_
#define ABSTRACTODEBASEDCONTRACTIONMODEL_

#include "AbstractOdeSystem.hpp"
#include "AbstractContractionModel.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"

/**
 *  Abstract base class for ODE-based contraction models. Inherits from AbstractOdeSystem
 *  and AbstractContractionModel and deals with the ODE solving.
 *
 *  Usage is either:
 *   // for the explicit electromechanics algorithm
 *   RunAndUpdate(); // solves and updates state variables to the solution
 *
 *  or
 *   // for the implicit electromechanics algorithm
 *   RunDoNotUpdate();
 *   // ..
 *   UpdateStateVariables(); // if keeping this solution
 */
class AbstractOdeBasedContractionModel : public AbstractOdeSystem, public AbstractContractionModel
{
protected:

    /** A second vector of state variables, where the results will go
     *  when RunDoNotUpdate() is called */
    std::vector<double> mTemporaryStateVariables;

    /** The time (at the next timestep) to be used in GetActiveTension if required */
    double mTime;

public:
    /**
     *  Constructor
     *  @param numStateVariables number of state variables
     */
    AbstractOdeBasedContractionModel(unsigned numStateVariables)
        : AbstractOdeSystem(numStateVariables),
          AbstractContractionModel(),
          mTime(0.0)
    {
        // note mTemporaryStateVariables not resized here - only resized if needed
    }


    /**
     *  Solves the ODEs, but doesn't update the state variables, instead keeps them in
     *  a temporary store. Call UpdateStateVariables() to save the new values. Call
     *  GetNextActiveTension() to get the active tension corresponding to the new values (if
     *  UpdateStateVariables() has not been called). Also saves the time (using endTime).
     *
     *  @param startTime start time
     *  @param endTime end time
     *  @param timeStep timestep for integrating ODEs
     */
    virtual void RunDoNotUpdate(double startTime, double endTime, double timeStep)
    {
        // save the state variables
        if(mTemporaryStateVariables.size()>0)
        {
            mTemporaryStateVariables.resize(mStateVariables.size());
        }
        mTemporaryStateVariables = mStateVariables;

        // solve (RunAndUpdate just calls SolveAndUpdateStateVariable() using an euler solver)
        RunAndUpdate(startTime, endTime, timeStep);

        // put the solution in mTemporaryStateVariables and return the state variables to its
        // original state
        for(unsigned i=0; i<mStateVariables.size(); i++)
        {
            double soln = mStateVariables[i];
            mStateVariables[i] = mTemporaryStateVariables[i];
            mTemporaryStateVariables[i] = soln;
        }

        // note that the end time was saved in RunAndUpdate(): mTime = endTime;
    }

    /**
     *  After RunDoNotUpdate() has been called, this call be used to update the state
     *  variables to the new (saved) values
     */
    void UpdateStateVariables()
    {
        // save the state variables
        for(unsigned i=0; i<mStateVariables.size(); i++)
        {
            mStateVariables[i] = mTemporaryStateVariables[i];
        }
    }

    /**
     *  Solves the ODEs and updates the state variable to the new solution
     *
     *  Alternative usage:
     *   RunDoNotUpdate();
     *   // ..
     *   UpdateStateVariables(); // if keeping this solution
     *
     *  @param startTime start time
     *  @param endTime end time
     *  @param timeStep timestep for integrating ODEs
     */
    void RunAndUpdate(double startTime, double endTime, double timeStep)
    {
        EulerIvpOdeSolver solver;
        solver.SolveAndUpdateStateVariable(this, startTime, endTime, timeStep);

        mTime = endTime;
    }
};


#endif /*ABSTRACTODEBASEDSTRETCHINDEPENDENTCONTRACTIONMODEL_*/
