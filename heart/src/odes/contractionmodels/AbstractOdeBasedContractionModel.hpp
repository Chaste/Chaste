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
        if (mTemporaryStateVariables.size() > 0)
        {
            mTemporaryStateVariables.resize(mStateVariables.size());
        }
        mTemporaryStateVariables = mStateVariables;

        // solve (RunAndUpdate just calls SolveAndUpdateStateVariable() using an euler solver)
        RunAndUpdate(startTime, endTime, timeStep);

        // put the solution in mTemporaryStateVariables and return the state variables to its
        // original state
        for (unsigned i=0; i<mStateVariables.size(); i++)
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
        for (unsigned i=0; i<mStateVariables.size(); i++)
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
