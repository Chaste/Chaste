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


#ifndef _ODESOLUTION_HPP_
#define _ODESOLUTION_HPP_

#include <vector>
#include <string>
#include <cassert>
#include <boost/shared_ptr.hpp>

#include "AbstractOdeSystemInformation.hpp"
#include "AbstractParameterisedSystem.hpp"

#ifdef CHASTE_CVODE
// CVODE headers
#include <nvector/nvector_serial.h>
#endif // CHASTE_CVODE

class AbstractOdeSystem; // Avoid cyclic include issues

/**
 * A class that that stores the output data from solving a system of ODEs, and allows us to save it to file.
 */
class OdeSolution
{
private:
    /** Variable for the number of timesteps. */
    unsigned mNumberOfTimeSteps;

    /** A vector of times at each timestep. */
    std::vector<double> mTimes;

    /** Solutions for each variable at each timestep. */
    std::vector<std::vector<double> > mSolutions;

    /** Derived quantities at each timestep. */
    std::vector<std::vector<double> > mDerivedQuantities;

    /** Parameters - these are currently assumed constant across all time */
    std::vector<double> mParameters;

    /** The ODE solver used to create these results */
    std::string mSolverName;

    /**
     * Information about the concrete ODE system class.
     *
     * Used to get names and units into the file output.
     */
    boost::shared_ptr<const AbstractOdeSystemInformation> mpOdeSystemInformation;

public:
    /**
     * Public constructor - ensures data is empty to start with.
     */
    OdeSolution();

    /**
     * Get the number of timesteps.
     *
     * @return #mNumberOfTimeSteps
     */
    unsigned GetNumberOfTimeSteps() const;

    /**
     * Set the number of timesteps.
     *
     * @param numTimeSteps the number of timesteps to use
     */
    void SetNumberOfTimeSteps(unsigned numTimeSteps);

    /**
     * Set the ODE system information
     *
     * @param pOdeSystemInfo  ODE system information (used to get the names and units of variables).
     */
    void SetOdeSystemInformation(boost::shared_ptr<const AbstractOdeSystemInformation> pOdeSystemInfo);

    /**
     * @return the values of a state variable, parameter or derived quantity with a given index in
     * the ODE system at each output timestep.
     * The index is that given by AbstractOdeSystemInformation::GetAnyVariableIndex, which for state
     * variables (the most common case) is equal to their index within the state variable vector.
     *
     * @param index  the index of the variable in the system
     */
    std::vector<double> GetVariableAtIndex(unsigned index) const;

    /**
     * @return the values of a state variable, parameter or derived quantity with a given name in
     * the ODE system for each output timestep.
     *
     * @param rName  the name of the variable to extract
     */
    std::vector<double> GetAnyVariable(const std::string& rName) const;

    /**
     * @return the times at which the solution to the ODE system is stored.
     *
     * @return #mTimes.
     */
    std::vector<double>& rGetTimes();

    /**
     * @return the times at which the solution to the ODE system is stored.
     *
     * @return #mTimes.
     */
    const std::vector<double>& rGetTimes() const;

    /**
     * @return the values of the solution to the ODE system at each timestep.
     *
     * @return #mSolutions.
     */
    std::vector<std::vector<double> >& rGetSolutions();

    /**
     * @return the values of the solution to the ODE system at each timestep.
     *
     * @return #mSolutions.
     */
    const std::vector<std::vector<double> >& rGetSolutions() const;


    /** Set the ODE solver used to create these results
     *  @param solverName solver used
     */
    void SetSolverName(std::string solverName)
    {
        mSolverName = solverName;
    }

    /** @return the ODE solver used to create these results */
    std::string GetSolverName()
    {
        return mSolverName;
    }

    /**
     * Calculate the derived quantities and store them and the current parameters for printing/accessing.
     *
     * @param pOdeSystem  the ODE system which was solved to generate this solution object
     */
    template<typename VECTOR>
    void CalculateDerivedQuantitiesAndParameters(AbstractParameterisedSystem<VECTOR>* pOdeSystem);

    /**
     * @return the derived quantities for this ODE system at each timestep.
     *
     * @param pOdeSystem  the ODE system which was solved to generate this solution object
     * @return  A vector of vectors of derived quantities for each time step.
     */
    std::vector<std::vector<double> >& rGetDerivedQuantities(AbstractParameterisedSystem<std::vector<double> >* pOdeSystem);

#ifdef CHASTE_CVODE
    /**
     * @return the derived quantities for this ODE system at each timestep.
     *
     * @param pOdeSystem  the ODE system which was solved to generate this solution object
     * @return  A std::vector of vectors of derived quantities for each time step.
     */
    std::vector<std::vector<double> >& rGetDerivedQuantities(AbstractParameterisedSystem<N_Vector>* pOdeSystem);
#endif //CHASTE_CVODE

    /**
     * This method currently assumes that #mParameters is constant through time.
     * This may not be the case when using modifiers.
     *
     * @param pOdeSystem  The ODE system which was solved to generate this solution object.
     * @return  A vector of the current system parameters.
     */
    template<typename VECTOR>
    std::vector<double>& rGetParameters(AbstractParameterisedSystem<VECTOR>* pOdeSystem);

    /**
     * Write the data to a file.
     *
     * @param directoryName  the directory in which to write the data to file
     * @param baseResultsFilename  the name of the file in which to write the data
     * @param timeUnits  name of the units of time used
     * @param stepsPerRow  the solution to the ODE system is written to file every
     *                    this number of timesteps (defaults to 1)
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     * @param precision the precision with which to write the data (i.e. exactly
     *    how many digits to display after the decimal point).  Defaults to 8.
     *    Must be between 2 and 20 (inclusive).
     * @param includeDerivedQuantities  whether to include parameters and derived quantities in the output.
     */
    void WriteToFile(std::string directoryName,
                     std::string baseResultsFilename,
                     std::string timeUnits,
                     unsigned stepsPerRow=1,
                     bool cleanDirectory=true,
                     unsigned precision=8,
                     bool includeDerivedQuantities=false);
};


#endif //_ODESOLUTION_HPP_
