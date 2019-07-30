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

#ifdef CHASTE_CVODE
#ifndef _ABSTRACTCVODECELL_HPP_
#define _ABSTRACTCVODECELL_HPP_

// Serialization headers
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/shared_ptr.hpp>

// Chaste headers
#include "AbstractOdeSystemInformation.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractIvpOdeSolver.hpp"
#include "AbstractCvodeSystem.hpp"
#include "AbstractCardiacCellInterface.hpp"
#include "VectorHelperFunctions.hpp"


/**
 * A cardiac cell that is designed to be simulated using CVODE.
 * It uses CVODE's vector type natively via AbstractCvodeSystem.
 *
 * Functionality is similar to that provided by AbstractCardiacCell and AbstractOdeSystem,
 * but not identical.  It also includes a direct interface to the CVODE solver, via the
 * Solve methods, since the CvodeAdaptor class may be a bit slower.
 *
 * Assumes that it will be solving stiff systems, so uses BDF/Newton.
 *
 * Various methods in this class just call methods on AbstractCvodeSystem, to
 * reduce compiler confusion when working with things in an inheritance tree.
 *
 * Any single cell work should generally use this class rather than AbstractCardiacCell.
 *
 * This class may also be faster for certain tissue problems (especially when using a
 * relatively large PDE timestep (~0.1ms)). BUT it will come with a memory overhead
 * as every node has to carry a CVODE solver object as it stores the internal
 * state of the solver, not just the state variables as ForwardEuler solver does.
 * So may not be appropriate for very large meshes.
 *
 */
class AbstractCvodeCell : public AbstractCardiacCellInterface, public AbstractCvodeSystem
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base classes.
        archive & boost::serialization::base_object<AbstractCvodeSystem>(*this);
        archive & boost::serialization::base_object<AbstractCardiacCellInterface>(*this);
        archive & mMaxDt;
    }

protected:
    /** The maximum timestep to use. */
    double mMaxDt;

public:
    /**
     * Create a new cardiac cell.
     *
     * @note subclasses @b must call Init() in their constructor after setting #mpSystemInfo.
     *
     * @param pSolver  not used for these cells (they're always solved with CVODE); may be empty
     * @param numberOfStateVariables  the size of the ODE system modelling this cell
     * @param voltageIndex  the index of the transmembrane potential within the vector of state variables
     * @param pIntracellularStimulus  the intracellular stimulus current
     */
    AbstractCvodeCell(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                      unsigned numberOfStateVariables,
                      unsigned voltageIndex,
                      boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Free the state variables, if they have been set.
     */
    virtual ~AbstractCvodeCell();

    /**
     * @return the current value of the transmembrane potential, as given
     * in our state variable vector.
     */
    double GetVoltage();

    /**
     * Set the transmembrane potential
     * @param voltage  new value
     */
    void SetVoltage(double voltage);

    /**
     * Set the maximum timestep to use for simulating this cell.
     *
     * @param maxDt  the maximum timestep
     */
    void SetMaxTimestep(double maxDt);

    /**
     * Set the \b maximum timestep to use for simulating this cell.
     *
     * As CVODE adaptively alters the timestep used, this method just sets the maximum timestep allowed (despite its name).
     * It is required as our base class AbstractCardiacCellInterface declares it as a pure virtual method.
     * Users using this class directly should call SetMaxTimestep instead of this method, for clearer code.
     *
     * @param maxDt  the maximum timestep
     */
    void SetTimestep(double maxDt);

    /**
     * @return The maximum timestep that is used by CVODE with this cell.
     */
    double GetTimestep();

    /**
     * Simulate this cell's behaviour between the time interval [tStart, tEnd],
     * updating the internal state variable values.
     *
     * The maximum time step to use is given by #mMaxDt, which defaults to
     * HeartConfig::Instance()->GetPrintingTimeStep() if unset.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    virtual void SolveAndUpdateState(double tStart, double tEnd);

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * and return state variable values.
     *
     * The maximum time step to use will be taken as #mMaxDt.  If this is unset
     * it is the same as tSamp, which defaults to
     * HeartConfig::Instance()->GetPrintingTimeStep().
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     * @param tSamp  sampling interval for returned results (defaults to HeartConfig printing time step)
     * @return solution object
     */
    OdeSolution Compute(double tStart, double tEnd, double tSamp=0.0);

    /**
     * Simulates this cell's behaviour between the time interval [tStart, tEnd],
     * but does not update the voltage.
     *
     * @param tStart  beginning of the time interval to simulate
     * @param tEnd  end of the time interval to simulate
     */
    void ComputeExceptVoltage(double tStart, double tEnd);

    /**
     * Set whether to clamp the voltage by setting its derivative to zero.
     * We need to ensure CVODE is re-initialised if this setting changes.
     * @param clamp  whether to clamp
     */
    void SetVoltageDerivativeToZero(bool clamp=true);

    /**
     * This just returns the number of state variables in the cell model.
     *
     * It is here because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface
     *
     * @return the number of state variables
     */
    unsigned GetNumberOfStateVariables() const;

    /**
     * This just returns the number of parameters in the cell model.
     *
     * It is here because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface
     *
     * @return the number of parameters
     */
    unsigned GetNumberOfParameters() const;

    /**
     * This just returns the state variables in the cell model.
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @return the state variables
     */
    std::vector<double> GetStdVecStateVariables();

    /**
     * Just calls AbstractCvodeSystem::rGetStateVariableNames().
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @return the state variable names in the cell's ODE system.
     */
    const std::vector<std::string>& rGetStateVariableNames() const;

    /**
     * This just sets the state variables in the cell model.
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @param rVariables  the state variables (to take a copy of).
     */
    void SetStateVariables(const std::vector<double>& rVariables);

    /**
     * This is also needed now just to show there is an alternative to the above method!
     * @param rVariables  the state variables (to take a copy of).
     */
    void SetStateVariables(const N_Vector& rVariables);

    /**
     * This just calls the method AbstractCvodeSystem::SetStateVariable
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @param index index of the state variable to be set
     * @param newValue new value of the state variable
     */
    void SetStateVariable(unsigned index, double newValue);

    /**
     * This just calls the method AbstractCvodeSystem::SetStateVariable
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @param rName name of the state variable to be set
     * @param newValue new value of the state variable
     */
    void SetStateVariable(const std::string& rName, double newValue);

    /**
     * This just calls the method AbstractCvodeSystem::GetAnyVariable
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @param rName variable name
     * @param time  the time at which to evaluate variable (only needed for derived quantities).
     * @return value of the variable at that time
     */
    double GetAnyVariable(const std::string& rName, double time=0.0);

    /**
     * This just calls AbstractCvodeSystem::GetParameter
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @param rParameterName  the name of a parameter to get the value of,
     * @return  the parameter's value.
     */
    double GetParameter(const std::string& rParameterName);

    /**
     * This is just here to show there is an alternative to the above!
     *
     * It just calls the base class method.
     *
     * @param parameterIndex  the index of the parameter vector entry to return
     * @return the parameter value
     */
    double GetParameter(unsigned parameterIndex);

    /**
     * This just calls AbstractCvodeSystem::SetParameter
     *
     * It is here (despite being inherited) because we seem to need to specify explicitly
     * which method in the parent classes we intend to implement
     * to take care of the pure definition in AbstractCardiacCellInterface.
     *
     * @param rParameterName  the parameter name to set the value of,
     * @param value  value to set it to.
     */
    void SetParameter(const std::string& rParameterName, double value);
};

CLASS_IS_ABSTRACT(AbstractCvodeCell)

#endif // _ABSTRACTCVODECELL_HPP_
#endif // CHASTE_CVODE
