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

#ifdef CHASTE_CVODE
#ifndef _CVODEADAPTOR_HPP_
#define _CVODEADAPTOR_HPP_

#include <vector>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractIvpOdeSolver.hpp"
#include "OdeSolution.hpp"

// CVODE headers
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>


/**
 * CVODE error handling function.
 *
 * Throw an Exception to report errors, rather than the CVODE approach of magic
 * return codes.
 */
void CvodeErrorHandler(int errorCode, const char *module, const char *function,
                       char *message, void* pData);
// Note: declared here since it's also used by AbstractCvodeCell.


/**
 * Data structure passed to CVODE calls, allowing our callback functions
 * to access the Chaste objects.
 */
typedef struct CvodeData_ {
    /** Working memory. */
    std::vector<realtype>* pY;
    /** The ODE system being solved. */
    AbstractOdeSystem* pSystem;
} CvodeData;

/**
 * The CVODE adaptor ODE solver class.
 *
 * Assumes that it will be solving stiff systems, so uses BDF/Newton.
 *
 * The timeStep parameters of the abstract class are here used to specify
 * *maximum* steps, since the solver is adaptive.
 *
 * Note that a call to Solve will initialise the CVODE solver, and free its
 * working memory when done.  There is thus a non-trivial overhead involved.
 *
 * \todo Add an option to just initialise once, and assume subsequent Solve
 *   calls are continuing from where we left off.
 */
class CvodeAdaptor : public AbstractIvpOdeSolver
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractIvpOdeSolver>(*this);
        archive & mRelTol;
        archive & mAbsTol;
        archive & mLastInternalStepSize;
        archive & mMaxSteps;
        archive & mCheckForRoots;
        // All other member variables given values on each call.
    }

    /** Pointer to the CVODE memory block. */
    void* mpCvodeMem;

    /** Initial conditions for the ODE solver. */
    N_Vector mInitialValues;

    /** The CVODE data structure. */
    CvodeData mData;

    /** Relative tolerance for the ODE solver. */
    double mRelTol;

    /** Absolute tolerance for the ODE solver. */
    double mAbsTol;

    /** The size of the previous timestep. */
    double mLastInternalStepSize;

    /**
     * The maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.
     */
    long int mMaxSteps;

    /** Whether to check for stopping events. */
    bool mCheckForRoots;

protected:

    /**
     * Set up the CVODE data structures needed to solve the given system.
     *
     * @param pOdeSystem  the ODE system being solved
     * @param rInitialY  initial conditions vector
     * @param startTime  when to simulate from
     * @param maxStep  maximum time step
     */
    void SetupCvode(AbstractOdeSystem* pOdeSystem,
                    std::vector<double>& rInitialY,
                    double startTime, double maxStep);

    /**
     * Free CVODE memory after solving a system of ODEs.
     */
    void FreeCvodeMemory();

    /**
     * Report an error from CVODE.
     *
     * This will (probably) never be called, since we supply an error handler function
     * which throws an exception.
     *
     * @param flag  error flag
     * @param msg  error message
     */
    void CvodeError(int flag, const char * msg);

public:

    /**
     * Default constructor.
     * Can optionally set relative and absolute tolerances.
     *
     * @param relTol the relative tolerance for the solver
     * @param absTol the absolute tolerance for the solver
     */
    CvodeAdaptor(double relTol=1e-4, double absTol=1e-6);

    /**
     * Set relative and absolute tolerances; both scalars.
     * If no parameters are given, tolerances will be reset to default values.
     *
     * @param relTol the relative tolerance for the solver
     * @param absTol the absolute tolerance for the solver
     */
    void SetTolerances(double relTol=1e-4, double absTol=1e-6);

    /**
     * Get the relative tolerance.
     */
    double GetRelativeTolerance();

    /**
     * Get the absolute tolerance.
     */
    double GetAbsoluteTolerance();

    /**
     * Get the last step size used internally by CVODE in the last Solve call
     */
    double GetLastStepSize();

    /**
     * Solve the given ODE system, returning the solution at sampling intervals.
     *
     * @param pOdeSystem  the ODE system to solve
     * @param rYValues  the initial state variable values
     *   (note: this vector will also be used as working memory)
     * @param startTime  the time to start solving at
     * @param endTime  the time to solve to
     * @param maxStep  the maximum time step to be taken by the adaptive solver
     * @param timeSampling  the interval at which to sample the solution
     * @return  the solution
     */
    OdeSolution Solve(AbstractOdeSystem* pOdeSystem,
                      std::vector<double>& rYValues,
                      double startTime,
                      double endTime,
                      double maxStep,
                      double timeSampling);

    /**
     * Solve the given ODE system, storing the final result in rYValues.
     *
     * @param pOdeSystem  the ODE system to solve
     * @param rYValues  the initial state variable values; will be filled in with
     *   the final values on return
     * @param startTime  the time to start solving at
     * @param endTime  the time to solve to
     * @param maxStep  the maximum time step to be taken by the adaptive solver
     */
    void Solve(AbstractOdeSystem* pOdeSystem,
               std::vector<double>& rYValues,
               double startTime,
               double endTime,
               double maxStep);

    /**
     * Make the solver check for stopping events using CVODE's rootfinding functionality.
     *
     * By default we do not check.
     */
    void CheckForStoppingEvents();

    /**
     * Change the maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.  Default is 500.
     *
     * @param numSteps  new maximum number of steps
     */
    void SetMaxSteps(long int numSteps);

    /**
     * Get the maximum number of steps to be taken by the solver
     * in its attempt to reach the next output time.
     */
    long int GetMaxSteps();

};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CvodeAdaptor)

#endif // _CVODEADAPTOR_HPP_
#endif // CHASTE_CVODE
