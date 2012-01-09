/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef BACKWARDEULERIVPODESOLVER_HPP_
#define BACKWARDEULERIVPODESOLVER_HPP_

#include "AbstractOneStepIvpOdeSolver.hpp"
#include "AbstractOdeSystemWithAnalyticJacobian.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractOneStepIvpOdeSolver.hpp"

/**
 * A concrete one step ODE solver class that employs the backward Euler
 * method. This numerical method is implicit and hence unconditionally stable.
 */
class BackwardEulerIvpOdeSolver  : public AbstractOneStepIvpOdeSolver
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the abstract IVP Solver, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // This calls serialize on the base class.
        archive & boost::serialization::base_object<AbstractOneStepIvpOdeSolver>(*this);
        //archive & mSizeOfOdeSystem; - this done in save and load construct now.
        archive & mNumericalJacobianEpsilon;
        archive & mForceUseOfNumericalJacobian;
    }

    /** The number of state variables in the ODE system. */
    unsigned mSizeOfOdeSystem;

    /** The epsilon to use in calculating the numerical Jacobian of the ODE system. */
    double mNumericalJacobianEpsilon;

    /**
     * Whether to force the solver to use the numerical Jacobian even if
     * the ODE system object provides an analytical Jacobian.
     */
    bool mForceUseOfNumericalJacobian;

    /*
     * NOTE: we use (unsafe) double pointers here rather than
     * std::vectors because using std::vectors would lead to a
     * slow down by a factor of about 4.
     */

    /** Working memory : residual vector */
    double* mResidual;

    /** Working memory : Jacobian matrix */
    double** mJacobian;

    /** Working memory : update vector */
    double* mUpdate;

    /**
     * Compute the current residual.
     *
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rCurrentGuess  current guess for the state at the next timestep
     */
    void ComputeResidual(AbstractOdeSystem* pAbstractOdeSystem,
                         double timeStep,
                         double time,
                         std::vector<double>& rCurrentYValues,
                         std::vector<double>& rCurrentGuess);

    /**
     * Compute the Jacobian of the ODE system.
     *
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rCurrentGuess  current guess for the state at the next timestep
     */
    void ComputeJacobian(AbstractOdeSystem* pAbstractOdeSystem,
                         double timeStep,
                         double time,
                         std::vector<double>& rCurrentYValues,
                         std::vector<double>& rCurrentGuess);

    /**
     * Solve a linear system of equations to update the
     * current guess for the solution to the ODE system at
     * the next timestep.
     * Used by the method CalculateNextYValue.
     */
    void SolveLinearSystem();

    /**
     * Compute the infinity/maximum norm of a vector.
     * Used by the method CalculateNextYValue.
     *
     * @param pVector  a pointer to a vector
     * @return the vector's norm.
     */
    double ComputeNorm(double* pVector);

    /**
     * Compute the Jacobian of the ODE system numerically.
     *
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rCurrentGuess  current guess for the state at the next timestep
     */
    void ComputeNumericalJacobian(AbstractOdeSystem* pAbstractOdeSystem,
                                  double timeStep,
                                  double time,
                                  std::vector<double>& rCurrentYValues,
                                  std::vector<double>& rCurrentGuess);

protected:

    /**
     * Calculate the solution to the ODE system at the next timestep.
     *
     * A usage example:
     *     BackwardEulerIvpOdeSolver mySolver;
     *     OdeSolution solution = mySolver.Solve(pMyOdeSystem, yInit, StartTime, EndTime, TimeStep, SamplingTime);
     *
     * @param pAbstractOdeSystem  the ODE system to solve
     * @param timeStep  dt
     * @param time  the current time
     * @param rCurrentYValues  the current (initial) state
     * @param rNextYValues  the state at the next timestep
     */
    void CalculateNextYValue(AbstractOdeSystem* pAbstractOdeSystem,
                             double timeStep,
                             double time,
                             std::vector<double>& rCurrentYValues,
                             std::vector<double>& rNextYValues);

public:

    /**
     * Constructor.
     *
     * @param sizeOfOdeSystem  the number of state variables in the ODE system
     */
    BackwardEulerIvpOdeSolver(unsigned sizeOfOdeSystem);

    /**
     * Destructor.
     */
    ~BackwardEulerIvpOdeSolver();

    /**
     * Set the epsilon to use in calculating the
     * numerical Jacobian of the ODE system.
     *
     * @param epsilon
     */
    void SetEpsilonForNumericalJacobian(double epsilon);

    /**
     * Force the solver to use the numerical Jacobian even if
     * the ODE system object provides an analytical Jacobian.
     */
    void ForceUseOfNumericalJacobian();

    /**
     * Public method used in archiving.
     *
     * @return the size of the system
     */
     unsigned GetSystemSize() const {return mSizeOfOdeSystem;};
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BackwardEulerIvpOdeSolver)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a BackwardEulerIvpOdeSolver instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const BackwardEulerIvpOdeSolver * t, const unsigned int file_version)
{
    const unsigned system_size = t->GetSystemSize();
    ar & system_size;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a BackwardEulerIvpOdeSolver instance (using existing constructor)
 *
 * NB this constructor allocates memory for the other member variables too.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, BackwardEulerIvpOdeSolver * t, const unsigned int file_version)
{
     unsigned ode_system_size;
     ar >> ode_system_size;
     ::new(t)BackwardEulerIvpOdeSolver(ode_system_size);
}
}
} // namespace ...

#endif /*BACKWARDEULERIVPODESOLVER_HPP_*/
