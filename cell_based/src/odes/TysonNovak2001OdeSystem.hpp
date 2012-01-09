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

#ifndef TYSONNOVAK2001ODESYSTEM_HPP_
#define TYSONNOVAK2001ODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystemWithAnalyticJacobian.hpp"

/**
 * Represents the Tyson & Novak (2001) system of ODEs.
 * [doi:10.1006/jtbi.2001.2293]
 */
class TysonNovak2001OdeSystem : public AbstractOdeSystemWithAnalyticJacobian
{
private:

    /**
     * Parameters for the Tyson & Novak (2001) model.
     */

    /** Dimensional parameter k_1. */
    double mK1;
    /** Dimensional parameter k_2'. */
    double mK2d;
    /** Dimensional parameter k_2''. */
    double mK2dd;
    /** Dimensional parameter k_2'''. */
    double mK2ddd;
    /** Dimensionless parameter [CycB]_threshold. */
    double mCycB_threshold;
    /** Dimensional parameter k_3'. */
    double mK3d;
    /** Dimensional parameter k_3''. */
    double mK3dd;
    /** Dimensional parameter k_4'. */
    double mK4d;
    /** Dimensional parameter k_4. */
    double mK4;
    /** Dimensionless parameter J_3. */
    double mJ3;
    /** Dimensionless parameter J_4. */
    double mJ4;
    /** Dimensional parameter k_5'. */
    double mK5d;
    /** Dimensional parameter k_5''. */
    double mK5dd;
    /** Dimensional parameter k_6. */
    double mK6;
    /** Dimensionless parameter J_5. */
    double mJ5;
    /** Dimensionless parameter n. */
    unsigned mN;
    /** Dimensional parameter k_7. */
    double mK7;
    /** Dimensional parameter k_8. */
    double mK8;
    /** Dimensionless parameter J_7. */
    double mJ7;
    /** Dimensionless parameter J_8. */
    double mJ8;
    /** Dimensionless parameter [Mad]. */
    double mMad;
    /** Dimensional parameter k_9. */
    double mK9;
    /** Dimensional parameter k_10. */
    double mK10;
    /** Dimensional parameter k_11. */
    double mK11;
    /** Dimensional parameter k_12'. */
    double mK12d;
    /** Dimensional parameter k_12''. */
    double mK12dd;
    /** Dimensional parameter k_12'''. */
    double mK12ddd;
    /** Dimensionless parameter K_eq. */
    double mKeq;
    /** Dimensional parameter k_13. */
    double mK13;
    /** Dimensional parameter k_14. */
    double mK14;
    /** Dimensional parameter k_15'. */
    double mK15d;
    /** Dimensional parameter k_15''. */
    double mK15dd;
    /** Dimensional parameter k_16'. */
    double mK16d;
    /** Dimensional parameter k_16''. */
    double mK16dd;
    /** Dimensionless parameter J_15. */
    double mJ15;
    /** Dimensionless parameter J_16. */
    double mJ16;
    /** Dimensional parameter mu. */
    double mMu;
    /** Dimensionless parameter m_star. */
    double mMstar;

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
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param stateVariables optional initial conditions for state variables (only used in archiving)
     */
    TysonNovak2001OdeSystem(std::vector<double> stateVariables=std::vector<double>());

    /**
     * Destructor.
     */
    ~TysonNovak2001OdeSystem();

    /**
     * Initialise parameter values.
     */
    void Init();

    /**
     * Compute the RHS of the Alarcon et al. (2004) system of ODEs.
     *
     * Returns a vector representing the RHS of the ODEs at each time step, y' = [y1' ... yn'].
     * An ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives (using Alarcons et al. (2004) system of equations).
     */
    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

    /**
     * Calculate whether the conditions for the cell cycle to finish have been met.
     * (Used by Chaste solvers to find whether or not to stop solving)
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS.
     *
     * @return whether or not stopping conditions have been met
     */
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    /**
     * Calculate whether the conditions for the cell cycle to finish have been met.
     * (Used by CVODE solver to find exact stopping position)
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS.
     *
     * @return How close we are to the root of the stopping condition
     */
    double CalculateRootFunction(double time, const std::vector<double>& rY);

    /**
     * Compute the Jacobian of the ODE system.
     *
     * @param rSolutionGuess initial guess for the solution vector.
     * @param jacobian the Jacobian of the ODE system.
     * @param time at which to calculate the Jacobian.
     * @param timeStep used to calculate the Jacobian.
     */
    virtual void AnalyticJacobian(const std::vector<double>& rSolutionGuess, double** jacobian, double time, double timeStep);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(TysonNovak2001OdeSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TysonNovak2001OdeSystem.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const TysonNovak2001OdeSystem * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/**
 * De-serialize constructor parameters and initialise a TysonNovak2001OdeSystem.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, TysonNovak2001OdeSystem * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new(t)TysonNovak2001OdeSystem(state_variables);
}
}
} // namespace ...

#endif /*TYSONNOVAK2001ODESYSTEM_HPP_*/
