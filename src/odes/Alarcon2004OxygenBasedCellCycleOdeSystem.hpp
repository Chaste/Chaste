/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef _ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_
#define _ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_

#include <cmath>

#include "AbstractOdeSystem.hpp"
#include "CryptCellMutationStates.hpp"

/**
 * Represents the Alarcon et al. (2004) system of ODEs (see ticket #461).
 * [doi:10.1016/j.jtbi.2004.04.016]
 *
 * Here the cell mutation state HEALTHY corresponds to a 'normal' state, while
 * LABELLED corresponds to a 'cancer' state.
 *
 * The variables are
 *
 *  0. x = Cdh1-APC complexes
 *  1. y = cyclin-CDK
 *  2. z = p27
 *  3. m = mass
 *  4. u = RBNP
 *  5. P = oxygen concentration
 */
class Alarcon2004OxygenBasedCellCycleOdeSystem : public AbstractOdeSystem
{
private:

    /**
     * Constants for the Alarcon et al. (2004) model
     */

    /** Dimensionless parameter a_1. */
    double ma1;
    /** Dimensionless parameter a_2. */
    double ma2;
    /** Dimensionless parameter a_3. */
    double ma3;
    /** Dimensionless parameter a_4. */
    double ma4;
    /** Dimensionless parameter b_3. */
    double mb3;
    /** Dimensionless parameter b_4. */
    double mb4;
    /** Dimensionless parameter c_1. */
    double mc1;
    /** Dimensionless parameter c_2. */
    double mc2;
    /** Dimensionless parameter d_1. */
    double md1;
    /** Dimensionless parameter d_2. */
    double md2;
    /** Dimensionless parameter J_3. */
    double mJ3;
    /** Dimensionless parameter J_4. */
    double mJ4;
    /** Dimensionless parameter eta. */
    double mEta;
    /** Dimensionless parameter m_star. */
    double mMstar;
    /** Dimensionless parameter B. */
    double mB;
    /** Dimensionless parameter x_THR. */
    double mxThreshold;
    /** Dimensionless parameter y_THR. */
    double myThreshold;

    /** The mutation state of the cell associated with this cell cycle model. */
    CryptCellMutationState mMutationState;

public:

    /**
     * Constructor.
     *
     * @param oxygenConcentration is a non-dimensional oxygen concentration value between 0 and 1.
     * @param rMutationState affects the ODE system and is given by CryptCellMutationStates.hpp
     */
    Alarcon2004OxygenBasedCellCycleOdeSystem(double oxygenConcentration, const CryptCellMutationState& rMutationState);

    /**
     * Destructor.
     */
    ~Alarcon2004OxygenBasedCellCycleOdeSystem();

    /**
     * Initialise parameter values.
     */
    void Init();

    /**
     * Set the mutation state of the cell.
     *
     * This should be called by the relevant cell cycle model before any solving
     * of the ODE system (as it is used to evaluate the Y derivatives).
     *
     * @param rMutationState the mutation state
     */
    void SetMutationState(const CryptCellMutationState& rMutationState);

    /**
     * Called by the archive function on the cell cycle model.
     *
     * @return mMutationState the mutation state of the cell defined by
     * CryptCellMutationStates.hpp
     */
    CryptCellMutationState& rGetMutationState();

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
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS.
     *
     * @return whether or not stopping conditions have been met
     */
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

};

#endif /*_ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_*/
