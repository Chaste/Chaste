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

#ifndef MIRAMS2010WNTODESYSTEM_HPP_
#define MIRAMS2010WNTODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"
#include "AbstractCellMutationState.hpp"

// Needed here to avoid serialization errors
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellLabel.hpp"
#include "WildTypeCellMutationState.hpp"

/**
 * Represents the Mirams et al. system of ODEs, based on Swat et al. (2004)
 * [doi:10.1093/bioinformatics/bth110]
 * and a simple Wnt model (unpublished)
 *
 * The variables are
 *
 * 6. b1 = Beta-Catenin (from 1st allele)
 * 7. b2 = Beta-Catenin (from 1st allele)
 * 8. WntLevel
 */
class Mirams2010WntOdeSystem : public AbstractOdeSystem
{
private:

    /**
     * Parameters for the Mirams et al. (2010) model
     */

    /** Dimensional parameter a. */
    double mA;
    /** Dimensional parameter b. */
    double mB;
    /** Dimensional parameter c. */
    double mC;
    /** Dimensional parameter d. */
    double mD;
    /** Dimensional parameter e. */
    double mE;
    /** Dimensional parameter f. */
    double mF;

    /** The mutation state of the cell (this affects the ODE system). */
    boost::shared_ptr<AbstractCellMutationState> mpMutationState;

    /** The Wnt level (this affects the ODE system). */
    double mWntLevel;

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
     * @param wntLevel is a non-dimensional Wnt value between 0 and 1. This sets up the Wnt pathway in its steady state.
     * @param pMutationState optional mutation state (affects the ODE system)
     * @param stateVariables optional initial conditions for state variables (only used in archiving)
     */
    Mirams2010WntOdeSystem(double wntLevel=0.0,
                           boost::shared_ptr<AbstractCellMutationState> pMutationState=boost::shared_ptr<AbstractCellMutationState>(),
                           std::vector<double> stateVariables=std::vector<double>());

    /**
     * Destructor.
     */
    ~Mirams2010WntOdeSystem();

    /**
     * Initialise parameter values.
     */
    void Init();

    /**
     * Set the mutation state of the cell.
     *
     * This should be called by the relevant cell-cycle model before any solving
     * of the ODE system (as it is used to evaluate the Y derivatives).
     *
     * @param pMutationState the mutation state.
     */
    void SetMutationState(boost::shared_ptr<AbstractCellMutationState> pMutationState);

    /**
     * Called by the archive function on the Wnt cell-cycle model.
     *
     * @return #mpMutationState
     */
    const boost::shared_ptr<AbstractCellMutationState> GetMutationState() const;

    /**
     * Compute the RHS of the WntCellCycle system of ODEs.
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
     * Get method for mWntLevel.
     */
    double GetWntLevel() const;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Mirams2010WntOdeSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a Mirams2010WntOdeSystem.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Mirams2010WntOdeSystem * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const double wnt_level = t->GetWntLevel();
    ar & wnt_level;

    const boost::shared_ptr<AbstractCellMutationState> p_mutation_state = t->GetMutationState();
    ar & p_mutation_state;

    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/**
 * De-serialize constructor parameters and initialise a Mirams2010WntOdeSystem.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Mirams2010WntOdeSystem * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    double wnt_level;
    ar & wnt_level;

    boost::shared_ptr<AbstractCellMutationState> p_mutation_state;
    ar & p_mutation_state;

    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new(t)Mirams2010WntOdeSystem(wnt_level, p_mutation_state, state_variables);
}
}
} // namespace ...

#endif /*MIRAMS2010WNTODESYSTEM_HPP_*/
