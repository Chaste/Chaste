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

#ifndef DELTANOTCHODESYSTEM_HPP_
#define DELTANOTCHODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"

/**
 * Represents the Delta-Notch ODE system described by Collier et al,
 * "Pattern formation by lateral inhibition with feedback: a mathematical
 * model of delta-notch intercellular signalling" (Journal of Theoretical
 * Biology 183:429-446, 1996).
 */
class DeltaNotchOdeSystem : public AbstractOdeSystem
{
private:

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

    ///\todo extract model parameters as member variables

public:

    /**
     * Default constructor.
     * 
     * @param meanDelta the average of the levels of Delta in the surrounding cells (defaults to zero)
     * @param stateVariables optional initial conditions for state variables (only used in archiving)
     */
    DeltaNotchOdeSystem(double meanDelta=0.0, std::vector<double> stateVariables=std::vector<double>());

    /**
     * Destructor.
     */
    ~DeltaNotchOdeSystem();
    
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
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchOdeSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a DeltaNotchOdeSystem.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const DeltaNotchOdeSystem * t, const BOOST_PFTO unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/**
 * De-serialize constructor parameters and initialise a DeltaNotchOdeSystem.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, DeltaNotchOdeSystem * t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new(t)DeltaNotchOdeSystem(0.0, state_variables);
}
}
} // namespace ...

#endif /*DELTANOTCHODESYSTEM_HPP_*/
