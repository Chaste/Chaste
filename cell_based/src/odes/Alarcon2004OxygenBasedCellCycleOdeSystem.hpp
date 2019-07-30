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

#ifndef _ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_
#define _ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>

#include "AbstractOdeSystem.hpp"

/**
 * Represents the Alarcon et al. (2004) system of ODEs (see ticket #461).
 * [doi:10.1016/j.jtbi.2004.04.016]
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

    /** The oxygen concentration (this affects the ODE system). */
    double mOxygenConcentration;

    /** Whether the cell associated with this cell cycle ODE system is labelled (this affects the ODE system). */
    bool mIsLabelled;

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
     * @param oxygenConcentration is a non-dimensional oxygen concentration value between 0 and 1
     * @param isLabelled whether the cell associated with this cell cycle ODE system is labelled (this affects the ODE system)
     * @param stateVariables optional initial conditions for state variables (only used in archiving)
     */
    Alarcon2004OxygenBasedCellCycleOdeSystem(double oxygenConcentration,
                                             bool isLabelled,
                                             std::vector<double> stateVariables=std::vector<double>());

    /**
     * Destructor.
     */
    ~Alarcon2004OxygenBasedCellCycleOdeSystem();

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
     *
     * @param time at which to calculate whether the stopping event has occurred
     * @param rY value of the solution vector used to evaluate the RHS.
     *
     * @return whether or not stopping conditions have been met
     */
    bool CalculateStoppingEvent(double time, const std::vector<double>& rY);

    /**
     * Set #mIsLabelled.
     *
     * @param isLabelled whether the cell associated with this cell cycle ODE system is labelled (this affects the ODE system)
     */
    void SetIsLabelled(bool isLabelled);

    /**
     * @return #mIsLabelled.
     */
    bool IsLabelled() const;

    /**
     * @return #mOxygenConcentration.
     */
    double GetOxygenConcentration() const;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Alarcon2004OxygenBasedCellCycleOdeSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an Alarcon2004OxygenBasedCellCycleOdeSystem.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Alarcon2004OxygenBasedCellCycleOdeSystem * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const double oxygen_concentration = t->GetOxygenConcentration();
    ar & oxygen_concentration;

    const bool is_labelled = t->IsLabelled();
    ar & is_labelled;

    const std::vector<double>& state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/**
 * De-serialize constructor parameters and initialise an Alarcon2004OxygenBasedCellCycleOdeSystem.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Alarcon2004OxygenBasedCellCycleOdeSystem * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    double oxygen_concentration;
    ar & oxygen_concentration;

    bool is_labelled;
    ar & is_labelled;

    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new(t)Alarcon2004OxygenBasedCellCycleOdeSystem(oxygen_concentration, is_labelled, state_variables);
}
}
} // namespace ...

#endif /*_ALARCON2004OXYGENBASEDCELLCYCLEODESYSTEM_HPP_*/
