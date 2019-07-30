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
#ifndef DYNAMICALLY_LOADABLE_LR91_HPP_
#define DYNAMICALLY_LOADABLE_LR91_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include "AbstractDynamicallyLoadableEntity.hpp"
#include <vector>

/**
 * This class represents the Luo-Rudy 1991 system of equations,
 * with support for being compiled into a .so and loaded at run-time.
 */
class DynamicallyLoadableLr91 : public AbstractCardiacCell, public AbstractDynamicallyLoadableEntity
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
        archive & boost::serialization::base_object<AbstractCardiacCell>(*this);
        archive & boost::serialization::base_object<AbstractDynamicallyLoadableEntity>(*this);
    }

    /* Constants for the model */

    /** membrane capcaitance*/
    static const double membrane_C;
    /** Faraday constant*/
    static const double membrane_F;
    /** Universal gas constant*/
    static const double membrane_R;
    /** Temeperature*/
    static const double membrane_T;
    /** Reversal potentila for background current*/
    static const double background_current_E_b;
    /** Maximal conductance for background current*/
    static const double background_current_g_b;
    /** Maximal conductance for sodium current*/
    static const double fast_sodium_current_g_Na;
    /** Intracellular potassium concentration*/
    static const double ionic_concentrations_Ki;
    /** Extracellular potassium concentration*/
    static const double ionic_concentrations_Ko;
    /** Intracellular sodium concentration*/
    static const double ionic_concentrations_Nai;
    /** Extracellular sodium concentration*/
    static const double ionic_concentrations_Nao;
    /** Maximal conductance for plateau current*/
    static const double plateau_potassium_current_g_Kp;
    /** Permeability ratio Na/K for potassium currents*/
    static const double time_dependent_potassium_current_PR_NaK;

    /** Another parameter, which is a function of the above */
    double fast_sodium_current_E_Na;

    /**
     *  Range-checking on the current values of the state variables. Make sure
     *  all gating variables have are within zero and one, and all concentrations
     *  are positive
     */
    void VerifyStateVariables();

public:
    /**
     * Constructor
     *
     * @param pSolver is a pointer to the ODE solver
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    DynamicallyLoadableLr91(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                            boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~DynamicallyLoadableLr91();

    /**
     * Fill in a vector representing the RHS of the Luo-Rudy 1991 system
     * of Odes at each time step, y' = [y1' ... yn'].
     * Some ODE solver will call this function repeatedly to solve for y = [y1 ... yn].
     *
     * @param time  the current time, in milliseconds
     * @param rY  current values of the state variables
     * @param rDY  to be filled in with derivatives
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

    /**
     * Returns the ionic current
     *
     * @param pStateVariables  optional state at which to evaluate the current
     * @return the total ionic current
     */
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);

    /**
     * Get the intracellular calcium concentration
     *
     * @return the intracellular calcium concentration
     */
    double GetIntracellularCalciumConcentration();
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(DynamicallyLoadableLr91)

namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a DynamicallyLoadableLr91 instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const DynamicallyLoadableLr91 * t, const unsigned int file_version)
{
    const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
    const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
    ar << p_solver;
    ar << p_stimulus;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a DynamicallyLoadableLr91 instance (using existing constructor).
 *
 * NB this constructor allocates memory for the other member variables too.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, DynamicallyLoadableLr91 * t, const unsigned int file_version)
{

    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
    ar >> p_solver;
    ar >> p_stimulus;
    ::new(t)DynamicallyLoadableLr91(p_solver, p_stimulus);
}
}
} // namespace ...

#endif // DYNAMICALLY_LOADABLE_LR91_HPP_
