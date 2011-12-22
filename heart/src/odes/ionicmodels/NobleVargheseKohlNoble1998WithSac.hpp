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

#ifndef _CML_noble_varghese_kohl_noble_1998_basic_with_sac_
#define _CML_noble_varghese_kohl_noble_1998_basic_with_sac_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <cassert>
#include "AbstractCardiacCell.hpp"
#include "Exception.hpp"
#include "AbstractStimulusFunction.hpp"
#include "OdeSystemInformation.hpp"


/**
 *  The Noble98 'Basic' model, but hand-altered to add a stretch activation channel ionic
 *  current, which is dependent on the stretch the cell is under (an additional member variable
 *  which can be set be the user/mechanics solver
 */
class CML_noble_varghese_kohl_noble_1998_basic_with_sac : public AbstractCardiacCell
{
    friend class boost::serialization::access;
    /**
     * Checkpointing.
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mStretch;
        archive & boost::serialization::base_object<AbstractCardiacCell>(*this);
    }
private:

    /** The stretch the cell is under - affects the stretch-activated-channel ionic current */
    double mStretch;

public:
    /**
     * Constructor.
     * @param pSolver  ODE solver to use
     * @param pIntracellularStimulus  intracellular stimulus to apply
     */
    CML_noble_varghese_kohl_noble_1998_basic_with_sac(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                                                      boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /** Destructor. */
    ~CML_noble_varghese_kohl_noble_1998_basic_with_sac();

    /** Compute the total ionic current
     * @param pStateVariables  state variable vector; if not given, use cell's internal state
     */
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);

    /**
     * Evaluate the RHS of this cell's ODE system.
     * @param var_environment__time  current simulation time
     * @param rY  state variable vector
     * @param rDY  will be filled in with the derivatives
     */
    void EvaluateYDerivatives(double var_environment__time,
                              const std::vector<double> &rY,
                              std::vector<double> &rDY);

    /**
     *  Set the stretch (overloaded)
     *  @param stretch stretch
     */
    void SetStretch(double stretch)
    {
        assert(stretch > 0.0);
        mStretch = stretch;
    }

    /**
     *  Get the stretch
     */
    double GetStretch()
    {
        return mStretch;
    }

    /**
     * Get the intracellular calcium concentration
     *
     * @return the intracellular calcium concentration
     */
    double GetIntracellularCalciumConcentration()
    {
        return mStateVariables[16];
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CML_noble_varghese_kohl_noble_1998_basic_with_sac)

namespace boost
{
    namespace serialization
    {
        /**
         * Avoid the need for a default constructor.
         * @param ar
         * @param t
         * @param fileVersion
         */
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const CML_noble_varghese_kohl_noble_1998_basic_with_sac * t, const unsigned int fileVersion)
        {
            const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
            const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
            ar << p_solver;
            ar << p_stimulus;
        }

        /**
         * Avoid the need for a default constructor.
         * @param ar
         * @param t
         * @param fileVersion
         */
        template<class Archive>
        inline void load_construct_data(
            Archive & ar, CML_noble_varghese_kohl_noble_1998_basic_with_sac * t, const unsigned int fileVersion)
        {
            boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
            ar >> p_solver;
            ar >> p_stimulus;
            ::new(t)CML_noble_varghese_kohl_noble_1998_basic_with_sac(p_solver, p_stimulus);
        }

    }

}

#endif
