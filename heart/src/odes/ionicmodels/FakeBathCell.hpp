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
#ifndef FAKEBATHCELL_HPP_
#define FAKEBATHCELL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * This class represents a fake cell for use within the bath of a bidomain simulation.
 *
 * \note Note that only a portion of the normal functionality of a cardiac cell is
 * actually redefined in this class.  If further calls to cardiac cells are later
 * added to the simulation process, additional overrides may need to be added here.
 */
class FakeBathCell : public AbstractCardiacCell
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive this cell.  Just calls the base class version.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCardiacCell>(*this);
        // In case we're loading an old archive, update mNumberOfStateVariables
        this->mNumberOfStateVariables = 1;
    }

public:
    /**
     * Constructor uses the same signature as normal cells, for convenience.
     *
     * @param pSolver  unused
     * @param pIntracellularStimulus  unused
     */
    FakeBathCell(boost::shared_ptr<AbstractIvpOdeSolver> pSolver,
                 boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor; does nothing.
     */
    ~FakeBathCell();

    /**
     * This method is pure in a base class, so we need it, but we never use it.
     * It has an empty body.
     *
     * @param time  unused
     * @param rY  unused
     * @param rDY  unused
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double> &rDY);

    /**
     * Fake cells have no transmembrane currents, so this method always returns 0.
     *
     * @param pStateVariables  unused
     * @return zero
     */
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);

    /**
     * There isn't really a cell here, so we override this method to do nothing.
     *
     * @param tStart  unused
     * @param tEnd  unused
     */
    void ComputeExceptVoltage(double tStart, double tEnd);

    /**
     * There is really no calcium here, so we override this method to return a dummy value (0)
     * Implementing this with a  dummy implementation is needed by mechanics
     *
     * @return always zero
     */
    double GetIntracellularCalciumConcentration();
};


#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(FakeBathCell)

namespace boost
{
namespace serialization
{
/**
 * Save the data needed to create a FakeBathCell.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const FakeBathCell * t, const unsigned int file_version)
{
    const boost::shared_ptr<AbstractIvpOdeSolver> p_solver = t->GetSolver();
    const boost::shared_ptr<AbstractStimulusFunction> p_stimulus = t->GetStimulusFunction();
    ar << p_solver;
    ar << p_stimulus;
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a FakeBathCell instance (using existing constructor)
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, FakeBathCell * t, const unsigned int file_version)
{

    boost::shared_ptr<AbstractIvpOdeSolver> p_solver;
    boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
    ar >> p_solver;
    ar >> p_stimulus;
    ::new(t)FakeBathCell(p_solver, p_stimulus);
}
}
} // namespace ...


#endif /*FAKEBATHCELL_HPP_*/
