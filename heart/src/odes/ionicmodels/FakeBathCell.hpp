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
     */
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);

    /**
     * There isn't really a cell here, so we override this method to do nothing.
     *
     * @param tStart  unused
     * @param tEnd  unused
     */
    void ComputeExceptVoltage(double tStart, double tEnd);
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
