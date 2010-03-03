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
#ifndef ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_
#define ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cfloat>

#include "AbstractOdeBasedCellCycleModel.hpp"
#include "Alarcon2004OxygenBasedCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "CellwiseData.hpp"
#include "Exception.hpp"

/**
 *  Oxygen-dependent cell cycle model.
 *
 * Note that this class uses C++'s default copying semantics, and so
 * doesn't implement a copy constructor or operator=.
 *
 * Note also that this model currently only works in 2D, since the
 * SolveOdeToTime() and GetDivideTime() methods involve instances of
 * CellwiseData<2>.
 */
class Alarcon2004OxygenBasedCellCycleModel : public AbstractOdeBasedCellCycleModel
{
    friend class boost::serialization::access;

private:

    /**
     * Fourth-order Runge Kutta solver for ODE system.
     */
    static RungeKutta4IvpOdeSolver msSolver;

    /**
     * The spatial dimension (needed by the templated class CellwiseData).
     */
    unsigned mDimension;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        assert(mpOdeSystem!=NULL);
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
        // Reference can be read or written into once mpOdeSystem has been set up
        // mpOdeSystem isn't set up by the first constructor, but is by the second
        // which is now utilised by the load_construct at the bottom of this file.
        archive & static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(mpOdeSystem)->rGetMutationState();
        archive & mDimension;
    }

public:

    /**
     * Default constructor, variables are set by abstract classes.
     *
     * @param dimension the spatial dimension (needed by the templated class CellwiseData)
     */
    Alarcon2004OxygenBasedCellCycleModel(unsigned dimension);

    /**
     * Copy constructor.
     *
     * Also copies our ODE system.
     *
     * @param rOtherModel the instance being copied.
     */
    Alarcon2004OxygenBasedCellCycleModel(const Alarcon2004OxygenBasedCellCycleModel& rOtherModel);



    /**
     * A private constructor for archiving.
     *
     * @param rParentProteinConcentrations a std::vector of doubles of the protein concentrations (see WntCellCycleOdeSystem)
     * @param rMutationState the mutation state of the cell (used by ODEs)
     * @param rDimension the spatial dimension
     */
    Alarcon2004OxygenBasedCellCycleModel(const std::vector<double>& rParentProteinConcentrations,
                                         const CryptCellMutationState& rMutationState,
                                         const unsigned& rDimension);

    /**
     * Resets the oxygen-based model to the start of the cell cycle
     * (this model does not cycle naturally). Cells are given a new
     * birth time and cell cycle proteins are reset. Note that the
     * oxygen concentration maintains its current value.
     *
     * Should only be called by the TissueCell Divide() method.
     */
    virtual void ResetForDivision();

    /**
     * Overridden builder method to create new copies of
     * this cell cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Initialise the cell cycle model at the start of a simulation.
     *
     * This overridden method sets up a new ODE system.
     */
    void Initialise();

    /**
     * Solve the ODEs up to the current time and return whether a stopping event occurred.
     *
     * @param currentTime the current time
     * @return whether a stopping event occured
     */
    bool SolveOdeToTime(double currentTime);

    /**
     * Get the time at which the ODE stopping event occured.
     *
     * @return the stopping event time
     */
    double GetOdeStopTime();

    /**
     * Get the spatial dimension.
     *
     * @return mDimension
     */
    unsigned GetDimension();

};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Alarcon2004OxygenBasedCellCycleModel)


namespace boost
{
namespace serialization
{
/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a Alarcon2004OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Alarcon2004OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a Alarcon2004OxygenBasedCellCycleModel instance.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Alarcon2004OxygenBasedCellCycleModel * t, const unsigned int file_version)
{
    /**
     * Invoke inplace constructor to initialise an instance of Alarcon2004OxygenBasedCellCycleModel.
     * It doesn't actually matter what values we pass to our standard constructor,
     * provided they are valid parameter values, since the state loaded later
     * from the archive will overwrite their effect in this case.
     */

    std::vector<double> state_vars;
    for (unsigned i=0; i<6; i++)
    {
        state_vars.push_back(0.0);
    }
    CryptCellMutationState mutation_state = HEALTHY;
    unsigned dimension = UINT_MAX;

    ::new(t)Alarcon2004OxygenBasedCellCycleModel(state_vars, mutation_state, dimension);
}
}
} // namespace ...

#endif /*ALARCON2004OXYGENBASEDCELLCYCLEMODEL_HPP_*/
