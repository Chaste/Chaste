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

#ifndef TYSONNOVAKCELLCYCLEMODEL_HPP_
#define TYSONNOVAKCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractOdeBasedCellCycleModel.hpp"
#include "TysonNovak2001OdeSystem.hpp"


/**
 *  Tyson-Novak 2001 cell-cycle model, taken from the version at  doi:10.1006/jtbi.2001.2293
 *
 *  Note that this is not a model for murine or human colonic-cell cycling, but is
 *  included in chaste as one of the most commonly known ODE based cell-cycle models.
 *
 *  Time taken to progress through the cycle is deterministic and given by
 *  an ODE system independent of external factors.
 */
class TysonNovakCellCycleModel : public AbstractOdeBasedCellCycleModel
{
private:

    friend class TestOdeBasedCellCycleModels;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
    }

public:

    /**
     * Default constructor.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    TysonNovakCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Initialise the cell-cycle model at the start of a simulation.
     *
     * This method will be called precisely once per cell set up in the initial
     * cell population. It is not called on cell division; use ResetForDivision(),
     * CreateCellCycleModel() and InitialiseDaughterCell() for that.
     *
     * By the time this is called, a CellPopulation will have been set up, so the model
     * can know where its cell is located in space. If relevant to the simulation,
     * the CellwiseData and/or other singletons will also have been initialised.
     */
    void Initialise();

    /**
     * Reset cell-cycle model by calling AbstractOdeBasedCellCycleModelWithStoppingEvent::ResetForDivision()
     * and setting initial conditions for protein concentrations.
     */
    void ResetForDivision();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Get the duration of the cell's S phase.
     */
    double GetSDuration();

    /**
     * Get the duration of the cell's G2 phase.
     */
    double GetG2Duration();

    /**
     * Get the duration of the cell's M phase.
     */
    double GetMDuration();

    /**
     * If the daughter cell type is stem, change it to transit.
     */
    void InitialiseDaughterCell();

    /**
     * Overridden GetAverageTransitCellCycleTime() method.
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     */
    double GetAverageStemCellCycleTime();

    /**
     * Overridden CanCellTerminallyDifferentiate() method.
     */
    bool CanCellTerminallyDifferentiate();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(TysonNovakCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(TysonNovakCellCycleModel)

#endif /*TYSONNOVAKCELLCYCLEMODEL_HPP_*/
