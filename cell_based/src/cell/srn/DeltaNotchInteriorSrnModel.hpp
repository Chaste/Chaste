/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef DELTANOTCHINTERIORSRNMODEL_HPP_
#define DELTANOTCHINTERIORSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "DeltaNotchInteriorOdeSystem.hpp"
#include "AbstractOdeSrnModel.hpp"

/**
 * A subclass of AbstractOdeSrnModel that includes a Delta-Notch ODE system in 
 * the sub-cellular reaction network. This class represents model describing 
 * cytoplasmic (interior, nuclear) concentration of Delta Notch. Together with 
 * DeltaNotchEdgeSrn, this class can model feedback between cytoplasmic and 
 * junctional concentration, and cell-cell communication via edge (junctional) 
 * concentrations.
 */
class DeltaNotchInteriorSrnModel : public AbstractOdeSrnModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the SRN model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSrnModel>(*this);
    }

protected:
    /**
     * Protected copy-constructor for use by CreateSrnModel(). The only way for 
     * external code to create a copy of a SRN model is by calling that method, 
     * to ensure that a model of the correct subclass is created. 
     * 
     * This copy-constructor helps subclasses to ensure that all member 
     * variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a 
     * daughter cell upon cell division. Note that the parent SRN model will 
     * have had ResetForDivision() called just before CreateSrnModel() is 
     * called, so performing an exact copy of the parent is suitable behaviour. 
     * Any daughter-cell-specific initialisation can be done in 
     * InitialiseDaughterCell().
     *
     * @param rModel  the SRN model to copy.
     */
    DeltaNotchInteriorSrnModel(const DeltaNotchInteriorSrnModel& rModel);

public:

    /**
     * Default constructor calls base class.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    DeltaNotchInteriorSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of
     * this SRN model.
     *
     * @return a copy of the current SRN model.
     */
    virtual AbstractSrnModel* CreateSrnModel();

    /**
     * Overriden ResetForDivision() method for custom cell division mitosis
     */
    virtual void ResetForDivision();

    /**
     * Initialise the SRN model at the start of a simulation.
     *
     * This overridden method sets up a new Delta-Notch ODE system for cell interior.
     */
    virtual void Initialise();

    /**
     * Overridden SimulateToTime() method for custom behaviour.
     *
     */
    virtual void SimulateToCurrentTime();

    /**
     * Updates model parameters, such as total edge concnetration of Delta/Notch, via processing data
     * from CellData()() object
     */
    void UpdateDeltaNotch();

    /**
     * @return the current Notch level in this cell.
     */
    double GetNotch();

    /**
     * Set the notch level in this edge
     * @param value
     */
    void SetNotch(double value);

    /**
     * @return the current Delta level in this cell.
     */
    double GetDelta();

    /**
     * Set the delta level in this edge
     * @param value
     */
    void SetDelta(double value);

    /**
     * @return the current total level of Delta in the edges.
     */
    double GetTotalEdgeDelta();
    /**
     * @return total edge notch
     */
    double GetTotalEdgeNotch();

    /**
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSrnModelParameters(out_stream& rParamsFile);

    /**
     * Sets how much of Delta/Notch is returned back to interior after a junction is shrunk
     * Override the method to reflect user-define assumptions
     * 
     * @param pShrunkEdgeSrn Pointer to an SRN model
     */
    virtual void AddShrunkEdgeToInterior(AbstractSrnModel* pShrunkEdgeSrn);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchInteriorSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(DeltaNotchInteriorSrnModel)


#endif /* DELTANOTCHINTERIORSRNMODEL_HPP_ */
