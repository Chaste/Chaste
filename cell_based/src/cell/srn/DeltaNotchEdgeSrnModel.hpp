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

#ifndef DELTANOTCHEDGESRNMODEL_HPP_
#define DELTANOTCHEDGESRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "DeltaNotchEdgeOdeSystem.hpp"
#include "AbstractOdeSrnModel.hpp"

/**
 * A subclass of AbstractOdeSrnModel that includes a Delta-Notch ODE system in 
 * the sub-cellular reaction network. This SRN model represents a membrane/
 * cortex of a single junction of a cell. This class of models can be used 
 * together with DeltaNotchInteriorSrn models. The ODE model used here is an 
 * attempt to use previous work (see DeltaNotchSrnModel class) for more detailed 
 * description of Delta-Notch interactions involving edge quantities (this or 
 * neighbour edge information) and potentially coupling with cytoplasmic 
 * concentrations (DeltaNotchInteriorSrn class).
 */
class DeltaNotchEdgeSrnModel : public AbstractOdeSrnModel
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
    DeltaNotchEdgeSrnModel(const DeltaNotchEdgeSrnModel& rModel);

public:

    /**
     * Default constructor calls base class.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    DeltaNotchEdgeSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this SRN model.
     *
     * @return a copy of the current SRN model.
     */
    virtual AbstractSrnModel* CreateSrnModel();

    /**
     * Initialise the SRN model at the start of a simulation.
     *
     * This overridden method sets up a new Delta-Notch ODE system.
     */
    virtual void Initialise();

    /**
     * This method is called when a new edge is created (e.g. after cell division or T1 swap)
     */
    virtual void InitialiseDaughterCell();

    /**
     * Overridden SimulateToTime() method for custom behaviour.
     * Updates parameters (such as neighbour or interior Delta/Notch) and
     * runs the simulation to current time
     */
    virtual void SimulateToCurrentTime();

    /**
     * Update the levels of Delta and Notch of neighbouring edge sensed by this edge
     * That is, fetch neighbour values from CellEdgeData object, storing the sensed information,
     * into this model
     */
    void UpdateDeltaNotch();

    /**
     * @return the current Notch level in this edge.
     */
    double GetNotch();

    /**
     * Set the Notch level in this edge.
     * 
     * @param value the new Notch level in this edge
     */
    void SetNotch(double value);

    /**
     * @return the current Delta level in this edge.
     */
    double GetDelta();

    /**
     * Set the Delta level in this edge.
     * 
     * @param value the new Delta level in this edge
     */
    void SetDelta(double value);

    /**
     * @return the current level of Delta in the neighbouring cell's edge.
     */
    double GetNeighbouringDelta() const;

    /**
     * The value of interior Delta is stored as parameters in this model, which is
     * retrieved by this method
     * 
     * @return the level of Delta in cell interior
     */
    double GetInteriorDelta() const;

    /**
     * The value of interior Notch is stored as parameters in this model, which is
     * retrieved by this method
     * 
     * @return the level of Notch in cell interior
     */
    double GetInteriorNotch() const;

    /**
     * Output SRN model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSrnModelParameters(out_stream& rParamsFile);

    /**
     * Adds Delta/Notch from the input SRN model to this model.
     * Override the method declared in AbstractSrnModel class
     * 
     * @param pOtherSrn Pointer to an SRN
     * @param scale scale factor
     */
    virtual void AddSrnQuantities(AbstractSrnModel* pOtherSrn,
                                  const double scale = 1.0);

    /**
     * Here we assume that when a neighbouring junctions shrinks, 25% of its Delta/Notch
     * concentration is added to this edge
     * Override the method declared in AbstractSrnModel class
     * 
     * @param pShrunkEdgeSrn Pointer to an SRN
     */
    virtual void AddShrunkEdgeSrn(AbstractSrnModel* pShrunkEdgeSrn);

    /**
     * Here we add Delta/Notch when junctions merge via common vertex deletion
     * Override the method declared in AbstractSrnModel class
     * 
     * @param pMergedEdgeSrn Pointer to an SRN
     */
    virtual void AddMergedEdgeSrn(AbstractSrnModel* pMergedEdgeSrn);

    /**
     * By default, Edge concentrations are split according to the lengths of two new edges , when an edge is split, 
     * relative to the original edge they are created from.
     * Override the method declared in AbstractSrnModel class
     * 
     * @param relativePosition position of the node splitting the two new edges relative to the original edge, used when splitting Edge concentrations
     */
    virtual void SplitEdgeSrn(const double relativePosition);
};

typedef boost::shared_ptr<DeltaNotchEdgeSrnModel> DeltaNotchEdgeSrnModelPtr;

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchEdgeSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(DeltaNotchEdgeSrnModel)

#endif  /* DELTANOTCHEDGESRNMODEL_HPP_ */
