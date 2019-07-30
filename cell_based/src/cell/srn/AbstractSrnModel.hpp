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

#ifndef ABSTRACTSRNMODEL_HPP_
#define ABSTRACTSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "Identifiable.hpp"

#include <boost/serialization/base_object.hpp>

#include <vector>

#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include "Cell.hpp"

class Cell; // Circular definition (cells need to know about SRN models and vice-versa)
typedef boost::shared_ptr<Cell> CellPtr;

/**
 * The AbstractSrnModel contains basic information to all sub-cellular reaction network (SRN) models.
 *
 * SRN models are noncopyable since cells are noncopyable.
 *
 * \todo #2752 Thoroughly document this class
 */
class AbstractSrnModel : public Identifiable
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the SRN model and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Make sure the SimulationTime singleton gets saved too
        SerializableSingleton<SimulationTime>* p_time_wrapper = SimulationTime::Instance()->GetSerializationWrapper();
        archive & p_time_wrapper;

        // DO NOT archive & mpCell; -- The SrnModel is only ever archived from the Cell
        // which knows this and it is handled in the load_construct of Cell.
        archive & mSimulatedToTime;
    }

    /**
     * Prevent copy-assignment of this class, or its subclasses.
     * Note that we do not define this method, therefore statements like "AbstractSrnModel new = old;" will not compile.
     * We do not inherit from boost::noncopyable because we *do* define a protected copy-constructor, for use by CreateSrnModel.
     *
     * @return the new SRN model.
     */
    AbstractSrnModel& operator=(const AbstractSrnModel&);

protected:

    /**
     * The cell that this model is associated with.
     */
    CellPtr mpCell;

    /**
     * The time the SRN model has been simulated to.
     */
    double mSimulatedToTime;

    /**
     * Protected copy-constructor for use by CreateSrnModel().  The only way for external code to create a copy of a SRN model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent SRN model will have had ResetForDivision() called just before CreateSrnModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the SRN model to copy.
     */
    AbstractSrnModel(const AbstractSrnModel& rModel);

public:

    /**
     * Sets up a new AbstractSrnModel, gives it a birth time of the
     * current simulation time (which is overwritten by some subclasses)
     */
    AbstractSrnModel();

    /**
     * Base class with virtual methods needs a virtual destructor. The destructor
     * does not delete mpCell. Instead, the cell takes responsibility for deleting
     * the SRN model when it is destroyed.
     */
    virtual ~AbstractSrnModel();

    /**
     * Gives the SRN model a pointer to its host cell.
     *
     * @param pCell pointer to the cell
     */
    void SetCell(CellPtr pCell);

    /**
     * Initialise the SRN model at the start of a simulation.
     *
     * This method will be called precisely once per cell set up in the initial
     * cell population. It is not called on cell division; use ResetForDivision(),
     * CreateSrnModel() and InitialiseDaughterCell() for that.
     *
     * By the time this is called, a CellPopulation will have been set up, so the model
     * can know where its cell is located in space. If relevant to the simulation,
     * any singletons will also have been initialised.
     */
    virtual void Initialise();

    /**
     * Initialise the new daughter cell's SRN model after a cell division.
     *
     * This is called by Cell::Divide() once the new cell object
     * has been fully created, to perform any initialisation of the
     * SRN which requires access to the cell.
     *
     * Note that much initialisation can be performed using the
     * combination of ResetForDivision() (called on the parent prior to
     * division) and CreateSrnModel() (called on the reset
     * parent to create the new cell-cycle model object).
     */
    virtual void InitialiseDaughterCell();

    /**
     * @return The cell which plays host to this SRN model.
     */
    CellPtr GetCell();

    /**
     * Set the time that the SRN simulation has run to.
     *
     * @param simulatedToTime the current simulation time of the SRN.
     */
    void SetSimulatedToTime(double simulatedToTime);

    /**
     * @return the time the SRN simulated has run to.
     */
    double GetSimulatedToTime() const;

    /**
     * Simulate the SRN to the current time.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     */
    virtual void SimulateToCurrentTime()=0;

    /**
     * Each SRN model must be able to be reset 'after' a cell division.
     *
     * Actually, this method is called from Cell::Divide() to
     * reset the SRN just before the daughter cell is created.
     * CreateSrnModel() can then clone our state to generate a
     * cell-cycle model instance for the daughter cell.
     */
    virtual void ResetForDivision();

    /**
     * Builder method to create new instances of the SRN model.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * This method is called by Cell::Divide() to create a SRN
     * model for the daughter cell.  Note that the parent SRN
     * model will have had ResetForDivision() called just before
     * CreateSrnModel() is called, so performing an exact copy of the
     * parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * Copy constructors are used to set all the member variables in the appropriate classes.
     *
     *  @return new srn model
     */
    virtual AbstractSrnModel* CreateSrnModel()=0;

    /**
     * Output SRN model used in the simulation to file and then call
     * OutputSrnParameters() to output all relevant parameters.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSrnModelInfo(out_stream& rParamsFile);

    /**
     * Output SRN model parameters to file.
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSrnModelParameters(out_stream& rParamsFile);
};

CLASS_IS_ABSTRACT(AbstractSrnModel)

#endif /* ABSTRACTSRNMODEL_HPP_ */
