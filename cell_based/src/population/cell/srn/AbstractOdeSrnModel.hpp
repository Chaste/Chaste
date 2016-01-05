/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef ABSTRACTODESRN_HPP_
#define ABSTRACTODESRN_HPP_

#include <vector>
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractSrnModel.hpp"
#include "CellCycleModelOdeHandler.hpp"
#include "SimulationTime.hpp"

//class AbstractSrnModel; ///\todo #2752 remove this commented code

/**
 * This class contains the abstract code for an ODE SRN model
 *  - based on AbstractOdeBasedCellCycleModel
 *
 * \todo #2752 Thoroughly document this class
 */
class AbstractOdeSrnModel : public AbstractSrnModel, public CellCycleModelOdeHandler
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the srn model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSrnModel>(*this);
        archive & boost::serialization::base_object<CellCycleModelOdeHandler>(*this);
        archive & mFinishedRunningOdes;
        archive & mInitialConditions;
        archive & mStateSize;
    }

protected:

    /**
     * Whether the SRN model is currently in a delay (not solving ODEs).
     */
    bool mFinishedRunningOdes;

    /**
     * The initial condition for the ODE state variables.
     */
    std::vector<double> mInitialConditions;

    /**
     * The number of state variables.
     */
    unsigned mStateSize;

    /*
     * Overridden Initialise method which here ets up the ODE system.
     *
     *Note we bring virtual funcs from AbstractSrnModel into derived namespace so overloading virtual works.
     */
    using AbstractSrnModel::Initialise;
    void Initialise(AbstractOdeSystem* pOdeSystem);

    /*
     * Overridden CreateSrnModel() method.
     *
     * Builder method to create new instances of the SRN model.
     * Each concrete subclass must implement this method to create an
     * instance of that subclass.
     *
     * This method is called by Cell::Divide() to create a SRN
     * model for the daughter cell.  Note that the parent SRN
     * model will have had ResetForDivision() called just before
     * CreateSrnModel() is called, so performing an exact copy of the
     * parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @return new SRN model
     *
     * Note bring virtual functions from AbstractSrnModel into derived namespace so overloading virtual works.
     */
    using AbstractSrnModel::CreateSrnModel;
    AbstractSrnModel* CreateSrnModel(AbstractOdeSrnModel* pModel);

public:
    /**
     * Create an AbstractOdeSrnModel.
     *
     * @param stateSize The number of state variables in the ODE system.
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers).
     */
    AbstractOdeSrnModel(unsigned stateSize,
                        boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Destructor.
     */
    virtual ~AbstractOdeSrnModel();

    /**
     * Here we solve the ODEs associated with the SRN.
     */
    virtual void SimulateToCurrentTime();

     /**
     * For a naturally cycling model this does not need to be overridden in the
     * subclasses. But most models should override this function and then
     * call AbstractOdeBasedCellCycleModel::ResetForDivision() from inside their version.
     */
    virtual void ResetForDivision();

    /**
     * Set mFinishedRunningOdes. Used in CreateSrnModel().
     *
     * @param finishedRunningOdes the new value of mFinishedRunningOdes
     */
    void SetFinishedRunningOdes(bool finishedRunningOdes);

    /**
     * Set mInitialConditions. Used in CreateSrnModel().
     *
     * @param initialConditions the new value of mInitialConditions
     */
    void SetInitialConditions(std::vector<double> initialConditions);

    /**
     * @return reference to the ODE system state values.
     */
    const std::vector<double>& GetStateVariables();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSrnModelParameters(out_stream& rParamsFile);
};

CLASS_IS_ABSTRACT(AbstractOdeSrnModel)

#endif /* ABSTRACTODESRN_HPP_ */
