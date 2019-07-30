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

#ifndef ABSTRACTODESRNMODEL_HPP_
#define ABSTRACTODESRNMODEL_HPP_

#include <vector>
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "AbstractSrnModel.hpp"
#include "CellCycleModelOdeHandler.hpp"
#include "SimulationTime.hpp"


/**
 * This class contains the abstract code for an ODE sub-cellular reaction network (SRN) model
 *  - based on AbstractOdeBasedCellCycleModel
 *
 * \todo #2752 Thoroughly document this class
 *
 */
class AbstractOdeSrnModel : public AbstractSrnModel, public CellCycleModelOdeHandler
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
        archive & boost::serialization::base_object<AbstractSrnModel>(*this);
        archive & boost::serialization::base_object<CellCycleModelOdeHandler>(*this);
        archive & mInitialConditions;
        archive & mStateSize;
    }

protected:

    /**
     * The initial condition for the ODE state variables.
     */
    std::vector<double> mInitialConditions;

    /**
     * The number of state variables.
     */
    unsigned mStateSize;

    using AbstractSrnModel::Initialise;
    /**
     * Overridden Initialise() method, which here sets up the ODE system.
     *
     * Note we bring virtual functions from AbstractSrnModel into derived namespace so overloading virtual works.
     *
     * @param pOdeSystem pointer to an ODE system
     */
    void Initialise(AbstractOdeSystem* pOdeSystem);

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
     * @param rModel  the SRN model to copy.
     */
    AbstractOdeSrnModel(const AbstractOdeSrnModel& rModel);

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
     * call AbstractSrnModel::ResetForDivision() from inside their version.
     */
    virtual void ResetForDivision();

    /**
     * Set mInitialConditions. Used in CreateSrnModel().
     *
     * @param initialConditions the new value of mInitialConditions
     */
    void SetInitialConditions(std::vector<double> initialConditions);

    /**
     * Outputs SRN model parameters to file. Virtual void so needs to be specified in child classes.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSrnModelParameters(out_stream& rParamsFile) = 0;
};

CLASS_IS_ABSTRACT(AbstractOdeSrnModel)

#endif /* ABSTRACTODESRNMODEL_HPP_ */
