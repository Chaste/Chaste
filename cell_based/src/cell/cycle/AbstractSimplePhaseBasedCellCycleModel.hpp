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

#ifndef ABSTRACTSIMPLEPHASEBASEDCELLCYCLEMODEL_HPP_
#define ABSTRACTSIMPLEPHASEBASEDCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractPhaseBasedCellCycleModel.hpp"

/**
 * This class contains all the functionality shared by 'simple' cell-cycle models,
 * where the duration of each cell cycle phase is determined when the cell-cycle
 * model is created. Note that whether or not the cell should actually divide may
 * still depend on further conditions in subclasses; for example, the cell may only
 * divide if the local concentration of a signalling molecule is sufficiently high/
 *
 * This class of cell-cycle models is distinct from 'ODE-based' cell-cycle models,
 * where the duration of one or more cell cycle phases are evaluated 'on the fly'
 * as the cell ages, according to a system of ordinary differential equations (ODEs)
 * governing (for example) the concentrations of key intracellular proteins.
 */
class AbstractSimplePhaseBasedCellCycleModel : public AbstractPhaseBasedCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractPhaseBasedCellCycleModel>(*this);
    }

protected:

    /**
     * Subclasses can override this function if they wish, this just
     * allocates the default values for each of the different cell
     * types' G1 durations as defined in AbstractPhaseBasedCellCycleModel.
     */
    virtual void SetG1Duration();

    /**
     * Protected copy-constructor for use by CreateCellCycleModel.
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    AbstractSimplePhaseBasedCellCycleModel(const AbstractSimplePhaseBasedCellCycleModel& rModel);

public:

    /**
     * Default constructor - creates an AbstractSimplePhaseBasedCellCycleModel.
     */
    AbstractSimplePhaseBasedCellCycleModel();

    /**
     * Destructor.
     */
    virtual ~AbstractSimplePhaseBasedCellCycleModel();

    /** See AbstractPhaseBasedCellCycleModel::ResetForDivision() */
    virtual void ResetForDivision();

    /**
     * Default UpdateCellCyclePhase() method for a simple cell-cycle model.
     *
     * Can be overridden if they should do something more subtle.
     */
    virtual void UpdateCellCyclePhase();

    /**
     * Set the new cell's G1 duration once it has been created after division.
     * The duration will be based on cell type.
     */
    void InitialiseDaughterCell();

    /** See AbstractPhaseBasedCellCycleModel::Initialise() */
    virtual void Initialise();

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

CLASS_IS_ABSTRACT(AbstractSimplePhaseBasedCellCycleModel)

#endif /*ABSTRACTSIMPLEPHASEBASEDCELLCYCLEMODEL_HPP_*/
