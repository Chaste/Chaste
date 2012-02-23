/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef ABSTRACTVANLEEUWEN2009WNTSWATCELLCYCLEMODEL_HPP_
#define ABSTRACTVANLEEUWEN2009WNTSWATCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>

#include <cfloat>

#include "AbstractOdeSystem.hpp"
#include "AbstractWntOdeBasedCellCycleModel.hpp"
#include "VanLeeuwen2009WntSwatCellCycleOdeSystem.hpp"
#include "AbstractCellMutationState.hpp"

/**
 * Wnt-dependent cell-cycle model.
 */
class AbstractVanLeeuwen2009WntSwatCellCycleModel : public AbstractWntOdeBasedCellCycleModel
{
private:

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
        archive & boost::serialization::base_object<AbstractWntOdeBasedCellCycleModel>(*this);
    }

    /**
     * Called by ::Initialise() and ::UpdateCellProliferativeType() only.
     * Updates mCellProliferativeType to match mpOdeSystem's beta-catenin levels
     *
     * This carries out the work for ::UpdateCellProliferativeType();
     * But does not check the current time so it can be used by the initialise method.
     */
    void ChangeCellProliferativeTypeDueToCurrentBetaCateninLevel();

    /**
     * Adjust any ODE parameters needed before solving until currentTime.
     *
     * @param currentTime  the time up to which the system will be solved.
     */
    void AdjustOdeParameters(double currentTime);

public:

    /**
     * Default constructor.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    AbstractVanLeeuwen2009WntSwatCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * See AbstractCellCycleModel::Initialise()
     *
     * In this case we set up a new ODE system for a daughter cell.
     */
    void Initialise();

    /**
     * @return the level of membrane bound beta-catenin. To be used in cell-cell adhesion calculations.
     */
    double GetMembraneBoundBetaCateninLevel();

    /**
     * @return the level of cytoplasmic beta-catenin (including ubiquitinated - awaiting degradation)
     */
    double GetCytoplasmicBetaCateninLevel();

    /**
     * @return the level of nuclear beta-catenin. To be used in transcription
     */
    double GetNuclearBetaCateninLevel();

    /**
     * Pure virtual method to be implemented in concrete classes, which
     * should should allocate the mOdeSystem variable using the appropriate
     * hypothesis (one or two).
     *
     * @param wntConcentration Wnt concentration
     * @param pMutationState Mutation state
     */
    virtual void InitialiseOdeSystem(double wntConcentration, boost::shared_ptr<AbstractCellMutationState> pMutationState)=0;

    /**
     * Outputs cell-cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

CLASS_IS_ABSTRACT(AbstractVanLeeuwen2009WntSwatCellCycleModel)

#endif /*ABSTRACTVANLEEUWEN2009WNTSWATCELLCYCLEMODEL_HPP_*/

