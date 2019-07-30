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

#ifndef SIMPLEWNTCELLCYCLEMODEL_HPP_
#define SIMPLEWNTCELLCYCLEMODEL_HPP_

#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "WntConcentration.hpp"

/**
 * Simple Wnt-dependent cell-cycle model.
 */
class SimpleWntCellCycleModel : public AbstractSimplePhaseBasedCellCycleModel
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
        archive & boost::serialization::base_object<AbstractSimplePhaseBasedCellCycleModel>(*this);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;

        archive & mUseCellProliferativeTypeDependentG1Duration;
        archive & mWntStemThreshold;
        archive & mWntTransitThreshold;
        archive & mWntLabelledThreshold;
    }

protected:

    /**
     * Whether to use different mean G1 durations for different cell types.
     * For use in SetG1Duration().
     */
    bool mUseCellProliferativeTypeDependentG1Duration;

    /**
     * Non-dimensionalized Wnt threshold, above which cells behave as stem cells.
     */
    double mWntStemThreshold;

    /**
     * Non-dimensionalized Wnt threshold, above which cells progress through the cell cycle.
     */
    double mWntTransitThreshold;

    /**
     * Non-dimensionalized Wnt threshold, above which labelled cells progress through the cell cycle.
     */
    double mWntLabelledThreshold;

    /**
     * @return the Wnt level experienced by the cell.
     */
    double GetWntLevel() const;

    /**
     * @return the type of Wnt concentration (LINEAR, RADIAL, EXPONENTIAL or NONE).
     * This affects how the cell cycle phase is updated.
     */
    WntConcentrationType GetWntType();

    /**
     * Stochastically set the G1 duration. The G1 duration is taken
     * from a normal distribution whose mean is the G1 duration given
     * in AbstractPhaseBasedCellCycleModel, for the cell type, and whose standard deviation
     * is 1.
     *
     * Called on cell creation at the start of a simulation, and for both
     * parent and daughter cells at cell division.
     */
    void SetG1Duration();

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
    SimpleWntCellCycleModel(const SimpleWntCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractPhaseBasedCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called.
     */
    SimpleWntCellCycleModel();

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    virtual void UpdateCellCyclePhase();

    /**
     * Overridden InitialiseDaughterCell() method.
     */
    virtual void InitialiseDaughterCell();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return the new cell-cycle model
     */
    virtual AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * @return mUseCellProliferativeTypeDependentG1Duration
     */
    bool GetUseCellProliferativeTypeDependentG1Duration() const;


    /**
     * Set whether Whether the duration of the G1 phase is dependent on cell type
     * @param useCellProliferativeTypeDependentG1Duration - boolean, defaults to true.
     */
    void SetUseCellProliferativeTypeDependentG1Duration(bool useCellProliferativeTypeDependentG1Duration=true);

    /**
     * Overridden CanCellTerminallyDifferentiate() method.
     * @return whether cell can terminally differentiate
     */
    virtual bool CanCellTerminallyDifferentiate();

    /**
     * @return mWntStemThreshold
     */
    double GetWntStemThreshold() const;

    /**
     * Set mWntStemThreshold.
     *
     * @param wntStemThreshold the value of mWntStemThreshold
     */
    void SetWntStemThreshold(double wntStemThreshold);

    /**
     * @return mWntTransitThreshold
     */
    double GetWntTransitThreshold() const;

    /**
     * Set mWntTransitThreshold.
     *
     * @param wntTransitThreshold the value of mWntTransitThreshold
     */
    void SetWntTransitThreshold(double wntTransitThreshold);

    /**
     * @return mWntLabelledThreshold
     */
    double GetWntLabelledThreshold() const;

    /**
     * Set mWntLabelledThreshold.
     *
     * @param wntLabelledThreshold the value of mWntLabelledThreshold
     */
    void SetWntLabelledThreshold(double wntLabelledThreshold);

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SimpleWntCellCycleModel)

#endif /*SIMPLEWNTCELLCYCLEMODEL_HPP_*/
