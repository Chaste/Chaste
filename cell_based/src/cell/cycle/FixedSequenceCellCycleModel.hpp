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

#ifndef FIXEDSEQUENCECELLCYCLEMODEL_HPP_
#define FIXEDSEQUENCECELLCYCLEMODEL_HPP_

#include "ExponentialG1GenerationalCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellCycleTimesGenerator.hpp"
#include "AbstractCellCycleModel.hpp"

/**
 * This FixedSequenceCellCycleModel generates a sequence of randomly distributed
 * cell cycle times at the beginning of the simulation rather than drawing random numbers
 * during the simulation.  The random seed for the generation of these random numbers
 * may be specified by the user and the cell cycle times are stored in a separate
 * singleton, the CellCycleTimesGenerator.  This gives the user control over the exact cell cycle
 * times used in the simulation while not interfering with other random events implemented in
 * the simulation algorithm.
 *
 * This cell cycle model is an extended version of the ExponentialG1GenerationalCellCycleModel.
 * Similar to ExponentialG1GenerationalCellCycleModel, the G1 duration is exponentially random
 * distributed, and the interface to FixedSequenceCellCycleModel and
 * ExponentialG1GenerationalCellCycleModel are similar.  In contrast to ExponentialG1GenerationalCellCycleModel
 * this cell cycle model does not treat stem cells and transit cells differently.
 *
 * If using this cell cycle model, the following line must be called before
 * starting the simulation:
 *
 * CellCycleTimesGenerator::Instance()->GenerateCellCycleTimeSequence();
 *
 * Further, the following line must be run at the end of the simulation:
 *
 * CellCycleTimesGenerator::Destroy();
 */
class FixedSequenceCellCycleModel : public ExponentialG1GenerationalCellCycleModel
{

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<ExponentialG1GenerationalCellCycleModel>(*this);
        // Make sure the CellCycleTimesGenerator singleton gets saved too
        SerializableSingleton<CellCycleTimesGenerator>* p_wrapper = CellCycleTimesGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
    }

protected:

    /**
     * Stochastically set the G1 duration following the sequence in CellCycleTimesGenerator. Called on cell creation at
     * the start of a simulation, and for both parent and daughter
     * cells at cell division.
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
    FixedSequenceCellCycleModel(const FixedSequenceCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    FixedSequenceCellCycleModel();

    /**
     * Overridden method to create new cell cycle models after division
     *
     * @return new cell cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Overridden SetRate method.
     *
     * @param rate  the new value of the rate parameter
     */
    void SetRate(double rate);

    /**
     * Overridden GetRate method.
     *
     * @returns rate  the value of the rate parameter
     */
    double GetRate();

    /**
     * Overridden SetStemCellG1Duration() method.
     *
     * Set mStemCellG1Duration to be a given value and also
     * set mRate to be the inverse of this value
     *
     * @param stemCellG1Duration  the new value of mStemCellG1Duration
     */
    void SetStemCellG1Duration(double stemCellG1Duration);

    /**
     * Overridden SetTransitCellG1Duration() method.
     *
     * Set mTransitCellG1Duration to be a given value and also
     * set mRate to be the inverse of this value
     *
     * @param transitCellG1Duration  the new value of mTransitCellG1Duration
     */
    void SetTransitCellG1Duration(double transitCellG1Duration);

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(FixedSequenceCellCycleModel)

#endif /*FIXEDSEQUENCECELLCYCLEMODEL_HPP_*/
