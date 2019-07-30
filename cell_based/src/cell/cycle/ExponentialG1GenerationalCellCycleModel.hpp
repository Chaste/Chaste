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

#ifndef EXPONENTIALG1GENERATIONALCELLCYCLEMODEL_HPP_
#define EXPONENTIALG1GENERATIONALCELLCYCLEMODEL_HPP_

#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A cell cycle model where the G1 duration is drawn from an exponential
 * random distribution. The rate parameter of this distribution is defined
 * by the member variable mRate.
 *
 * The default value of mRate is set to the inverse of the default value
 * of mTransitCellG1Duration in the constructor.
 *
 * mRate may be set using the method SetRate(), which also sets
 * mTransitCellG1Duration and mStemCellG1Duration to take the inverse of
 * mRate's new value. Similarly, the overridden methods SetTransitCellG1Duration()
 * and SetStemCellG1Duration() also set mRate to take the inverse of the
 * new value of mTransitCellG1Duration and mStemCellG1Duration, respectively.
 */
class ExponentialG1GenerationalCellCycleModel : public AbstractSimpleGenerationalCellCycleModel
{
    friend class TestSimpleCellCycleModels;

private:

    /**
     * The rate parameter of the exponential distribution, often denoted by lambda.
     * This parameter is initialised to 1.0/mTransitCellG1Duration. The latter variable is inherited from the abstract cell cycle model.
     * This initialisation ensures that cells will have the standard G1 duration (as defined in the abstract cell cycle model classes) on average.
     */
    double mRate;

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
        archive & boost::serialization::base_object<AbstractSimpleGenerationalCellCycleModel>(*this);

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
        archive & mRate;
    }

protected:

    /**
     * Stochastically set the G1 duration following an exponential distribution. Called on cell creation at
     * the start of a simulation, and for both parent and daughter
     * cells at cell division.
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
    ExponentialG1GenerationalCellCycleModel(const ExponentialG1GenerationalCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    ExponentialG1GenerationalCellCycleModel();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Return the rate parameter (lambda) of the exponential distribution.
     *
     * @return mRate, the rate parameter of the distribution
     */
    virtual double GetRate();

    /**
     * Set the rate parameter of the exponential distribution. For consistency,
     * this function also resets the internal values for mTransitCellG1Duration
     * and mStemCellG1Duration to the average value of the exponential random variable,
     * i.e. 1.0/rate. This ensures that GetAverageTransitCellCycleTime and GetAverageStemCellCycleTime
     * returns correct values.
     *
     * @param rate
     */
    virtual void SetRate(double rate);

    /**
     * Overridden SetStemCellG1Duration() method.
     *
     * Set mStemCellG1Duration to be a given value and also
     * set mRate to be the inverse of this value
     *
     * @param stemCellG1Duration  the new value of mStemCellG1Duration
     */
    virtual void SetStemCellG1Duration(double stemCellG1Duration);

    /**
     * Overridden SetTransitCellG1Duration() method.
     *
     * Set mTransitCellG1Duration to be a given value and also
     * set mRate to be the inverse of this value
     *
     * @param transitCellG1Duration  the new value of mTransitCellG1Duration
     */
    virtual void SetTransitCellG1Duration(double transitCellG1Duration);

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(ExponentialG1GenerationalCellCycleModel)

#endif /*EXPONENTIALG1GENERATIONALCELLCYCLEMODEL_HPP_*/
