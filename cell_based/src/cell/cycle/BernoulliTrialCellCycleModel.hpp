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

#ifndef BERNOULLITRIALCELLCYCLEMODEL_HPP_
#define BERNOULLITRIALCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"

/**
 * Simple cell-cycle model where mature non-differentiated cells have a specified probability of
 * dividing per hour.
 *
 * The class includes two parameters: the first, mDivisionProbability, defines the probability
 * of dividing per hour; the second, mMinimumDivisionAge, defines a minimum age at which cells
 * may divide.
 */
class BernoulliTrialCellCycleModel : public AbstractCellCycleModel
{
private:

    friend class boost::serialization::access;

    /**
     * Boost Serialization method for archiving/checkpointing
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mDivisionProbability;
        archive & mMinimumDivisionAge;
    }

protected:
    /**
     * Probability of dividing per hour.
     * Defaults to 0.1.
     */
    double mDivisionProbability;

    /**
     * Minimum age of a cell at which it may divide.
     * Defaults to 1 hour.
     */
    double mMinimumDivisionAge;

    /**
     * Protected copy-constructor for use by CreateCellCycleModel.
     *
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass
     * is created. This copy-constructor helps subclasses to ensure that all
     * member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a
     * daughter cell upon cell division. Note that the parent cell cycle model
     * will have had ResetForDivision() called just before CreateCellCycleModel()
     * is called, so performing an exact copy of the parent is suitable behaviour.
     * Any daughter-cell-specific initialisation can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    BernoulliTrialCellCycleModel(const BernoulliTrialCellCycleModel& rModel);

public:

    /**
     * Constructor.
     */
    BernoulliTrialCellCycleModel();

    /**
     * Overridden ReadyToDivide() method.
     *
     * If the cell's age is greater than mMinimumDivisionAge, then we draw a uniform
     * random number r ~ U[0,1]. If r < mDivisionProbability*dt, where dt is the
     * simulation time step, then the cell is ready to divide and we return true.
     * Otherwise, the cell is not yet ready to divide and we return false.
     *
     * @return whether the cell is ready to divide.
     */
    virtual bool ReadyToDivide();

    /**
     * Overridden builder method to create new instances of
     * the cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Set the value of mDivisionProbability.
     *
     * @param divisionProbability the new value of mDivisionProbability
     */
    void SetDivisionProbability(double divisionProbability);

    /**
     * Get mDivisionProbability.
     *
     * @return mDivisionProbability
     */
    double GetDivisionProbability();

    /**
     * Set the value of mMinimumDivisionAge.
     *
     * @param minimumDivisionAge the new value of mMinimumDivisionAge
     */
    void SetMinimumDivisionAge(double minimumDivisionAge);

    /**
     * Get mMinimumDivisionAge.
     *
     * @return mMinimumDivisionAge
     */
    double GetMinimumDivisionAge();

    /**
     * Overridden GetAverageTransitCellCycleTime() method.
     *
     * @return the average cell cycle time for cells of transit proliferative type
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     *
     * @return the average cell cycle time for cells of stem proliferative type
     */
    double GetAverageStemCellCycleTime();

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BernoulliTrialCellCycleModel)

#endif // BERNOULLITRIALCELLCYCLEMODEL_HPP_
