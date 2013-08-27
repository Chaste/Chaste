/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef CONTACTINHIBITIONCELLCYCLEMODEL_HPP_
#define CONTACTINHIBITIONCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"

/**
 * Simple stress-based cell-cycle model.
 *
 * A simple stress-dependent cell-cycle model that inherits from
 * AbstractSimpleCellCycleModel. The duration of G1 phase depends
 * on the local stress, interpreted here as deviation from target
 * volume (or area/length in 2D/1D).
 *
 * This model allows for quiescence imposed by transient periods
 * of high stress, followed by relaxation.
 *
 * Note that in this cell cycle model, quiescence is implemented
 * by extending the G1 phase. If a cell is compressed during G2
 * or S phases then it will still divide, and thus cells whose
 * volumes are smaller than the given threshold may still divide.
 */
class ContactInhibitionCellCycleModel : public AbstractSimpleCellCycleModel
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
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
        archive & mQuiescentVolumeFraction;
        archive & mEquilibriumVolume;
        archive & mCurrentQuiescentDuration;
        archive & mCurrentQuiescentOnsetTime;
    }

    /**
     * The fraction of the cells' equilibrium volume in G1 phase below which these cells are quiescent.
     */
    double mQuiescentVolumeFraction;

    /**
     * The cell equilibrium volume while in G1 phase.
     */
    double mEquilibriumVolume;

    /**
     * The time when the current period of quiescence began.
     */
    double mCurrentQuiescentOnsetTime;

    /**
     * How long the current period of quiescence has lasted.
     * Has units of hours.
     */
    double mCurrentQuiescentDuration;

public:

    /**
     * Constructor.
     */
    ContactInhibitionCellCycleModel();

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    void UpdateCellCyclePhase();

    /**
     * Overridden builder method to create new instances of
     * the cell-cycle model.
     *
     * @return new cell-cycle model
     *
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * @param quiescentVolumeFraction
     */
    void SetQuiescentVolumeFraction(double quiescentVolumeFraction);

    /**
     * @return mQuiescentVolumeFraction
     */
    double GetQuiescentVolumeFraction();

    /**
     * @param equilibriumVolume
     */
    void SetEquilibriumVolume(double equilibriumVolume);

    /**
     * @return mEquilibriumVolume
     */
    double GetEquilibriumVolume();

    /**
    * Set method for mCurrentQuiescentDuration.
    *
    * @param currentQuiescentDuration the new value of mCurrentQuiescentDuration
    */
    void SetCurrentQuiescentDuration(double currentQuiescentDuration);

    /**
     * @return mCurrentQuiescentDuration
     */
    double GetCurrentQuiescentDuration();

    /**
    * Set method for mCurrentQuiescentOnsetTime.
    *
    * @param currentQuiescentOnsetTime the new value of mCurrentQuiescentOnsetTime
    */
    void SetCurrentQuiescentOnsetTime(double currentQuiescentOnsetTime);

    /**
     * @return mCurrentQuiescentOnsetTime
     */
    double GetCurrentQuiescentOnsetTime();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ContactInhibitionCellCycleModel)

#endif // CONTACTINHIBITIONCELLCYCLEMODEL_HPP_
