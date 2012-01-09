/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
