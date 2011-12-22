/*

Copyright (C) University of Oxford, 2005-2011

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

#ifndef FIXEDDURATIONGENERATIONBASEDCELLCYCLEMODEL_HPP_
#define FIXEDDURATIONGENERATIONBASEDCELLCYCLEMODEL_HPP_

#include "AbstractSimpleGenerationBasedCellCycleModel.hpp"

/**
 *  Fixed cell-cycle model.
 *
 *  Cell cycle time is deterministic for stem and transit cells (with values
 *  StemCellG1Duration + SG2MDuration
 *  and TransitCellG1Duration + SG2MDuration (values found in AbstractCellCycleModel))
 */
class FixedDurationGenerationBasedCellCycleModel : public AbstractSimpleGenerationBasedCellCycleModel
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
        archive & boost::serialization::base_object<AbstractSimpleGenerationBasedCellCycleModel>(*this);
    }

public:

    /**
     * Default constructor. Note that mBirthTime is set in
     * AbstractCellCycleModel() and mG1Duration is set in
     * AbstractSimpleCellCycleModel().
     */
    FixedDurationGenerationBasedCellCycleModel();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(FixedDurationGenerationBasedCellCycleModel)

#endif /*FIXEDDURATIONGENERATIONBASEDCELLCYCLEMODEL_HPP_*/
