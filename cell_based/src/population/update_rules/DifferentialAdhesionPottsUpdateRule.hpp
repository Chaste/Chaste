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
#ifndef DIFFERENTIALADHESIONUPDATERULE_HPP_
#define DIFFERENTIALADHESIONUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AdhesionPottsUpdateRule.hpp"
#include "PottsBasedCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

/**
 * An adhesion update rule for use in cell-based simulations
 * using the cellular Potts model. This rule implements differential adhesion
 * between unlabelled and labelled cells
 */
template<unsigned DIM>
class DifferentialAdhesionPottsUpdateRule : public AdhesionPottsUpdateRule<DIM>
{
friend class TestPottsUpdateRules;

private:

    /**
     * LabelledCell-LabelledCell adhesion energy parameter.
     * Set to the default value 0.1 in the constructor.
     * \todo provide units
     */
    double mLabelledCellLabelledCellAdhesionEnergyParameter;

    /**
     * LablledCell-cell adhesion energy parameter.
     * Set to the default value 0.1 in the constructor.
     * \todo provide units
     */
    double mLabelledCellCellAdhesionEnergyParameter;

    /**
     * LabelledCell-boundary adhesion energy parameter.
     * Set to the default value 0.2 in the constructor.
     * \todo provide units
     */
    double mLabelledCellBoundaryAdhesionEnergyParameter;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AdhesionPottsUpdateRule<DIM> >(*this);
        archive & mLabelledCellLabelledCellAdhesionEnergyParameter;
        archive & mLabelledCellCellAdhesionEnergyParameter;
        archive & mLabelledCellBoundaryAdhesionEnergyParameter;
    }

public:

    /**
     * Constructor.
     */
    DifferentialAdhesionPottsUpdateRule();

    /**
     * Destructor.
     */
    ~DifferentialAdhesionPottsUpdateRule();

    /**
     * Overridden GetCellCellAdhesionEnergy method to implement differential adhesion.
     *
     * @param pCellA pointer to the 1st cell
     * @param pCellB pointer to the 2nd cell
     *
     * @return The cell cell interaction adhesion energy between the two cells
     */
    virtual double GetCellCellAdhesionEnergy(CellPtr pCellA, CellPtr pCellB);

    /**
     * Overridden GetCellBoundaryAdhesionEnergy method to implement differential adhesion.
     *
     * @param pCell pointer to the cell
     *
     * @return Cell boundary interaction adhesion energy for the cell
     */
    virtual double GetCellBoundaryAdhesionEnergy(CellPtr pCell);

    /**
     * @return mLabelledCellLabelledCellAdhesionEnergyParameter
     */
    double GetLabelledCellLabelledCellAdhesionEnergyParameter();

    /**
     * @return mLabelledCellCellAdhesionEnergyParameter
     */
    double GetLabelledCellCellAdhesionEnergyParameter();

    /**
     * @return mLabelledCellBoundaryAdhesionEnergyParameter
     */
    double GetLabelledCellBoundaryAdhesionEnergyParameter();

    /**
     * Set mLabelledCellLabelledCellAdhesionEnergyParameter.
     *
     * @param labelledCellLabelledCellAdhesionEnergyParameter the new value of mLabelledCellLabelledCellAdhesionEnergyParameter
     */
    void SetLabelledCellLabelledCellAdhesionEnergyParameter(double labelledCellLabelledCellAdhesionEnergyParameter);

    /**
     * Set mLabelledCellCellAdhesionEnergyParameter.
     *
     * @param labelledCellCellAdhesionEnergyParameter the new value of mLabelledCelldCellAdhesionEnergyParameter
     */
    void SetLabelledCellCellAdhesionEnergyParameter(double labelledCellCellAdhesionEnergyParameter);

    /**
     * Set mLabelledCellBoundaryAdhesionEnergyParameter.
     *
     * @param labelledCellBoundaryAdhesionEnergyParameter the new value of mLabelledCellBoundaryAdhesionEnergyParameter
     */
    void SetLabelledCellBoundaryAdhesionEnergyParameter(double labelledCellBoundaryAdhesionEnergyParameter);

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DifferentialAdhesionPottsUpdateRule)

#endif /*DIFFERENTIALADHESIONUPDATERULE_HPP_*/
