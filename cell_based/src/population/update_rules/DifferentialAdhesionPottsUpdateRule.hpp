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
#ifndef DIFFERENTIALADHESIONUPDATERULE_HPP_
#define DIFFERENTIALADHESIONUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AdhesionPottsUpdateRule.hpp"
#include "PottsBasedCellPopulation.hpp"

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
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
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
    virtual ~DifferentialAdhesionPottsUpdateRule();

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
