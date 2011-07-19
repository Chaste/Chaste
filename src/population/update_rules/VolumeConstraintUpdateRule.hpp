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
#ifndef VOLUMECONSTRAINTUPDATERULE_HPP_
#define VOLUMECONSTRAINTUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractPottsUpdateRule.hpp"
#include "PottsBasedCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"

/**
 * A volume constraint update rule class for use in Potts based simulations.
 *
 * Note this currently assumes cells don't grow, i.e the target Volume is constant
 * for each cell over time.
 *
 */
template<unsigned DIM>
class VolumeConstraintUpdateRule  : public AbstractPottsUpdateRule<DIM>
{
friend class TestPottsUpdateRules;

private:

	/**
	 * Cell deformation energy parameter. Has units of ?
	 */
	double mDeformationEnergyParameter;

    /**
     * Non-dimensional target Volume of a mature (fully-grown) Cell, given in number of lattice sites
     */
    double mMatureCellTargetVolume;


    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'

    	///\todo implement archiving (see #1685)
        //archive & boost::serialization::base_object<AbstractPottsUpdateRule<DIM> >(*this);
        //archive & mDeformationEnergyParameter;
    	//archive & mMatureCellTargetVolume;
    }

public:

    /**
     * Constructor.
     */
    VolumeConstraintUpdateRule();

    /**
     * Destructor.
     */
    ~VolumeConstraintUpdateRule();

    /**
	 * Overridden EvaluateHamiltonianContribution() method
	 *
	 * Uses  sum_elements alpha (V_i - V_i^T)^2
	 *
	 * @param currentNodeIndex The index of the current node/lattice site
	 * @param targetNodeIndex The index of the target node/lattice site
	 * @param rCellPopulation The cell population
	 *
	 * @return The difference in the Hamiltonian with the current configuration and
	 * the configuration with the target node having the same spin as the current node.
	 */
    double EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                           unsigned targetNodeIndex,
                                           AbstractCellPopulation<2>& rCellPopulation);

    /**
     * @return mDeformationEnergyParameter
     */
    double GetDeformationEnergyParameter();

    /**
     * Set mDeformationEnergyParameter.
     *
     * @param deformationEnergyParameter the new value of mDeformationEnergyParameter
     */
    void SetDeformationEnergyParameter(double deformationEnergyParameter);


    /**
     * @return mMatureCellTargetVolume
     */
    double GetMatureCellTargetVolume() const;

    /**
     * Set mMatureCellTargetVolume.
     *
     * @param matureCellTargetVolume the new value of mMatureCellTargetVolume
     */
    void SetMatureCellTargetVolume(double matureCellTargetVolume);

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VolumeConstraintUpdateRule)

#endif /*VOLUMECONSTRAINTUPDATERULE_HPP_*/

