/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "ImmersedBoundarySingleCellMigrationForce.hpp"

template <unsigned DIM>
SingleCellMigrationForce<DIM>::SingleCellMigrationForce()
		: AbstractImmersedBoundaryForce<DIM>(),
          mStrength(1.0 * 1e2),
          mElementIndex(0)
{
	/**
	 * This force class applies a force of strength mStrength and in the direction mDirection.
	 * The default the strength is 100 and default direction is x_unit.
	 */
	unit_vector<double> x_unit(DIM, 0);
	mDirection = x_unit;
}

template <unsigned DIM>
SingleCellMigrationForce<DIM>::~SingleCellMigrationForce()
{
}

template <unsigned DIM>
void SingleCellMigrationForce<DIM>::AddImmersedBoundaryForceContribution(std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
																		 ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)
{
    // Data common across the entire cell population
    double intrinsicSpacingSquared = rCellPopulation.GetIntrinsicSpacing() * rCellPopulation.GetIntrinsicSpacing();

    // Loop over all elements ( <DIM, DIM> )
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryElementIterator elem_it = rCellPopulation.rGetMesh().GetElementIteratorBegin();
         elem_it != rCellPopulation.rGetMesh().GetElementIteratorEnd();
         ++elem_it)
    {
        CalculateForcesOnElement(*elem_it, rCellPopulation, intrinsicSpacingSquared);
    }

    // Loop over all laminas ( <DIM-1, DIM> )
    for (typename ImmersedBoundaryMesh<DIM, DIM>::ImmersedBoundaryLaminaIterator lam_it = rCellPopulation.rGetMesh().GetLaminaIteratorBegin();
         lam_it != rCellPopulation.rGetMesh().GetLaminaIteratorEnd();
         ++lam_it)
    {
        CalculateForcesOnElement(*lam_it, rCellPopulation, intrinsicSpacingSquared);
    }
}

template <unsigned DIM>
template <unsigned ELEMENT_DIM>
void ImmersedBoundaryLinearMembraneForce<DIM>::CalculateForcesOnElement(ImmersedBoundaryElement<ELEMENT_DIM, DIM>& rElement,
                                                                        ImmersedBoundaryCellPopulation<DIM>& rCellPopulation,
                                                                        double intrinsicSpacingSquared)
{
    // Get index and number of nodes of current element
    unsigned elem_idx = rElement.GetIndex();
    unsigned num_nodes = rElement.GetNumNodes();

    // Set the vector that represents the force on each node
    c_vector<double, DIM> force = mStrength * mDirection;

    if (elem_idx == mElementIndex)
    {
        // Add the contributions of springs adjacent to each node
        for (unsigned node_idx = 0; node_idx < num_nodes; node_idx++) {

            // Add the aggregate force contribution to the node
            rElement.GetNode(node_idx)->AddAppliedForceContribution(force);
        }
    }
}

template<unsigned DIM>
void SingleCellMigrationForce<DIM>::OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
    *rParamsFile << "\t\t\t<ElementIndex>" << mElementIndex << "</ElementIndex>\n";

    // Call method on direct parent class
    AbstractImmersedBoundaryForce<DIM>::OutputImmersedBoundaryForceParameters(rParamsFile);
}

template <unsigned DIM>
void SingleCellMigrationForce<DIM>::SetElementIndex(unsigned element_index)
{
    mElementIndex = element_index;
}

template <unsigned DIM>
unsigned SingleCellMigrationForce<DIM>::GetElementIndex() const
{
    return mElementIndex;
}

template <unsigned DIM>
void SingleCellMigrationForce<DIM>::SetStrength(double strength)
{
	mStrength = strength;
}

template<unsigned DIM>
double SingleCellMigrationForce<DIM>::GetStrength() const
{
	return mStrength;
}

template <unsigned DIM>
void SingleCellMigrationForce<DIM>::SetDirection(c_vector<double, DIM> direction)
{
	// normalize the direction
	direction *= (1.0 / norm_2(direction));

	mDirection = direction;
}

template<unsigned DIM>
c_vector<double, DIM> SingleCellMigrationForce<DIM>::GetDirection() const
{
	return mDirection;
}

// Explicit Instantiation
template class SingleCellMigrationForce<1>;
template class SingleCellMigrationForce<2>;
template class SingleCellMigrationForce<3>;

// Serializaion for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SingleCellMigrationForce)
