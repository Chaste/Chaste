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

#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "CellLabel.hpp"

template<unsigned DIM>
NagaiHondaDifferentialAdhesionForce<DIM>::NagaiHondaDifferentialAdhesionForce()
    : NagaiHondaForce<DIM>(),
      mNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(1.0),
      mNagaiHondaLabelledCellCellAdhesionEnergyParameter(1.0),
      mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(1.0)
{
}

template<unsigned DIM>
NagaiHondaDifferentialAdhesionForce<DIM>::~NagaiHondaDifferentialAdhesionForce()
{
}

template<unsigned DIM>
double NagaiHondaDifferentialAdhesionForce<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA,
                                                                      Node<DIM>* pNodeB,
                                                                      VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        unsigned element_index = *(shared_elements.begin());

        // Get cell associated with this element
        CellPtr p_cell = rVertexCellPopulation.GetCellUsingLocationIndex(element_index);

        if (p_cell->template HasCellProperty<CellLabel>())
        {
            // This cell is labelled
            return this->GetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter();
        }
        else
        {
            // This cell is not labelled
            return this->GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
        }
    }
    else
    {
        // Work out the number of labelled cells: 0,1 or 2
        unsigned num_labelled_cells = 0;
        for (std::set<unsigned>::iterator iter = shared_elements.begin();
             iter != shared_elements.end();
             ++iter)
        {
            unsigned element_index = *(iter);

            // Get cell associated with this element
            CellPtr p_cell = rVertexCellPopulation.GetCellUsingLocationIndex(element_index);

            if (p_cell->template HasCellProperty<CellLabel>())
            {
                num_labelled_cells++;
            }
        }

        if (num_labelled_cells == 2)
        {
            // Both cells are labelled
            return this->GetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter();
        }
        else if (num_labelled_cells == 1)
        {
            // One cell is labelled
            return this->GetNagaiHondaLabelledCellCellAdhesionEnergyParameter();
        }
        else
        {
            // Neither cell is labelled
            assert(num_labelled_cells == 0);
            return this->GetNagaiHondaCellCellAdhesionEnergyParameter();
        }
    }
}

template<unsigned DIM>
double NagaiHondaDifferentialAdhesionForce<DIM>::GetNagaiHondaLabelledCellCellAdhesionEnergyParameter()
{
    return mNagaiHondaLabelledCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double NagaiHondaDifferentialAdhesionForce<DIM>::GetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter()
{
    return mNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double NagaiHondaDifferentialAdhesionForce<DIM>::GetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaDifferentialAdhesionForce<DIM>::SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(double labelledCellCellAdhesionEnergyParameter)
{
    mNagaiHondaLabelledCellCellAdhesionEnergyParameter = labelledCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaDifferentialAdhesionForce<DIM>::SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(double labelledCellLabelledCellAdhesionEnergyParameter)
{
    mNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter = labelledCellLabelledCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaDifferentialAdhesionForce<DIM>::SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(double labelledCellBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter = labelledCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void NagaiHondaDifferentialAdhesionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Output member variables
    *rParamsFile << "\t\t\t<NagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter>" << mNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter << "</NagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter> \n";
    *rParamsFile << "\t\t\t<NagaiHondaLabelledCellCellAdhesionEnergyParameter>" << mNagaiHondaLabelledCellCellAdhesionEnergyParameter << "</NagaiHondaLabelledCellCellAdhesionEnergyParameter> \n";
    *rParamsFile << "\t\t\t<NagaiHondaLabelledCellBoundaryAdhesionEnergyParameter>" << mNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter << "</NagaiHondaLabelledCellBoundaryAdhesionEnergyParameter> \n";

    // Call method on direct parent class
    NagaiHondaForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class NagaiHondaDifferentialAdhesionForce<1>;
template class NagaiHondaDifferentialAdhesionForce<2>;
template class NagaiHondaDifferentialAdhesionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaDifferentialAdhesionForce)
