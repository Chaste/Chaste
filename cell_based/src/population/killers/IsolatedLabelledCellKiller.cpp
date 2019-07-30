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

#include "IsolatedLabelledCellKiller.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellLabel.hpp"

template<unsigned DIM>
IsolatedLabelledCellKiller<DIM>::IsolatedLabelledCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(pCellPopulation) == nullptr)
    {
        EXCEPTION("IsolatedLabelledCellKiller only works with a VertexBasedCellPopulation.");
    }
}

template<unsigned DIM>
void IsolatedLabelledCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    MutableVertexMesh<DIM, DIM>& vertex_mesh = static_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation)->rGetMesh();

    unsigned num_labelled_cells = this->mpCellPopulation->GetCellPropertyRegistry()->template Get<CellLabel>()->GetCellCount();

    // If there is more than one labelled cell...
    if (num_labelled_cells > 1)
    {
        // Iterate over cell population
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
             cell_iter != this->mpCellPopulation->End();
             ++cell_iter)
        {
            // Only consider cells with the CellLabel property
            if (cell_iter->template HasCellProperty<CellLabel>())
            {
                // Get the element index corresponding to this cell
                unsigned elem_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

                // Get the set of neighbouring element indices
                std::set<unsigned> neighbouring_elem_indices = vertex_mesh.GetNeighbouringElementIndices(elem_index);

                // Check if any of the corresponding cells have the CellLabel property...
                unsigned num_labelled_neighbours = 0;
                for (std::set<unsigned>::iterator elem_iter = neighbouring_elem_indices.begin();
                     elem_iter != neighbouring_elem_indices.end();
                     ++elem_iter)
                {
                    if (this->mpCellPopulation->GetCellUsingLocationIndex(*elem_iter)->template HasCellProperty<CellLabel>())
                    {
                        num_labelled_neighbours++;
                    }
                }

                // ...and if none do, then kill this cell
                if (num_labelled_neighbours == 0)
                {
                    cell_iter->Kill();
                }
            }
        }
    }
}

template<unsigned DIM>
void IsolatedLabelledCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    // There are no member variables, so just call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class IsolatedLabelledCellKiller<1>;
template class IsolatedLabelledCellKiller<2>;
template class IsolatedLabelledCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(IsolatedLabelledCellKiller)
