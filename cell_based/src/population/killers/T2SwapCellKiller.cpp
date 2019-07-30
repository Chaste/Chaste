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

#include "T2SwapCellKiller.hpp"

template<unsigned DIM>
T2SwapCellKiller<DIM>::T2SwapCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
    // Throw an exception if the population is not a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("A T2SwapCellKiller should only be used together with a VertexBasedCellPopulation.");
    }
}

template<unsigned DIM>
void T2SwapCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    /*
     * This killer is different to other killers: it does not only check and label
     * cells for apoptosis or death, it actually carries out vertex rearrangements
     * and removes elements from the vertex mesh.
     *
     * We start with carrying out T2 swaps. Get the mesh and an element map.
     * The static_cast will work since we already know it's a VertexBasedCellPopulation.
     */
    MutableVertexMesh<DIM,DIM>& mesh = static_cast<MutableVertexMesh<DIM,DIM>&>(this->mpCellPopulation->rGetMesh());
    VertexBasedCellPopulation<DIM>* p_vertex_population = static_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation);
    VertexElementMap element_map(mesh.GetNumAllElements());

    bool recheck_mesh = true;
    while (recheck_mesh == true)
    {
        // Note that whenever we call CheckForT2Swaps(), the element indices must run from zero up to mElements.size()-1
        recheck_mesh = mesh.CheckForT2Swaps(element_map);
        /*
         * There might have maximally one T2 swap happened above, where a vertex element was removed from the
         * mesh but the associated cell is still there. Here we check whether a new cell
         * underwent a T2 swap and label it as dead as well as record its location and ID.
         */
        for (unsigned elem_index = 0; elem_index < element_map.Size(); elem_index++)
        {
            CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(elem_index);
            if (element_map.IsDeleted(elem_index) && !(p_cell->IsDead()))
            {
                p_vertex_population->AddLocationOfT2Swap(mesh.GetLastT2SwapLocation());
                p_vertex_population->AddCellIdOfT2Swap(p_cell->GetCellId());
                p_cell->Kill();

                // There can't have been more than one new cell death, so leave the for loop here.
                break;
            }
        }
    }

}

template<unsigned DIM>
void T2SwapCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class T2SwapCellKiller<1>;
template class T2SwapCellKiller<2>;
template class T2SwapCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(T2SwapCellKiller)
