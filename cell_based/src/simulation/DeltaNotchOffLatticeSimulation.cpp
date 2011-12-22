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

#include "DeltaNotchOffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CellwiseData.hpp"
#include "DeltaNotchCellCycleModel.hpp"

template<unsigned DIM>
DeltaNotchOffLatticeSimulation<DIM>::DeltaNotchOffLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                                                  bool deleteCellPopulationInDestructor,
                                                                  bool initialiseCells)
    : OffLatticeSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells)
{
}

template<unsigned DIM>
DeltaNotchOffLatticeSimulation<DIM>::~DeltaNotchOffLatticeSimulation()
{
}

template<unsigned DIM>
void DeltaNotchOffLatticeSimulation<DIM>::PostSolve()
{
    this->mrCellPopulation.Update();
    CellwiseData<DIM>::Instance()->ReallocateMemory();

    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
    {
        // Create and initialize vector of mean Delta concentrations
        unsigned num_cells = this->mrCellPopulation.GetNumRealCells();
        std::vector<double> mean_delta(num_cells);
        std::vector<unsigned> num_neighbours(num_cells);
        for (unsigned i=0; i<num_cells; i++)
        {
            mean_delta[i] = 0.0;
            num_neighbours[i] = 0;
        }
        if (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&(this->mrCellPopulation)))
        {
            ///\todo needs testing

            // Iterate over all springs (i.e. neighbouring cells)
            MeshBasedCellPopulation<DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<DIM>*>(&(this->mrCellPopulation));
            for (typename MeshBasedCellPopulation<DIM>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
                 spring_iterator != p_static_cast_cell_population->SpringsEnd();
                 ++spring_iterator)
            {
                // Get the index of each node on this edge
                unsigned node_a_index = spring_iterator.GetNodeA()->GetIndex();
                unsigned node_b_index = spring_iterator.GetNodeB()->GetIndex();

                // Get the cell associated with each node
                CellPtr p_cell_a = this->mrCellPopulation.GetCellUsingLocationIndex(node_a_index);
                CellPtr p_cell_b = this->mrCellPopulation.GetCellUsingLocationIndex(node_b_index);

                // Get the Delta concentration at each cell
                double delta_a = static_cast<DeltaNotchCellCycleModel*>(p_cell_a->GetCellCycleModel())->GetDelta();
                double delta_b = static_cast<DeltaNotchCellCycleModel*>(p_cell_b->GetCellCycleModel())->GetDelta();

                // Add each cell's neighbour's Delta concentration to mean_delta
                mean_delta[node_a_index] += delta_b;
                mean_delta[node_b_index] += delta_a;

                // Increment the number of neighbours found for each cell
                num_neighbours[node_a_index]++;
                num_neighbours[node_b_index]++;
            }
        }
        else // assume a NodeBasedCellPopulation
        {
            // Get a set of all neighbouring node pairs
            for (unsigned node_a_index=0; node_a_index<this->mrCellPopulation.GetNumNodes(); node_a_index++)
            {
                for (unsigned node_b_index=0; node_b_index<this->mrCellPopulation.GetNumNodes(); node_b_index++)
                {
                     bool neighbours = (norm_2(this->mrCellPopulation.GetNode(node_a_index)->rGetLocation() - this->mrCellPopulation.GetNode(node_b_index)->rGetLocation()) <= 1.5);
                     if ((node_a_index!= node_b_index) && (neighbours == true))
                     {
                        // Get the cell associated with each node
                        CellPtr p_cell_a = this->mrCellPopulation.GetCellUsingLocationIndex(node_a_index);
                        CellPtr p_cell_b = this->mrCellPopulation.GetCellUsingLocationIndex(node_b_index);

                        // Get the Delta concentration at each cell
                        //double delta_a = static_cast<DeltaNotchCellCycleModel*>(p_cell_a->GetCellCycleModel())->GetDelta();
                        double delta_b = static_cast<DeltaNotchCellCycleModel*>(p_cell_b->GetCellCycleModel())->GetDelta();

                        // Add each cell's neighbour's Delta concentration to mean_delta
                        mean_delta[node_a_index] += delta_b;

                        // Increment the number of neighbours found for each cell
                        num_neighbours[node_a_index]++;
                    }
                }
            }
        }

        // Make mean_delta a vector of mean (rather than total) concentrations and store in CellwiseData
        for (unsigned i=0; i<num_cells; i++)
        {
            if (num_neighbours[i] != 0)
            {
                mean_delta[i] /= num_neighbours[i];
            }
        }

        // Loop over cells
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
             cell_iter != this->mrCellPopulation.End();
             ++cell_iter)
        {
            // Also store this cell's Delta and Notch concentrations
            unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            double this_delta = static_cast<DeltaNotchCellCycleModel*>(cell_iter->GetCellCycleModel())->GetDelta();
            double this_notch = static_cast<DeltaNotchCellCycleModel*>(cell_iter->GetCellCycleModel())->GetNotch();
            CellwiseData<DIM>::Instance()->SetValue(mean_delta[node_index], node_index, 0);
            CellwiseData<DIM>::Instance()->SetValue(this_delta, node_index, 1);
            CellwiseData<DIM>::Instance()->SetValue(this_notch, node_index, 2);
        }
    }
    else // assume a VertexBasedCellPopulation
    {
        // Loop over cells
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mrCellPopulation.Begin();
             cell_iter != this->mrCellPopulation.End();
             ++cell_iter)
        {
            // Get the vertex element corresponding to this cell
            unsigned element_index = this->mrCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            VertexElement<DIM,DIM>* p_element = static_cast<VertexBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->GetElement(element_index);

            // Get a set of neighbouring element indices
            std::set<unsigned> neighbour_indices;
            for (unsigned local_index=0; local_index<p_element->GetNumNodes(); local_index++)
            {
                Node<DIM>* p_node = p_element->GetNode(local_index);
                std::set<unsigned> elements = p_node->rGetContainingElementIndices();
                std::set<unsigned> all_elements;
                std::set_union(neighbour_indices.begin(), neighbour_indices.end(),
                               elements.begin(), elements.end(),
                               std::inserter(all_elements, all_elements.begin()));

                neighbour_indices = all_elements;
            }
            neighbour_indices.erase(element_index);
            unsigned num_neighbours = neighbour_indices.size();

            // Compute mean Delta concentration, averaged over neighbours
            double mean_delta = 0.0;
            for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                 iter != neighbour_indices.end();
                 ++iter)
            {
                CellPtr p_cell = this->mrCellPopulation.GetCellUsingLocationIndex(*iter);
                double delta = static_cast<DeltaNotchCellCycleModel*>(p_cell->GetCellCycleModel())->GetDelta();
                if (num_neighbours != 0)
                {
                    mean_delta += delta/num_neighbours;
                }
            }

            // Store this value in CellwiseData
            CellwiseData<DIM>::Instance()->SetValue(mean_delta, p_element->GetIndex(), 0);

            // Also store this cell's Delta and Notch concentrations
            double this_delta = static_cast<DeltaNotchCellCycleModel*>(cell_iter->GetCellCycleModel())->GetDelta();
            double this_notch = static_cast<DeltaNotchCellCycleModel*>(cell_iter->GetCellCycleModel())->GetNotch();
            CellwiseData<DIM>::Instance()->SetValue(this_delta, p_element->GetIndex(), 1);
            CellwiseData<DIM>::Instance()->SetValue(this_notch, p_element->GetIndex(), 2);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DeltaNotchOffLatticeSimulation<1>;
template class DeltaNotchOffLatticeSimulation<2>;
template class DeltaNotchOffLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DeltaNotchOffLatticeSimulation)
