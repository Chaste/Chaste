/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "NodeBasedCellPopulationWithParticles.hpp"
#include "VtkMeshWriter.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"

template<unsigned DIM>
NodeBasedCellPopulationWithParticles<DIM>::NodeBasedCellPopulationWithParticles(NodesOnlyMesh<DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh)
    : NodeBasedCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh, false)
{
    if (!locationIndices.empty())
    {
        // Create a set of node indices corresponding to particles
        std::set<unsigned> node_indices;
        std::set<unsigned> location_indices;
        std::set<unsigned> particle_indices;

        for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
                node_iter != rMesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_indices.insert(node_iter->GetIndex());
        }
        for (unsigned i=0; i<locationIndices.size(); i++)
        {
            location_indices.insert(locationIndices[i]);
        }

        std::set_difference(node_indices.begin(), node_indices.end(),
                            location_indices.begin(), location_indices.end(),
                            std::inserter(particle_indices, particle_indices.begin()));

        // This method finishes and then calls Validate()
        SetParticles(particle_indices);
    }
    else
    {
        for (typename NodesOnlyMesh<DIM>::NodeIterator node_iter = rMesh.GetNodeIteratorBegin();
                node_iter != rMesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            (*node_iter).SetIsParticle(false);
        }
        NodeBasedCellPopulationWithParticles::Validate();
    }
}

template<unsigned DIM>
NodeBasedCellPopulationWithParticles<DIM>::NodeBasedCellPopulationWithParticles(NodesOnlyMesh<DIM>& rMesh)
    : NodeBasedCellPopulation<DIM>(rMesh)
{
}

template<unsigned DIM>
bool NodeBasedCellPopulationWithParticles<DIM>::IsParticle(unsigned index)
{
    return this->GetNode(index)->IsParticle();
}

template<unsigned DIM>
std::set<unsigned> NodeBasedCellPopulationWithParticles<DIM>::GetParticleIndices()
{
    std::set<unsigned> particle_indices;

    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        if (node_iter->IsParticle())
        {
            particle_indices.insert(node_iter->GetIndex());
        }
    }

    return particle_indices;
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::SetParticles(const std::set<unsigned>& rParticleIndices)
{
    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        if (rParticleIndices.find(node_iter->GetIndex()) != rParticleIndices.end())
        {
            node_iter->SetIsParticle(true);
        }
    }
    NodeBasedCellPopulationWithParticles::Validate();
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::UpdateParticlePositions(double dt)
{
    // Initialise vector of forces on particles
    std::vector<c_vector<double, DIM> > drdt(this->GetNumNodes());
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i] = zero_vector<double>(DIM);
    }

    // Calculate forces on particles
    double damping_constant = this->GetDampingConstantNormal();
    for (unsigned i=0; i<drdt.size(); i++)
    {
        drdt[i]=this->GetNode(i)->rGetAppliedForce()/damping_constant;
    }

    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        if (node_iter->IsParticle())
        {
            ChastePoint<DIM> new_point(node_iter->rGetLocation() + dt*drdt[node_iter->GetIndex()]);
            node_iter->SetPoint(new_point);
        }
    }
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::UpdateParticlesAfterReMesh(NodeMap& rMap)
{
}

template<unsigned DIM>
CellPtr NodeBasedCellPopulationWithParticles<DIM>::AddCell(CellPtr pNewCell, const c_vector<double,DIM>& rCellDivisionVector, CellPtr pParentCell)
{
    assert(pNewCell);

    // Add new cell to cell population
    CellPtr p_created_cell = AbstractCentreBasedCellPopulation<DIM>::AddCell(pNewCell, rCellDivisionVector, pParentCell);
    assert(p_created_cell == pNewCell);

    // Then set the new cell radius in the NodesOnlyMesh
    unsigned node_index = this->GetLocationIndexUsingCell(p_created_cell);
    this->GetNode(node_index)->SetRadius(0.5);
    this->GetNode(node_index)->SetIsParticle(false);

    // Return pointer to new cell
    return p_created_cell;
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::Validate()
{
    std::map<unsigned, bool> validated_nodes;
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
            node_iter != this->mrMesh.GetNodeIteratorEnd();
            ++node_iter)
    {
        validated_nodes[node_iter->GetIndex()] = node_iter->IsParticle();
    }

    // Look through all of the cells and record what node they are associated with.
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->GetLocationIndexUsingCell((*cell_iter));

        // If the node attached to this cell is labelled as a particle, then throw an error
        if (this->GetNode(node_index)->IsParticle())
        {
            EXCEPTION("Node " << node_index << " is labelled as a particle and has a cell attached");
        }
        validated_nodes[node_index] = true;
    }

    for (std::map<unsigned, bool>::iterator map_iter = validated_nodes.begin();
            map_iter != validated_nodes.end();
            map_iter++)
    {
        if (!map_iter->second)
        {
            EXCEPTION("Node " << map_iter->first << " does not appear to be a particle or has a cell associated with it");
        }
    }
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::UpdateNodeLocations(double dt)
{
    // First update particle positions
    UpdateParticlePositions(dt);
    // Then call the base class method
    AbstractCentreBasedCellPopulation<DIM>::UpdateNodeLocations(dt);
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::WriteVtkResultsToFile(const std::string& rDirectory)
{
///\todo #2441 - change this method to make use of the CellWriter functionality
///(see MeshBasedCellPopulation::WriteVtkResultsToFile, for example)
#ifdef CHASTE_VTK
    unsigned num_timesteps = SimulationTime::Instance()->GetTimeStepsElapsed();
    std::stringstream time;
    time << num_timesteps;

    VtkMeshWriter<DIM, DIM> mesh_writer(rDirectory, "results_"+time.str(), false);

    unsigned num_cells = this->GetNumNodes();
    std::vector<double> particles(num_cells);
    std::vector<double> cell_types(num_cells);
    std::vector<double> cell_ancestors(num_cells);
    std::vector<double> cell_mutation_states(num_cells);
    std::vector<double> cell_ages(num_cells);
    std::vector<double> cell_cycle_phases(num_cells);
    std::vector<double> cell_radii(num_cells);

    ///\todo #1975 - deal with possibility of information stored in CellData

    // Loop over nodes
    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        Node<DIM>* p_node = this->GetNode(node_index);

        if (!this->IsParticle(node_index))
        {
            CellPtr cell_iter = this->GetCellUsingLocationIndex(node_index);
            if (this-> template HasWriter<CellAncestorWriter>())
            {
                double ancestor_index = (cell_iter->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)cell_iter->GetAncestor();
                cell_ancestors[node_index] = ancestor_index;
            }
            if (this-> template HasWriter<CellProliferativeTypesWriter>())
            {
                double cell_type = cell_iter->GetCellProliferativeType()->GetColour();
                cell_types[node_index] = cell_type;
            }
            if (this-> template HasWriter<CellMutationStatesCountWriter>())
            {
                double mutation_state = cell_iter->GetMutationState()->GetColour();
                cell_mutation_states[node_index] = mutation_state;
            }
            if (this-> template HasWriter<CellAgesWriter>())
            {
                double age = cell_iter->GetAge();
                cell_ages[node_index] = age;
            }
            if (this-> template HasWriter<CellProliferativePhasesWriter>())
            {
                double cycle_phase = cell_iter->GetCellCycleModel()->GetCurrentCellCyclePhase();
                cell_cycle_phases[node_index] = cycle_phase;
            }
            if (this-> template HasWriter<CellVolumesWriter>())
            {
                double cell_radius = p_node->GetRadius();
                cell_radii[node_index] = cell_radius;
            }
        }
        else
        {
            particles[node_index] = (double)(this->IsParticle(node_index));
            if (this-> template HasWriter<CellAncestorWriter>())
            {
                cell_ancestors[node_index] = -2.0;
            }
            if (this-> template HasWriter<CellProliferativeTypesWriter>())
            {
                cell_types[node_index] = -2.0;
            }
            if (this-> template HasWriter<CellMutationStatesCountWriter>())
            {
                cell_mutation_states[node_index] = -2.0;
            }
            if (this-> template HasWriter<CellAgesWriter>())
            {
                cell_ages[node_index] = -2.0;
            }
            if (this-> template HasWriter<CellProliferativePhasesWriter>())
            {
                cell_cycle_phases[node_index] = -2.0;
            }
        }
    }

    mesh_writer.AddPointData("Non-particles", particles);
    if (this-> template HasWriter<CellProliferativeTypesWriter>())
    {
        mesh_writer.AddPointData("Cell types", cell_types);
    }
    if (this-> template HasWriter<CellAncestorWriter>())
    {
        mesh_writer.AddPointData("Ancestors", cell_ancestors);
    }
    if (this-> template HasWriter<CellMutationStatesCountWriter>())
    {
        mesh_writer.AddPointData("Mutation states", cell_mutation_states);
    }
    if (this-> template HasWriter<CellAgesWriter>())
    {
        mesh_writer.AddPointData("Ages", cell_ages);
    }
    if (this-> template HasWriter<CellProliferativePhasesWriter>())
    {
        mesh_writer.AddPointData("Cycle phases", cell_cycle_phases);
    }
    if (this-> template HasWriter<CellVolumesWriter>())
    {
        mesh_writer.AddPointData("Cell radii", cell_radii);
    }
    ///\todo #1975 - deal with possibility of information stored in CellData

    mesh_writer.WriteFilesUsingMesh(static_cast<NodesOnlyMesh<DIM>& >((this->mrMesh)));

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << num_timesteps;
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#endif //CHASTE_VTK
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

// Explicit instantiation
template class NodeBasedCellPopulationWithParticles<1>;
template class NodeBasedCellPopulationWithParticles<2>;
template class NodeBasedCellPopulationWithParticles<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithParticles)
