/*

Copyright (C) University of Oxford, 2005-2012

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

#include "NodeBasedCellPopulationWithParticles.hpp"
#include "CellwiseData.hpp"
#include "VtkMeshWriter.hpp"

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

			for (unsigned i=0; i<this->GetNumNodes(); i++)
			{
				node_indices.insert(this->GetNode(i)->GetIndex());
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
			this->mIsParticle = std::vector<bool>(this->GetNumNodes(), false);
			NodeBasedCellPopulationWithParticles::Validate();
		}
}

template<unsigned DIM>
NodeBasedCellPopulationWithParticles<DIM>::NodeBasedCellPopulationWithParticles(NodesOnlyMesh<DIM>& rMesh)
    : NodeBasedCellPopulation<DIM>(rMesh)
{
}


template<unsigned DIM>
std::vector<bool>& NodeBasedCellPopulationWithParticles<DIM>::rGetParticles()
{
    return this->mIsParticle;
}

template<unsigned DIM>
bool NodeBasedCellPopulationWithParticles<DIM>::IsParticle(unsigned index)
{
    return this->mIsParticle[index];
}

template<unsigned DIM>
std::set<unsigned> NodeBasedCellPopulationWithParticles<DIM>::GetParticleIndices()
{
    std::set<unsigned> particle_indices;
    for (unsigned i=0; i<this->mIsParticle.size(); i++)
    {
        if (this->mIsParticle[i])
        {
            particle_indices.insert(i);
        }
    }
    return particle_indices;
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::SetParticles(const std::set<unsigned>& rParticleIndices)
{
    // Reinitialise all entries of mIsParticle to false
    this->mIsParticle = std::vector<bool>(this->mrMesh.GetNumNodes(), false);

    // Update mIsParticle
    for (std::set<unsigned>::iterator iter=rParticleIndices.begin(); iter!=rParticleIndices.end(); ++iter)
    {
        this->mIsParticle[*iter] = true;
    }

    NodeBasedCellPopulationWithParticles::Validate();
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::UpdateParticlePositions(const std::vector< c_vector<double, DIM> >& rNodeForces,
																		double dt)
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
    drdt[i]=rNodeForces[i]/damping_constant;
    }

    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
         node_iter != this->mrMesh.GetNodeIteratorEnd();
         ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        if (this->mIsParticle[node_index])
        {
            ChastePoint<DIM> new_point(node_iter->rGetLocation() + dt*drdt[node_index]);
            this->mrMesh.SetNode(node_index, new_point, false);
        }
    }
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::UpdateParticlesAfterReMesh(NodeMap& rMap)
{
    // Copy mIsParticle to a temporary vector
    std::vector<bool> particles_before_remesh = mIsParticle;

    // Reinitialise mIsParticle
    mIsParticle.clear();
    mIsParticle.resize(this->GetNumNodes());

    // Update mIsParticle using the node map
    for (unsigned old_index=0; old_index<rMap.Size(); old_index++)
    {
        if (!rMap.IsDeleted(old_index))
        {
            unsigned new_index = rMap.GetNewIndex(old_index);
            mIsParticle[new_index] = particles_before_remesh[old_index];
        }
    }
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
    this->mrMesh.SetCellRadius(node_index, 0.5);

    // Update size of mIsParticle if necessary
    if (this->GetNumNodes() > this->mIsParticle.size())
    {
        this->mIsParticle.resize(this->GetNumNodes());
        this->mIsParticle[node_index] = false;
    }

    // Return pointer to new cell
    return p_created_cell;
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::Validate()
{

    // Get a list of all the nodes that are particles
    std::vector<bool> validated_node = mIsParticle;
    assert(mIsParticle.size()==this->GetNumNodes());

    // Look through all of the cells and record what node they are associated with.
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        unsigned node_index = this->mCellLocationMap[(*cell_iter).get()];

        // If the node attached to this cell is labelled as a particle, then throw an error
        if (mIsParticle[node_index])
        {
            EXCEPTION("Node " << node_index << " is labelled as a particle and has a cell attached");
        }
        validated_node[node_index] = true;
    }

    for (unsigned i=0; i<validated_node.size(); i++)
    {
        if (!validated_node[i])
        {
            EXCEPTION("Node " << i << " does not appear to be a particle or has a cell associated with it");
        }
    }
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)
{
    // First update particle positions
    UpdateParticlePositions(rNodeForces, dt);
    // Then call the base class method
    AbstractCentreBasedCellPopulation<DIM>::UpdateNodeLocations(rNodeForces, dt);
}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::WriteVtkResultsToFile()
{
#ifdef CHASTE_VTK
    std::stringstream time;
    time << SimulationTime::Instance()->GetTimeStepsElapsed();
    VtkMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results_"+time.str(), false);

    unsigned num_nodes = this->GetNumNodes();
    std::vector<double> particles(num_nodes);
    std::vector<double> cell_types(num_nodes);
    std::vector<double> cell_ancestors(num_nodes);
    std::vector<double> cell_mutation_states(num_nodes);
    std::vector<double> cell_ages(num_nodes);
    std::vector<double> cell_cycle_phases(num_nodes);
    std::vector<double> cell_radii(num_nodes);
    std::vector<std::vector<double> > cellwise_data;

    if (CellwiseData<DIM>::Instance()->IsSetUp())
    {
        CellwiseData<DIM>* p_data = CellwiseData<DIM>::Instance();
        unsigned num_variables = p_data->GetNumVariables();
        for (unsigned var=0; var<num_variables; var++)
        {
            std::vector<double> cellwise_data_var(num_nodes);
            cellwise_data.push_back(cellwise_data_var);
        }
    }

    // Loop over nodes

    for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = this->mrMesh.GetNodeIteratorBegin();
            node_iter != this->mrMesh.GetNodeIteratorEnd();
            ++node_iter)
       {
        unsigned node_index = node_iter->GetIndex();

        if (!this->IsParticle(node_index))
        {
        	CellPtr cell_iter = this->GetCellUsingLocationIndex(node_index);
        	if (this->mOutputCellAncestors)
			{
				double ancestor_index = (cell_iter->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)cell_iter->GetAncestor();
				cell_ancestors[node_index] = ancestor_index;
			}
			if (this->mOutputCellProliferativeTypes)
			{
				double cell_type = cell_iter->GetCellCycleModel()->GetCellProliferativeType();
				cell_types[node_index] = cell_type;
			}
			if (this->mOutputCellMutationStates)
			{
				double mutation_state = cell_iter->GetMutationState()->GetColour();
				cell_mutation_states[node_index] = mutation_state;
			}
			if (this->mOutputCellAges)
			{
				double age = cell_iter->GetAge();
				cell_ages[node_index] = age;
			}
			if (this->mOutputCellCyclePhases)
			{
				double cycle_phase = cell_iter->GetCellCycleModel()->GetCurrentCellCyclePhase();
				cell_cycle_phases[node_index] = cycle_phase;
			}
			if (this->mOutputCellVolumes)
			{
				double cell_radius = this->mrMesh.GetCellRadius(node_index);
				cell_radii[node_index] = cell_radius;
			}
			if (CellwiseData<DIM>::Instance()->IsSetUp())
			{
				CellwiseData<DIM>* p_data = CellwiseData<DIM>::Instance();
				unsigned num_variables = p_data->GetNumVariables();
				for (unsigned var=0; var<num_variables; var++)
				{
					cellwise_data[var][node_index] = p_data->GetValue(cell_iter, var);
				}
			}
        }
        else
        {
        	particles[node_index] = (double)(this->IsParticle(node_index));
			if (this->mOutputCellAncestors)
			{
				cell_ancestors[node_index] = -2.0;
			}
			if (this->mOutputCellProliferativeTypes)
			{
				cell_types[node_index] = -2.0;
			}
			if (this->mOutputCellMutationStates)
			{
				cell_mutation_states[node_index] = -2.0;
			}
			if (this->mOutputCellAges)
			{
				cell_ages[node_index] = -2.0;
			}
			if (this->mOutputCellCyclePhases)
			{
				cell_cycle_phases[node_index] = -2.0;
			}
        }
       }

    mesh_writer.AddCellData("Non-particles", particles);
    if (this->mOutputCellProliferativeTypes)
    {
        mesh_writer.AddPointData("Cell types", cell_types);
    }
    if (this->mOutputCellAncestors)
    {
        mesh_writer.AddPointData("Ancestors", cell_ancestors);
    }
    if (this->mOutputCellMutationStates)
    {
        mesh_writer.AddPointData("Mutation states", cell_mutation_states);
    }
    if (this->mOutputCellAges)
    {
        mesh_writer.AddPointData("Ages", cell_ages);
    }
    if (this->mOutputCellCyclePhases)
    {
        mesh_writer.AddPointData("Cycle phases", cell_cycle_phases);
    }
    if (this->mOutputCellVolumes)
    {
        mesh_writer.AddPointData("Cell radii", cell_radii);
    }
    if (CellwiseData<DIM>::Instance()->IsSetUp())
    {
        for (unsigned var=0; var<cellwise_data.size(); var++)
        {
            std::stringstream data_name;
            data_name << "Cellwise data " << var;
            std::vector<double> cellwise_data_var = cellwise_data[var];
            mesh_writer.AddPointData(data_name.str(), cellwise_data_var);
        }
    }

    mesh_writer.WriteFilesUsingMesh(this->mrMesh);

    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
    *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
    *(this->mpVtkMetaFile) << SimulationTime::Instance()->GetTimeStepsElapsed();
    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
#endif //CHASTE_VTK

}

template<unsigned DIM>
void NodeBasedCellPopulationWithParticles<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class NodeBasedCellPopulationWithParticles<1>;
template class NodeBasedCellPopulationWithParticles<2>;
template class NodeBasedCellPopulationWithParticles<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithParticles)
