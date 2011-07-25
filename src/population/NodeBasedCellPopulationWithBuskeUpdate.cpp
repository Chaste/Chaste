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
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"

// The header file below must be included in any file that uses Petsc
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

template<unsigned DIM>
NodeBasedCellPopulationWithBuskeUpdate<DIM>::NodeBasedCellPopulationWithBuskeUpdate(NodesOnlyMesh<DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh)
    : NodeBasedCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh)
{
}

template<unsigned DIM>
NodeBasedCellPopulationWithBuskeUpdate<DIM>::NodeBasedCellPopulationWithBuskeUpdate(NodesOnlyMesh<DIM>& rMesh)
    : NodeBasedCellPopulation<DIM>(rMesh)
{
    // No Validate() because the cells are not associated with the cell population yet in archiving
}

template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::UpdateNodeLocations(const std::vector< c_vector<double, DIM> >& rNodeForces, double dt)
{


	// Declare solver and give the size of the system and timestep
	unsigned system_size = rNodeForces.size()*DIM;

	OdeLinearSystemSolver solver(system_size, dt);

	// Set up the matrix
	Mat& r_matrix = solver.rGetLhsMatrix();

	// Initial condition
	Vec initial_condition = PetscTools::CreateAndSetVec(system_size, 0.0);

	// Then an rGetForceVector for RHS
	Vec& r_vector = solver.rGetForceVector();

	// Iterate over all nodes associated with real cells to construct the matrix A.
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
		 cell_iter != this->End();
		 ++cell_iter)
	{
		// Get index of node associated with cell
		unsigned node_index = this->mCellLocationMap[(*cell_iter).get()];

		// Get damping constant for node
		double damping_const = this->GetDampingConstant(node_index);


		// This is currently just the identity matrix will need to add the area contribution see #1770
		for (unsigned i=0; i<DIM; i++)
		{
			PetscMatTools::SetElement(r_matrix, DIM*node_index+i, DIM*node_index+i, damping_const);
		}

		// Add current positions to initial_conditions and RHS vector
		c_vector<double, DIM> current_location = this->GetNode(node_index)->rGetLocation();
		c_vector<double, DIM> forces = rNodeForces[node_index];

		for (unsigned i=0; i<DIM; i++)
		{
			PetscVecTools::SetElement(initial_condition, DIM*node_index+i, current_location(i));
			PetscVecTools::SetElement(r_vector, DIM*node_index+i, forces(i));
		}
	}
	PetscMatTools::Finalise(r_matrix);

	solver.SetInitialConditionVector(initial_condition);

	// Solve to get solution at next timestep
	Vec soln_next_timestep = solver.SolveOneTimeStep();

	ReplicatableVector soln_next_timestep_repl(soln_next_timestep);

	// Iterate over all nodes associated with real cells to update the node locations
	for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        // Get index of node associated with cell
        unsigned node_index = this->mCellLocationMap[(*cell_iter).get()];

        c_vector<double, DIM> new_node_location;

        // Get new node location
		for (unsigned i=0; i<DIM; i++)
		{
			new_node_location(i) = soln_next_timestep_repl[DIM*node_index+i];
		}

        // Create ChastePoint for new node location
        ChastePoint<DIM> new_point(new_node_location);

        // Move the node
        this->SetNode(node_index, new_point);
    }

	// Tidy up
	VecDestroy(initial_condition);
}

template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Currently no specific parameters to output all come from parent classes

    // Call method on direct parent class
	NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}



/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class NodeBasedCellPopulationWithBuskeUpdate<1>;
template class NodeBasedCellPopulationWithBuskeUpdate<2>;
template class NodeBasedCellPopulationWithBuskeUpdate<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithBuskeUpdate)
