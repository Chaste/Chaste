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
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"

#include "ReplicatableVector.hpp"
#include "OdeLinearSystemSolver.hpp"

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
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::UpdateNodeLocations(double dt)
{
    // Declare solver and give the size of the system and timestep
    unsigned system_size = this->GetNumNodes()*DIM;

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
        unsigned global_node_index = this->GetLocationIndexUsingCell((*cell_iter));

        // Get the local index using the mesh
        unsigned node_index = this->rGetMesh().SolveNodeMapping(global_node_index);

        Node<DIM>* p_node_i = this->GetNode(global_node_index);

        // Get the location of this node
        const c_vector<double, DIM>& r_node_i_location = p_node_i->rGetLocation();

        // Get the radius of this cell
        double radius_of_cell_i = p_node_i->GetRadius();

        // Get damping constant for node
        double damping_const = this->GetDampingConstant(global_node_index);

        // loop over neighbours to add contribution

        // Get the set of node indices corresponding to this cell's neighbours
        std::set<unsigned> neighbouring_node_indices = this->GetNeighbouringNodeIndices(global_node_index);

        for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
             iter != neighbouring_node_indices.end();
             ++iter)
        {
            unsigned neighbour_node_global_index = *iter;

            unsigned neighbour_node_index = this->rGetMesh().SolveNodeMapping(neighbour_node_global_index);

            // Calculate Aij
            double Aij = 0.0;

            Node<DIM>* p_node_j = this->GetNode(neighbour_node_global_index);

            // Get the location of this node
            const c_vector<double, DIM>& r_node_j_location = p_node_j->rGetLocation();

            // Get the unit vector parallel to the line joining the two nodes (assuming no periodicities etc.)
            c_vector<double, DIM> unit_vector = r_node_j_location - r_node_i_location;

            // Calculate the distance between the two nodes
            double dij = norm_2(unit_vector);

            unit_vector /= dij;

            // Get the radius of the cell corresponding to this node
            double radius_of_cell_j = p_node_j->GetRadius();

            if (dij < radius_of_cell_i + radius_of_cell_j)
            {
                // ...then compute the adhesion force and add it to the vector of forces...
                double xij = 0.5*(radius_of_cell_i*radius_of_cell_i - radius_of_cell_j*radius_of_cell_j + dij*dij)/dij;

                Aij = M_PI*(radius_of_cell_i*radius_of_cell_i - xij*xij);

                // This is contribution from the sum term in (A7)
                for (unsigned i=0; i<DIM; i++)
                {
                    PetscMatTools::AddToElement(r_matrix, DIM*neighbour_node_index+i, DIM*neighbour_node_index+i, -damping_const*Aij);
                    PetscMatTools::AddToElement(r_matrix, DIM*node_index+i, DIM*node_index+i, damping_const*Aij);
                }
            }
        }

        // This is the standard contribution (i.e. not in the sum) in (A7)
        for (unsigned i=0; i<DIM; i++)
        {
            PetscMatTools::AddToElement(r_matrix, DIM*node_index+i, DIM*node_index+i, damping_const);
        }

        // Add current positions to initial_conditions and RHS vector

        // Note that we define these vectors before setting them as otherwise the profiling build will break (see #2367)
        c_vector<double, DIM> current_location;
        c_vector<double, DIM> forces;
        current_location = this->GetNode(global_node_index)->rGetLocation();
        forces = this->GetNode(global_node_index)->rGetAppliedForce();

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
        unsigned global_node_index = this->GetLocationIndexUsingCell((*cell_iter));

        unsigned node_index = this->rGetMesh().SolveNodeMapping(global_node_index);

        c_vector<double, DIM> new_node_location;

        // Get new node location
        for (unsigned i=0; i<DIM; i++)
        {
            new_node_location(i) = soln_next_timestep_repl[DIM*node_index+i];
        }

        // Create ChastePoint for new node location
        ChastePoint<DIM> new_point(new_node_location);

        // Move the node
        this->SetNode(global_node_index, new_point);
    }

    // Tidy up
    PetscTools::Destroy(initial_condition);
}

template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    // Currently no specific parameters to output all come from parent classes

    // Call method on direct parent class
    NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

// Explicit instantiation
template class NodeBasedCellPopulationWithBuskeUpdate<1>;
template class NodeBasedCellPopulationWithBuskeUpdate<2>;
template class NodeBasedCellPopulationWithBuskeUpdate<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithBuskeUpdate)
