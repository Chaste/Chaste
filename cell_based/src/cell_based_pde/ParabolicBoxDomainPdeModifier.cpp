/*

Copyright (c) 2005-2024, University of Oxford.
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

#include "ParabolicBoxDomainPdeModifier.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "VtkMeshWriter.hpp"
#include "MutableMesh.hpp"

template<unsigned DIM>
ParabolicBoxDomainPdeModifier<DIM>::ParabolicBoxDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde,
                                                                  boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition,
                                                                  bool isNeumannBoundaryCondition,
                                                                  boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid,
                                                                  double stepSize,
                                                                  Vec solution)
    : AbstractBoxDomainPdeModifier<DIM>(pPde,
                                        pBoundaryCondition,
                                        isNeumannBoundaryCondition,
                                        pMeshCuboid,
                                        stepSize,
                                        solution),
      mMoveSolutionWithCells(false)
{
}

template<unsigned DIM>
ParabolicBoxDomainPdeModifier<DIM>::~ParabolicBoxDomainPdeModifier()
{
}

template<unsigned DIM>
void ParabolicBoxDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer(rCellPopulation);

    this->UpdateCellPdeElementMap(rCellPopulation);

    // When using a PDE mesh which doesn't coincide with the cells, we must set up the source terms before solving the PDE.
    // Pass in already updated CellPdeElementMap to speed up finding cells.
    this->SetUpSourceTermsForAveragedSourcePde(this->mpFeMesh, &this->mCellPdeElementMap);

    // Use SimpleLinearParabolicSolver as averaged Source PDE
    SimpleLinearParabolicSolver<DIM,DIM> solver(this->mpFeMesh,
                                                boost::static_pointer_cast<AbstractLinearParabolicPde<DIM,DIM> >(this->GetPde()).get(),
                                                p_bcc.get());

    ///\todo Investigate more than one PDE time step per spatial step
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();
    double dt = p_simulation_time->GetTimeStep();
    solver.SetTimes(current_time,current_time + dt);
    solver.SetTimeStep(dt);


    // Use previous solution as the initial condition
    Vec previous_solution = this->mSolution;
    if (mMoveSolutionWithCells)
    {
        // interpolate solution from cells movement onto fe mesh
        previous_solution = InterpolateSolutionFromCellMovement(rCellPopulation);
    }
    solver.SetInitialCondition(previous_solution);

    // Note that the linear solver creates a vector, so we have to keep a handle on the old one
    // in order to destroy it
    this->mSolution = solver.Solve();
    PetscTools::Destroy(previous_solution);

    // Now copy solution to cells
    this->UpdateCellData(rCellPopulation);

    // Finally, if needed store the locations of cells to be used as old loactions in the next timestep
    if (mMoveSolutionWithCells)
    {
        /*
         * If required, store the current locations of cell centres. Note that we need to
         * use a std::map between cells and locations, rather than (say) a std::vector with
         * location indices corresponding to cells, since once we call UpdateCellLocations()
         * the location index of each cell may change. This is especially true in the case
         * of a CaBasedCellPopulation.
         */
        mOldCellLocations.clear();
        for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
            mOldCellLocations[*cell_iter] = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        }
    } 
}

template<unsigned DIM>
void ParabolicBoxDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractBoxDomainPdeModifier<DIM>::SetupSolve(rCellPopulation,outputDirectory);

    // Copy the cell data to mSolution (this is the initial condition)
    SetupInitialSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);

    // If needed store the locations of cells to be used as old locations at the next timestep.
    if (mMoveSolutionWithCells)
    {
        /*
         * If required, store the current locations of cell centres. Note that we need to
         * use a std::map between cells and locations, rather than (say) a std::vector with
         * location indices corresponding to cells, since once we call UpdateCellLocations()
         * the location index of each cell may change. This is especially true in the case
         * of a CaBasedCellPopulation.
         */
        mOldCellLocations.clear();
        for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
        {
            mOldCellLocations[*cell_iter] = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        }
    }
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ParabolicBoxDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    this->ConstructBoundaryConditionsContainerHelper(rCellPopulation,p_bcc);
   
    return p_bcc;
}

template<unsigned DIM>
void ParabolicBoxDomainPdeModifier<DIM>::SetupInitialSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Specify homogeneous initial conditions based upon the values stored in CellData.
    // Note need all the CellDataValues to be the same.

    double initial_condition = rCellPopulation.Begin()->GetCellData()->GetItem(this->mDependentVariableName);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        double initial_condition_at_cell = cell_iter->GetCellData()->GetItem(this->mDependentVariableName);
        UNUSED_OPT(initial_condition_at_cell);
        assert(fabs(initial_condition_at_cell - initial_condition)<1e-12);
    }

    // Initialise mSolution
    this->mSolution = PetscTools::CreateAndSetVec(this->mpFeMesh->GetNumNodes(), initial_condition);
}

template<unsigned DIM>
Vec ParabolicBoxDomainPdeModifier<DIM>::InterpolateSolutionFromCellMovement(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Store the PDE solution in an accessible form
    ReplicatableVector solution_repl(this->mSolution);

    Vec interpolated_solution = PetscTools::CreateAndSetVec(this->mpFeMesh->GetNumNodes(), 0.0);

    // Store max radius so can speed up interpolation.
    // TODO replace this with more general exclusion based on dimensions of tissue.
    double max_radius = 0.0;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        double radius = norm_2(r_position_of_cell.rGetLocation());
        
        if (max_radius < radius)
        {
            max_radius = radius;
        }                
    }

    // Create mesh from cell centres so can interpolate tissue velocity onto mpFeMesh
    std::vector<Node<DIM>*> temp_nodes;
    std::vector<c_vector<double,DIM>> cell_displacements;
    unsigned cell_index = 0;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {
        c_vector<double,DIM> position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        cell_index++;

        // Only use cells that were in both this and the previous timestep.
        if (mOldCellLocations.find(*cell_iter) != mOldCellLocations.end())
        {
            //PRINT_VARIABLE(mOldCellLocations[*cell_iter]);
            c_vector<double,DIM> displacement = position_of_cell - mOldCellLocations[*cell_iter];
       
            cell_displacements.push_back(displacement);
            temp_nodes.push_back(new Node<DIM>(cell_index, position_of_cell));
        }
    }
    MutableMesh<DIM,DIM> cell_mesh(temp_nodes);

    // Make the deformed mesh. Based on the displacement of cells. 
    TetrahedralMesh<DIM, DIM>* p_deformed_mesh = new TetrahedralMesh<DIM,DIM>();
    this->GenerateAndReturnFeMesh(this->mpMeshCuboid,this->mStepSize,p_deformed_mesh);
    
    for (typename TetrahedralMesh<DIM, DIM>::NodeIterator node_iter = p_deformed_mesh->GetNodeIteratorBegin();
         node_iter != p_deformed_mesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        c_vector<double, DIM> node_location = node_iter->rGetLocation();

        
        c_vector<double, DIM> new_node_location = node_location;
        
        if (norm_2(node_location) <= max_radius)
        {
            // Find the element in the cell mesh that contains this node.
            try
            {
                unsigned elem_index = cell_mesh.GetContainingElementIndex(node_location, false);
        
                // Now do the interpolation
                Element<DIM,DIM>* p_element = cell_mesh.GetElement(elem_index);
                c_vector<double,DIM+1> weights = p_element->CalculateInterpolationWeights(node_location);

                c_vector<double,DIM> interpolated_cell_displacement = zero_vector<double>(DIM);
                for (unsigned i=0; i<DIM+1; i++)
                {
                    c_vector<double,DIM> nodal_value = cell_displacements[p_element->GetNodeGlobalIndex(i)];
                    interpolated_cell_displacement += nodal_value * weights(i);
                }
                new_node_location = node_location + interpolated_cell_displacement;
                node_iter->rGetModifiableLocation() = new_node_location;
            }
            catch (Exception&) // not_in_mesh
            {
                //Don't do anything as these FE nodes are outside of the Cell mesh.               
            }
        }   
    }

    // Loop over nodes of the mpFeMesh and get solution values from the deformed mesh
    for (typename TetrahedralMesh<DIM,DIM>::NodeIterator node_iter = this->mpFeMesh->GetNodeIteratorBegin();
            node_iter != this->mpFeMesh->GetNodeIteratorEnd();
            ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();

        const ChastePoint<DIM>& node_location = node_iter->rGetLocation();

        // Find the element in the deformed mesh that contains this FE node.
        std::set<unsigned> elements_to_check = node_iter->rGetContainingElementIndices();
        try
        { 
            unsigned elem_index = p_deformed_mesh->GetContainingElementIndex(node_location, false, elements_to_check);
       
            // Now do the interpolation
            Element<DIM,DIM>* p_element = p_deformed_mesh->GetElement(elem_index);
            c_vector<double,DIM+1> weights;

            weights = p_element->CalculateInterpolationWeights(node_location);
        
            double solution_at_node = 0.0;
            for (unsigned i=0; i<DIM+1; i++)
            {
                double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(i)];
                solution_at_node += nodal_value * weights(i);
            }
            PetscVecTools::SetElement(interpolated_solution, node_index, solution_at_node);
        }
        catch (Exception &e)
        {
            // Handy debug code to work out why the node is not in any element

            // output the cell mesh
            std::ostringstream time_string1;
            time_string1 << SimulationTime::Instance()->GetTimeStepsElapsed();
            std::string results_file1 = "cell_mesh_" + this->mDependentVariableName + "_" + time_string1.str();
            VtkMeshWriter<DIM, DIM>* p_vtk_mesh_writer1 = new VtkMeshWriter<DIM, DIM>(this->mOutputDirectory, results_file1, false);
            p_vtk_mesh_writer1->AddPointData(this->mDependentVariableName, cell_displacements);
            p_vtk_mesh_writer1->WriteFilesUsingMesh(cell_mesh);
            delete p_vtk_mesh_writer1; 

            // output the deformed mesh
            std::ostringstream time_string;
            time_string << SimulationTime::Instance()->GetTimeStepsElapsed();
            std::string results_file = "deformed_mesh_" + this->mDependentVariableName + "_" + time_string.str();
            VtkMeshWriter<DIM, DIM>* p_vtk_mesh_writer = new VtkMeshWriter<DIM, DIM>(this->mOutputDirectory, results_file, false);
            p_vtk_mesh_writer->WriteFilesUsingMesh(*p_deformed_mesh);
            delete p_vtk_mesh_writer;  

            assert(0);
        }
    }   

    // Tidy Up
    delete p_deformed_mesh;

    return interpolated_solution;
}

template<unsigned DIM>
void ParabolicBoxDomainPdeModifier<DIM>::SetMoveSolutionWithCells(bool moveSolutionWithCells)
{
    mMoveSolutionWithCells = moveSolutionWithCells;
}

template<unsigned DIM>
bool ParabolicBoxDomainPdeModifier<DIM>::GetMoveSolutionWithCells()
{
    return mMoveSolutionWithCells;
}

template<unsigned DIM>
void ParabolicBoxDomainPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractBoxDomainPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ParabolicBoxDomainPdeModifier<1>;
template class ParabolicBoxDomainPdeModifier<2>;
template class ParabolicBoxDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ParabolicBoxDomainPdeModifier)
