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

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "TetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "CellBasedParabolicPdeSolver.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "AveragedSourcePde.hpp"

template<unsigned DIM>
ParabolicGrowingDomainPdeModifier<DIM>::ParabolicGrowingDomainPdeModifier(ParabolicPdeAndBoundaryConditions<DIM>* pPdeAndBcs)
    : AbstractGrowingDomainPdeModifier<DIM>(),
      mpPdeAndBcs(pPdeAndBcs)
{
    assert(DIM==2);
}

template<unsigned DIM>
ParabolicGrowingDomainPdeModifier<DIM>::~ParabolicGrowingDomainPdeModifier()
{
    // It we have used this modifier, then we will have created a solution vector
    if (this->mSolution)
    {
        PetscTools::Destroy(this->mSolution);
    }
}


template<unsigned DIM>
void ParabolicGrowingDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    this->GenerateFeMesh(rCellPopulation);

    // Set up boundary conditions
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer();

    // Construct the solution vector from cell data (takes care of cells dividing);
    UpdateSolutionVector(rCellPopulation);

    // Use CellBasedParabolicPdeSolver as cell wise PDE
    CellBasedParabolicPdeSolver<DIM> solver(this->mpFeMesh, mpPdeAndBcs->GetPde(), p_bcc.get());

    ///\todo Investigate more than one PDE timestep per spatial step (#2687)
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();
    double dt = p_simulation_time->GetTimeStep();
    solver.SetTimes(current_time,current_time + dt);
    solver.SetTimeStep(dt);

    // Use previous solution as the initial condition
    solver.SetInitialCondition(this->mSolution);

    this->mSolution = solver.Solve();

    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ParabolicGrowingDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Temporarily cache the variable name until we create a ParaolicPdeAndBcs object and move it to the Abstract class
    ///\todo this comment doesn't quite make sense (#2687)
    this->mCachedDependentVariableName = mpPdeAndBcs->rGetDependentVariableName();

    // Cache the Output Directory
    this->mOutputDirectory = outputDirectory;

    // Setup a finite element mesh on which to save the initial condition
    this->GenerateFeMesh(rCellPopulation);

    // Copy the cell data to mSolution (this is the initial condition)
    UpdateSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ParabolicGrowingDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer()
{
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    if (mpPdeAndBcs->IsNeumannBoundaryCondition())
    {
        // Impose any Neumann boundary conditions
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
             elem_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
             ++elem_iter)
        {
            p_bcc->AddNeumannBoundaryCondition(*elem_iter, mpPdeAndBcs->GetBoundaryCondition());
        }
    }
    else
    {
        // Impose any Dirichlet boundary conditions
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
                         node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
                         ++node_iter)
        {
            p_bcc->AddDirichletBoundaryCondition(*node_iter, mpPdeAndBcs->GetBoundaryCondition());
        }
    }

    return p_bcc;
}


template<unsigned DIM>
void ParabolicGrowingDomainPdeModifier<DIM>::UpdateSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Clear (if it's not the first time) and resize the solution vector
    if (this->mSolution)
    {
        PetscTools::Destroy(this->mSolution);
    }
    this->mSolution = PetscTools::CreateAndSetVec(this->mpFeMesh->GetNumNodes(), 0.0);

    // Loop over nodes and get appropriate solution value from CellData
    for (typename TetrahedralMesh<DIM,DIM>::NodeIterator node_iter = this->mpFeMesh->GetNodeIteratorBegin();
          node_iter != this->mpFeMesh->GetNodeIteratorEnd();
          ++node_iter)
    {
        unsigned node_index = node_iter->GetIndex();
        double solution_at_node;

        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) != NULL) // Intel compiler wants the "!= NULL"
        {
            // Cells correspond to nodes in the Center of the vertex element
            // nodes on vertices have averaged values from containing cells

            unsigned num_vertex_nodes = rCellPopulation.GetNumNodes();
            if (node_index >= num_vertex_nodes)
            {
                // Offset to relate elements in vertex mesh to nodes in tetrahedral mesh
                assert(node_index-num_vertex_nodes < num_vertex_nodes);

                CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index-num_vertex_nodes);
                solution_at_node = p_cell->GetCellData()->GetItem(mpPdeAndBcs->rGetDependentVariableName());
            }
            else
            {
                ///\todo Work out a better way to do the nodes not associated with cells (#2687)
                if (node_iter->IsBoundaryNode() && !mpPdeAndBcs->IsNeumannBoundaryCondition())
                {
                    // We need to impose the Dirichlet boundaries again here as not represented in cell data
                    solution_at_node = mpPdeAndBcs->GetBoundaryCondition()->GetValue(node_iter->rGetLocation());
                }
                else
                {
                    assert(node_index<num_vertex_nodes);
                    Node<DIM>* p_vertex_node = rCellPopulation.rGetMesh().GetNode(node_index);

                    // Average over data from containing elements (cells)
                    std::set<unsigned> containing_elelments  = p_vertex_node->rGetContainingElementIndices();

                    solution_at_node = 0.0;

                    for (std::set<unsigned>::iterator index_iter = containing_elelments.begin();
                          index_iter != containing_elelments.end();
                         ++index_iter)
                    {

                        assert(*index_iter<num_vertex_nodes);
                        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*index_iter);
                        solution_at_node += p_cell->GetCellData()->GetItem(mpPdeAndBcs->rGetDependentVariableName());
                    }
                    solution_at_node /= containing_elelments.size();
                }
            }
        }
        else if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation) != NULL ||
                 dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation) != NULL) // Intel compiler wants the "!= NULL"

        {
            // Simple 1-1 correspondence between cells and nodes in the finite element mesh
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);
            solution_at_node = p_cell->GetCellData()->GetItem(mpPdeAndBcs->rGetDependentVariableName());
        }
        else
        {
            // To do this would need to iterate over cells instead of loop over nodes.
            assert(dynamic_cast<CaBasedCellPopulation<DIM>*>(&rCellPopulation) != NULL); // Intel compiler wants the "!= NULL"
            EXCEPTION("ParabolicGrowingDomainPde Modifier doesn't work with CaBasedCellPopulations yet");
        }

        PetscVecTools::SetElement(this->mSolution,node_index,solution_at_node);
    }
}

template<unsigned DIM>
void ParabolicGrowingDomainPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractGrowingDomainPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ParabolicGrowingDomainPdeModifier<1>;
template class ParabolicGrowingDomainPdeModifier<2>;
template class ParabolicGrowingDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ParabolicGrowingDomainPdeModifier)

