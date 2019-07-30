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

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "CellBasedParabolicPdeSolver.hpp"
#include "AveragedSourceParabolicPde.hpp"

template<unsigned DIM>
ParabolicGrowingDomainPdeModifier<DIM>::ParabolicGrowingDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde,
                                                                          boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition,
                                                                          bool isNeumannBoundaryCondition,
                                                                          Vec solution)
    : AbstractGrowingDomainPdeModifier<DIM>(pPde,
                                            pBoundaryCondition,
                                            isNeumannBoundaryCondition,
                                            solution)
{
}

template<unsigned DIM>
ParabolicGrowingDomainPdeModifier<DIM>::~ParabolicGrowingDomainPdeModifier()
{
}

template<unsigned DIM>
void ParabolicGrowingDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    this->GenerateFeMesh(rCellPopulation);

    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer();

    // Construct the solution vector from cell data (takes care of cells dividing);
    UpdateSolutionVector(rCellPopulation);

    // Use CellBasedParabolicPdeSolver as cell wise PDE
    CellBasedParabolicPdeSolver<DIM> solver(this->mpFeMesh,
                                            boost::static_pointer_cast<AbstractLinearParabolicPde<DIM,DIM> >(this->mpPde).get(),
                                            p_bcc.get());

    ///\todo Investigate more than one PDE time step per spatial step
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();
    double dt = p_simulation_time->GetTimeStep();
    solver.SetTimes(current_time,current_time + dt);
    solver.SetTimeStep(dt);

    // Use previous solution as the initial condition
    Vec previous_solution = this->mSolution;
    solver.SetInitialCondition(previous_solution);

    // Note that the linear solver creates a vector, so we have to keep a handle on the old one
    // in order to destroy it
    this->mSolution = solver.Solve();
    PetscTools::Destroy(previous_solution);
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ParabolicGrowingDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractGrowingDomainPdeModifier<DIM>::SetupSolve(rCellPopulation, outputDirectory);

    if (boost::dynamic_pointer_cast<AveragedSourceParabolicPde<DIM> >(this->mpPde))
    {
        EXCEPTION("ParabolicGrowingDomainPdeModifier cannot be used with an AveragedSourceParabolicPde. Use a ParabolicBoxDomainPdeModifier instead.");
    }

    // Setup a finite element mesh on which to save the initial condition
    this->GenerateFeMesh(rCellPopulation);

    // Copy the cell data to mSolution (this is the initial condition)
    UpdateSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ParabolicGrowingDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer()
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    if (this->IsNeumannBoundaryCondition())
    {
        // Impose any Neumann boundary conditions
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
             elem_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
             ++elem_iter)
        {
            p_bcc->AddNeumannBoundaryCondition(*elem_iter, this->mpBoundaryCondition.get());
        }
    }
    else
    {
        // Impose any Dirichlet boundary conditions
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
             node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
             ++node_iter)
        {
            p_bcc->AddDirichletBoundaryCondition(*node_iter, this->mpBoundaryCondition.get());
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

    std::string& variable_name = this->mDependentVariableName;

    for (typename TetrahedralMesh<DIM,DIM>::NodeIterator node_iter = this->mpFeMesh->GetNodeIteratorBegin();
         node_iter != this->mpFeMesh->GetNodeIteratorEnd();
         ++node_iter)
    {
        // Loop over nodes of the finite element mesh and get appropriate solution values from CellData
        for (typename TetrahedralMesh<DIM,DIM>::NodeIterator node_iter = this->mpFeMesh->GetNodeIteratorBegin();
             node_iter != this->mpFeMesh->GetNodeIteratorEnd();
             ++node_iter)
        {
            unsigned node_index = node_iter->GetIndex();
            bool dirichlet_bc_applies = (node_iter->IsBoundaryNode()) && (!(this->IsNeumannBoundaryCondition()));
            double boundary_value = this->GetBoundaryCondition()->GetValue(node_iter->rGetLocation());

            double solution_at_node = rCellPopulation.GetCellDataItemAtPdeNode(node_index, variable_name, dirichlet_bc_applies, boundary_value);

            PetscVecTools::SetElement(this->mSolution, node_index, solution_at_node);
        }
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

