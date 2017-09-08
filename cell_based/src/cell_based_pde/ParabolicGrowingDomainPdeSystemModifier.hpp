/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef PARABOLICGROWINGDOMAINPDESYSTEMMODIFIER_HPP_
#define PARABOLICGROWINGDOMAINPDESYSTEMMODIFIER_HPP_

#include "AbstractGrowingDomainPdeSystemModifier.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "CellBasedParabolicPdeSystemSolver.hpp"
#include "AveragedSourceParabolicPde.hpp"

/**
 * A modifier class in which a linear parabolic PDE coupled to a cell-based simulation
 * is solved on a growing domain. The value of the dependent variable at each cell is
 * stored and updated in a CellData item.
 *
 * At each time step, the finite element mesh used to solve the PDE numerically
 * is defined by the spatial domain associated with the cell population. The
 * precise definition of this domain is implemented in the method
 * GetTetrahedralMeshForPdeModifier(), which is overridden for each cell population
 * class and is used in the AbstractGrowingDomainPdeSystemModifier method GenerateFeMesh()
 * that is inherited by this class.
 *
 * Examples of PDEs in the source folder that can be solved using this class are
 * CellwiseSourceParabolicPde and UniformSourceParabolicPde.
 */
template<unsigned DIM, unsigned PROBLEM_DIM=1>
class ParabolicGrowingDomainPdeSystemModifier : public AbstractGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>
{
    friend class TestParabolicGrowingDomainPdeSystemModifier;

public:

    /**
     * Constructor.
     *
     * @param pPdeSystem A shared pointer to a linear PDE system object (defaults to NULL)
     * @param pBoundaryConditions A vector of shared pointers to abstract boundary conditions
     *     (defaults to NULL, corresponding to a constant boundary condition with value zero)
     * @param isNeumannBoundaryCondition Whether the boundary condition is Neumann (defaults to true)
     * @param solution solution vector (defaults to NULL)
     */
    ParabolicGrowingDomainPdeSystemModifier(
        boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > pPdeSystem=boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> >(),
        std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions=std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > >(),
        bool isNeumannBoundaryCondition=true,
        Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~ParabolicGrowingDomainPdeSystemModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to construct the boundary conditions container for the PDE.
     *
     * @return the full boundary conditions container
     */
    virtual std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > ConstructBoundaryConditionsContainer();

    /**
     * Helper method to copy the CellData to the PDE solution.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

/*
 * As this class is templated over PROBLEM_DIM, we put the implementation
 * in the header file to avoid explicit instantiation.
 */

template <unsigned DIM, unsigned PROBLEM_DIM>
ParabolicGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::ParabolicGrowingDomainPdeSystemModifier(
    boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > pPdeSystem,
    std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions,
    bool isNeumannBoundaryCondition,
    Vec solution)
    : AbstractGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>(pPdeSystem,
                                                              pBoundaryConditions,
                                                              isNeumannBoundaryCondition,
                                                              solution)
{
}

template <unsigned DIM, unsigned PROBLEM_DIM>
ParabolicGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::~ParabolicGrowingDomainPdeSystemModifier()
{
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void ParabolicGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    this->GenerateFeMesh(rCellPopulation);

    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > p_bcc = ConstructBoundaryConditionsContainer();

    // Construct the solution vector from cell data (takes care of cells dividing);
    UpdateSolutionVector(rCellPopulation);

    // Use CellBasedParabolicPdeSystemSolver as cell wise PDE
    CellBasedParabolicPdeSystemSolver<DIM, PROBLEM_DIM> solver(this->mpFeMesh,
                                                               boost::static_pointer_cast<AbstractLinearParabolicPdeSystem<DIM,DIM,PROBLEM_DIM> >(this->mpPdeSystem).get(),
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

template <unsigned DIM, unsigned PROBLEM_DIM>
void ParabolicGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::SetupSolve(rCellPopulation, outputDirectory);

    if (boost::dynamic_pointer_cast<AveragedSourceParabolicPde<DIM> >(this->mpPdeSystem))
    {
        EXCEPTION("ParabolicGrowingDomainPdeSystemModifier cannot be used with an AveragedSourceParabolicPde. Use a ParabolicBoxDomainPdeSystemModifier instead.");
    }

    // Setup a finite element mesh on which to save the initial condition
    this->GenerateFeMesh(rCellPopulation);

    // Copy the cell data to mSolution (this is the initial condition)
    UpdateSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template <unsigned DIM, unsigned PROBLEM_DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > ParabolicGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::ConstructBoundaryConditionsContainer()
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM>(false));

    if (this->IsNeumannBoundaryCondition())
    {
        // Impose any Neumann boundary conditions
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = this->mpFeMesh->GetBoundaryElementIteratorBegin();
             elem_iter != this->mpFeMesh->GetBoundaryElementIteratorEnd();
             ++elem_iter)
        {
            // Loop over PDEs
            for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
            {
                p_bcc->AddNeumannBoundaryCondition(*elem_iter, this->mpBoundaryConditions[pde_index].get(), pde_index);
            }
        }
    }
    else
    {
        // Impose any Dirichlet boundary conditions
        for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
             node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
             ++node_iter)
        {
            // Loop over PDEs
            for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
            {
                p_bcc->AddDirichletBoundaryCondition(*node_iter, this->mpBoundaryConditions[pde_index].get(), pde_index);
            }
        }
    }

    return p_bcc;
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void ParabolicGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::UpdateSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Clear (if it's not the first time) and resize the solution vector
    if (this->mSolution)
    {
        PetscTools::Destroy(this->mSolution);
    }
    this->mSolution = PetscTools::CreateAndSetVec(this->mpFeMesh->GetNumNodes(), 0.0);

    std::string& variable_name = this->mDependentVariableNames[0];

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

template <unsigned DIM, unsigned PROBLEM_DIM>
void ParabolicGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

#endif /*PARABOLICGROWINGDOMAINPDESYSTEMMODIFIER_HPP_*/
