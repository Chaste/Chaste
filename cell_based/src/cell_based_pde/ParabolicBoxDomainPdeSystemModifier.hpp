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

#ifndef PARABOLICBOXDOMAINPDESYSTEMMODIFIER_HPP_
#define PARABOLICBOXDOMAINPDESYSTEMMODIFIER_HPP_

#include "AbstractBoxDomainPdeSystemModifier.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "LinearParabolicPdeSystemSolver.hpp"

/**
 * A modifier class in which a linear parabolic PDE coupled to a cell-based simulation
 * is solved on a coarse domain.
 *
 * The finite element mesh used to solve the PDE numerically is a fixed tessellation of
 * a cuboid (box), which must be supplied to the constructor. The value of the dependent
 * variable is interpolated between coarse mesh nodes to obtain a value at each cell,
 * which is stored and updated in a CellData item.
 *
 * At each time step the boundary condition supplied to the constructor may be imposed
 * either on the boundary of the box domain, or on the boundary of the cell population
 * (which is assumed to lie within the box domain). This choice can be made using the
 * AbstractBoxDomainPdeSystemModifier method SetBcsOnBoxBoundary(), which is inherited by this
 * class.
 *
 * Examples of PDEs in the source folder that can be solved using this class are
 * AveragedSourceParabolicPde and UniformSourceParabolicPde.
 */
template<unsigned DIM, unsigned PROBLEM_DIM=1>
class ParabolicBoxDomainPdeSystemModifier : public AbstractBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>
{
    friend class TestParabolicBoxDomainPdeSystemModifier;

public:

    /**
     * Constructor.
     *
     * @param pPdeSystem A shared pointer to a linear PDE system object (defaults to NULL)
     * @param pBoundaryConditions A vector of shared pointers to abstract boundary conditions
     *     (defaults to NULL, corresponding to a constant boundary condition with value zero)
     * @param isNeumannBoundaryCondition Whether the boundary condition is Neumann (defaults to true)
     * @param pMeshCuboid A shared pointer to a ChasteCuboid specifying the outer boundary for the FE mesh (defaults to NULL)
     * @param stepSize step size to be used in the FE mesh (defaults to 1.0, i.e. the default cell size)
     * @param solution solution vector (defaults to NULL)
     */
    ParabolicBoxDomainPdeSystemModifier(
        boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > pPdeSystem=boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> >(),
        std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions=std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > >(),
        bool isNeumannBoundaryCondition=true,
        boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid=boost::shared_ptr<ChasteCuboid<DIM> >(),
        double stepSize=1.0,
        Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~ParabolicBoxDomainPdeSystemModifier();

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
     * @param rCellPopulation reference to the cell population
     *
     * @return the full boundary conditions container
     */
    virtual std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Helper method to initialise the PDE solution using the CellData.
     *
     * Here we assume a homogeneous initial consition.
     *
     * @param rCellPopulation reference to the cell population
     */
    void SetupInitialSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

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
ParabolicBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>::ParabolicBoxDomainPdeSystemModifier(
    boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > pPdeSystem,
    std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions,
    bool isNeumannBoundaryCondition,
    boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid,
    double stepSize,
    Vec solution)
    : AbstractBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>(pPdeSystem,
                                              pBoundaryConditions,
                                              isNeumannBoundaryCondition,
                                              pMeshCuboid,
                                              stepSize,
                                              solution)
{
}

template <unsigned DIM, unsigned PROBLEM_DIM>
ParabolicBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>::~ParabolicBoxDomainPdeSystemModifier()
{
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void ParabolicBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > p_bcc = ConstructBoundaryConditionsContainer(rCellPopulation);

    this->UpdateCellPdeElementMap(rCellPopulation);

    // When using a PDE mesh which doesn't coincide with the cells, we must set up the source terms before solving the PDE.
    // Pass in already updated CellPdeElementMap to speed up finding cells.
    this->SetUpSourceTermsForAveragedSourcePde(this->mpFeMesh, &this->mCellPdeElementMap);

    // Use LinearParabolicPdeSystemSolver as averaged Source PDE
    LinearParabolicPdeSystemSolver<DIM,DIM,PROBLEM_DIM> solver(this->mpFeMesh,
                                                            boost::static_pointer_cast<AbstractLinearParabolicPdeSystem<DIM,DIM,PROBLEM_DIM> >(this->GetPdeSystem()).get(),
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
void ParabolicBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>::SetupSolve(rCellPopulation,outputDirectory);

    // Copy the cell data to mSolution (this is the initial condition)
    SetupInitialSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template <unsigned DIM, unsigned PROBLEM_DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > ParabolicBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>::ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,PROBLEM_DIM>(false));

    if (!this->mSetBcsOnBoxBoundary)
    {
        EXCEPTION("Boundary conditions cannot yet be set on the cell population boundary for a ParabolicBoxDomainPdeSystemModifier");
    }
    else // Apply BC at boundary nodes of box domain FE mesh
    {
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
    }

    return p_bcc;
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void ParabolicBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>::SetupInitialSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Specify homogeneous initial conditions based upon the values stored in CellData.
    // Note need all the CellDataValues to be the same.

    double initial_condition = rCellPopulation.Begin()->GetCellData()->GetItem(this->mDependentVariableNames[0]);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        double initial_condition_at_cell = cell_iter->GetCellData()->GetItem(this->mDependentVariableNames[0]);
        UNUSED_OPT(initial_condition_at_cell);
        assert(fabs(initial_condition_at_cell - initial_condition)<1e-12);
    }

    // Initialise mSolution
    this->mSolution = PetscTools::CreateAndSetVec(this->mpFeMesh->GetNumNodes(), initial_condition);
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void ParabolicBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractBoxDomainPdeSystemModifier<DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

#endif /*PARABOLICBOXDOMAINPDESYSTEMMODIFIER_HPP_*/
