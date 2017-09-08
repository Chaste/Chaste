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

#ifndef ELLIPTICGROWINGDOMAINPDESYSTEMMODIFIER_HPP_
#define ELLIPTICGROWINGDOMAINPDESYSTEMMODIFIER_HPP_

#include "AbstractGrowingDomainPdeSystemModifier.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "CellBasedEllipticPdeSystemSolver.hpp"
#include "AveragedSourceEllipticPde.hpp"

/**
 * A modifier class in which a linear elliptic PDE coupled to a cell-based simulation
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
 * CellwiseSourceEllipticPde and UniformSourceEllipticPde.
 */
template<unsigned DIM, unsigned PROBLEM_DIM=1>
class EllipticGrowingDomainPdeSystemModifier : public AbstractGrowingDomainPdeSystemModifier<DIM, PROBLEM_DIM>
{
    friend class TestEllipticGrowingDomainPdeSystemModifier;

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
    EllipticGrowingDomainPdeSystemModifier(
        boost::shared_ptr<AbstractLinearPdeSystem<DIM, DIM, PROBLEM_DIM> > pPdeSystem=boost::shared_ptr<AbstractLinearPdeSystem<DIM, DIM, PROBLEM_DIM> >(),
        std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions=std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > >(),
        bool isNeumannBoundaryCondition=true,
        Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~EllipticGrowingDomainPdeSystemModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to construct the boundary conditions container for the PDE.
     *
     * @return the full boundary conditions container
     */
    virtual std::shared_ptr<BoundaryConditionsContainer<DIM, DIM, PROBLEM_DIM> > ConstructBoundaryConditionsContainer();

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
EllipticGrowingDomainPdeSystemModifier<DIM, PROBLEM_DIM>::EllipticGrowingDomainPdeSystemModifier(
    boost::shared_ptr<AbstractLinearPdeSystem<DIM, DIM, PROBLEM_DIM> > pPdeSystem,
    std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions,
    bool isNeumannBoundaryCondition,
    Vec solution)
    : AbstractGrowingDomainPdeSystemModifier<DIM, PROBLEM_DIM>(pPdeSystem,
                                                              pBoundaryConditions,
                                                              isNeumannBoundaryCondition,
                                                              solution)
{
}

template <unsigned DIM, unsigned PROBLEM_DIM>
EllipticGrowingDomainPdeSystemModifier<DIM, PROBLEM_DIM>::~EllipticGrowingDomainPdeSystemModifier()
{
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void EllipticGrowingDomainPdeSystemModifier<DIM, PROBLEM_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    this->GenerateFeMesh(rCellPopulation);

    // If the solution at the previous timestep exists...
    PetscInt previous_solution_size = 0;
    if (this->mSolution)
    {
        VecGetSize(this->mSolution, &previous_solution_size);
    }

    // ...then record whether it is the correct size...
    bool is_previous_solution_size_correct = (previous_solution_size == (int)this->mpFeMesh->GetNumNodes());

    // ...and if it is, store it as an initial guess for the PDE solver
    Vec initial_guess;
    if (is_previous_solution_size_correct)
    {
        // This Vec is copied by the solver's Solve() method, so must be deleted here too
        VecDuplicate(this->mSolution, &initial_guess);
        VecCopy(this->mSolution, initial_guess);
        PetscTools::Destroy(this->mSolution);
    }

    // Add the BCs to the BCs container
    std::shared_ptr<BoundaryConditionsContainer<DIM, DIM, PROBLEM_DIM> > p_bcc = this->ConstructBoundaryConditionsContainer();

    // Use CellBasedEllipticPdeSystemSolver as cell wise PDE
    CellBasedEllipticPdeSystemSolver<DIM, PROBLEM_DIM> solver(this->mpFeMesh,
                                                              boost::static_pointer_cast<AbstractLinearEllipticPdeSystem<DIM,DIM,PROBLEM_DIM> >(this->GetPdeSystem()).get(),
                                                              p_bcc.get());

    // If we have an initial guess, use this when solving the system...
    if (is_previous_solution_size_correct)
    {
        this->mSolution = solver.Solve(initial_guess);
        PetscTools::Destroy(initial_guess);
    }
    else // ...otherwise do not supply one
    {
        // The solver creates a Vec, so we have to keep a handle on the old one to destroy it
        Vec old_solution_copy = this->mSolution;

        this->mSolution = solver.Solve();

        // On the first go round the vector has yet to be initialised, so we don't destroy it
        if (old_solution_copy != nullptr)
        {
            PetscTools::Destroy(old_solution_copy);
        }
    }

    this->UpdateCellData(rCellPopulation);
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void EllipticGrowingDomainPdeSystemModifier<DIM, PROBLEM_DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
    if (boost::dynamic_pointer_cast<AveragedSourceEllipticPde<DIM> >(this->GetPdeSystem()))
    {
        EXCEPTION("EllipticGrowingDomainPdeSystemModifier cannot be used with an AveragedSourceEllipticPde. Use an EllipticBoxDomainPdeSystemModifier instead.");
    }

    AbstractGrowingDomainPdeSystemModifier<DIM, PROBLEM_DIM>::SetupSolve(rCellPopulation, outputDirectory);

    // Call these methods to solve the PDE on the initial step and output the results
    UpdateAtEndOfTimeStep(rCellPopulation);
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template <unsigned DIM, unsigned PROBLEM_DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM, DIM, PROBLEM_DIM> > EllipticGrowingDomainPdeSystemModifier<DIM, PROBLEM_DIM>::ConstructBoundaryConditionsContainer()
{
    std::shared_ptr<BoundaryConditionsContainer<DIM, DIM, PROBLEM_DIM> > p_bcc(new BoundaryConditionsContainer<DIM, DIM, PROBLEM_DIM>(false));

    // To be well-defined, elliptic PDE problems on growing domains require Dirichlet boundary conditions
    assert(!(this->IsNeumannBoundaryCondition()));
    for (typename TetrahedralMesh<DIM, DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
         node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
         ++node_iter)
    {
        // Loop over PDEs
        for (unsigned pde_index=0; pde_index<PROBLEM_DIM; pde_index++)
        {
            p_bcc->AddDirichletBoundaryCondition(*node_iter, this->mpBoundaryConditions[pde_index].get(), pde_index);
        }
    }

    return p_bcc;
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void EllipticGrowingDomainPdeSystemModifier<DIM, PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractGrowingDomainPdeSystemModifier<DIM, PROBLEM_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

#endif /*ELLIPTICGROWINGDOMAINPDESYSTEMMODIFIER_HPP_*/
