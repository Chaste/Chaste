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

#include "EllipticGrowingDomainPdeModifier.hpp"
#include "CellBasedEllipticPdeSolver.hpp"
#include "AveragedSourceEllipticPde.hpp"

template<unsigned DIM>
EllipticGrowingDomainPdeModifier<DIM>::EllipticGrowingDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde,
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
EllipticGrowingDomainPdeModifier<DIM>::~EllipticGrowingDomainPdeModifier()
{
}

template<unsigned DIM>
void EllipticGrowingDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = this->ConstructBoundaryConditionsContainer();

    // Use CellBasedEllipticPdeSolver as cell wise PDE
    CellBasedEllipticPdeSolver<DIM> solver(this->mpFeMesh,
                                           boost::static_pointer_cast<AbstractLinearEllipticPde<DIM,DIM> >(this->GetPde()).get(),
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

template<unsigned DIM>
void EllipticGrowingDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    if (boost::dynamic_pointer_cast<AveragedSourceEllipticPde<DIM> >(this->GetPde()))
    {
        EXCEPTION("EllipticGrowingDomainPdeModifier cannot be used with an AveragedSourceEllipticPde. Use an EllipticBoxDomainPdeModifier instead.");
    }

    AbstractGrowingDomainPdeModifier<DIM>::SetupSolve(rCellPopulation, outputDirectory);

    // Call these methods to solve the PDE on the initial step and output the results
    UpdateAtEndOfTimeStep(rCellPopulation);
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > EllipticGrowingDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer()
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    // To be well-defined, elliptic PDE problems on growing domains require Dirichlet boundary conditions
    assert(!(this->IsNeumannBoundaryCondition()));
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
         node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
         ++node_iter)
    {
        p_bcc->AddDirichletBoundaryCondition(*node_iter, this->mpBoundaryCondition.get());
    }

    return p_bcc;
}

template<unsigned DIM>
void EllipticGrowingDomainPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractGrowingDomainPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class EllipticGrowingDomainPdeModifier<1>;
template class EllipticGrowingDomainPdeModifier<2>;
template class EllipticGrowingDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EllipticGrowingDomainPdeModifier)
