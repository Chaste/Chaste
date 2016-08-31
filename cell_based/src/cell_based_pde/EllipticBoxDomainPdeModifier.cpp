/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "EllipticBoxDomainPdeModifier.hpp"
#include "PdeAndBoundaryConditions.hpp"
#include "SimpleLinearEllipticSolver.hpp"

template<unsigned DIM>
EllipticBoxDomainPdeModifier<DIM>::EllipticBoxDomainPdeModifier(boost::shared_ptr<PdeAndBoundaryConditions<DIM> > pPdeAndBcs,
                                                                ChasteCuboid<DIM>* pMeshCuboid,
                                                                double stepSize)
    : AbstractBoxDomainPdeModifier<DIM>(pPdeAndBcs, pMeshCuboid, stepSize)
{
    assert(DIM == 2);
}

template<unsigned DIM>
EllipticBoxDomainPdeModifier<DIM>::~EllipticBoxDomainPdeModifier()
{
    // Destroy the most recent solution vector
    if (this->mSolution != NULL)
    {
        PetscTools::Destroy(this->mSolution);
    }
}

template<unsigned DIM>
void EllipticBoxDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Set up boundary conditions
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer();

    this->UpdateCellPdeElementMap(rCellPopulation);

    // When using a PDE mesh which doesn't coincide with the cells, we must set up the source terms before solving the PDE.
    // Pass in already updated CellPdeElementMap to speed up finding cells.
    this->mpPdeAndBcs->SetUpSourceTermsForAveragedSourcePde(this->mpFeMesh, &this->mCellPdeElementMap);

    // Use SimpleLinearEllipticSolver as Averaged Source PDE
    SimpleLinearEllipticSolver<DIM,DIM> solver(this->mpFeMesh,
                                               static_cast<AbstractLinearEllipticPde<DIM,DIM>*>(this->mpPdeAndBcs->GetPde()),
                                               p_bcc.get());

    ///\todo Use initial guess when solving the system (#2687)
    Vec old_solution_copy = this->mSolution;
    this->mSolution = solver.Solve();

    // Note that the linear solver creates a vector, so we have to keep a handle on the old one
    // in order to destroy it.
    ///\todo #2687 This will change when initial guess is used.
    /// On the first go round the vector has yet to be initialised, so we don't destroy it.
    if (old_solution_copy != NULL)
    {
        PetscTools::Destroy(old_solution_copy);
    }
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void EllipticBoxDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractBoxDomainPdeModifier<DIM>::SetupSolve(rCellPopulation,outputDirectory);

    // Call these  methods to solve the PDE on the initial step and Output the results.
    UpdateAtEndOfTimeStep(rCellPopulation);
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > EllipticBoxDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer()
{
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    // To be well-defined, elliptic PDE problems on box domains require at least some Dirichlet boundary conditions
    ///\todo Replace this assertion with an exception in the PdeAndBoundaryConditions constructor
    assert(!(this->mpPdeAndBcs->IsNeumannBoundaryCondition()));

    for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = this->mpFeMesh->GetBoundaryNodeIteratorBegin();
         node_iter != this->mpFeMesh->GetBoundaryNodeIteratorEnd();
         ++node_iter)
    {
        p_bcc->AddDirichletBoundaryCondition(*node_iter, this->mpPdeAndBcs->GetBoundaryCondition());
    }

    return p_bcc;
}

template<unsigned DIM>
void EllipticBoxDomainPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractBoxDomainPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class EllipticBoxDomainPdeModifier<1>;
template class EllipticBoxDomainPdeModifier<2>;
template class EllipticBoxDomainPdeModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EllipticBoxDomainPdeModifier)

