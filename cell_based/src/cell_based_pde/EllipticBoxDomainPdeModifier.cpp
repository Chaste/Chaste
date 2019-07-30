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

#include "EllipticBoxDomainPdeModifier.hpp"
#include "SimpleLinearEllipticSolver.hpp"

template<unsigned DIM>
EllipticBoxDomainPdeModifier<DIM>::EllipticBoxDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde,
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
                                        solution)
{
}

template<unsigned DIM>
EllipticBoxDomainPdeModifier<DIM>::~EllipticBoxDomainPdeModifier()
{
}

template<unsigned DIM>
void EllipticBoxDomainPdeModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Set up boundary conditions
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer(rCellPopulation);

    this->UpdateCellPdeElementMap(rCellPopulation);

    // When using a PDE mesh which doesn't coincide with the cells, we must set up the source terms before solving the PDE.
    // Pass in already updated CellPdeElementMap to speed up finding cells.
    this->SetUpSourceTermsForAveragedSourcePde(this->mpFeMesh, &this->mCellPdeElementMap);

    // Use SimpleLinearEllipticSolver as Averaged Source PDE
    ///\todo allow other PDE classes to be used with this modifier
    SimpleLinearEllipticSolver<DIM,DIM> solver(this->mpFeMesh,
                                               boost::static_pointer_cast<AbstractLinearEllipticPde<DIM,DIM> >(this->GetPde()).get(),
                                               p_bcc.get());

    ///\todo Use solution at previous time step as an initial guess for Solve()
    Vec old_solution_copy = this->mSolution;
    this->mSolution = solver.Solve();
    if (old_solution_copy != nullptr)
    {
        PetscTools::Destroy(old_solution_copy);
    }

    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void EllipticBoxDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractBoxDomainPdeModifier<DIM>::SetupSolve(rCellPopulation,outputDirectory);

    // Call these  methods to solve the PDE on the initial step and output the results
    UpdateAtEndOfTimeStep(rCellPopulation);
    this->UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > EllipticBoxDomainPdeModifier<DIM>::ConstructBoundaryConditionsContainer(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    std::shared_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    // To be well-defined, elliptic PDE problems on box domains require at least some Dirichlet boundary conditions
    ///\todo Replace this assertion with an exception in the constructor
    assert(!(this->IsNeumannBoundaryCondition()));

    if (!this->mSetBcsOnBoxBoundary)
    {
        // Get the set of coarse element indices that contain cells
        std::set<unsigned> coarse_element_indices_in_map;
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            coarse_element_indices_in_map.insert(this->mCellPdeElementMap[*cell_iter]);
        }

        // Find the node indices associated with elements whose indices are NOT in the set coarse_element_indices_in_map
        std::set<unsigned> coarse_mesh_boundary_node_indices;
        for (unsigned i=0; i<this->mpFeMesh->GetNumElements(); i++)
        {
            if (coarse_element_indices_in_map.find(i) == coarse_element_indices_in_map.end())
            {
                Element<DIM,DIM>* p_element = this->mpFeMesh->GetElement(i);
                for (unsigned j=0; j<DIM+1; j++)
                {
                    unsigned node_index = p_element->GetNodeGlobalIndex(j);
                    coarse_mesh_boundary_node_indices.insert(node_index);
                }
            }
        }

        // Apply boundary condition to the nodes in the set coarse_mesh_boundary_node_indices
        for (std::set<unsigned>::iterator iter = coarse_mesh_boundary_node_indices.begin();
             iter != coarse_mesh_boundary_node_indices.end();
             ++iter)
        {
            p_bcc->AddDirichletBoundaryCondition(this->mpFeMesh->GetNode(*iter), this->mpBoundaryCondition.get(), 0, false);
        }
    }
    else // Apply BC at boundary nodes of box domain FE mesh
    {
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

