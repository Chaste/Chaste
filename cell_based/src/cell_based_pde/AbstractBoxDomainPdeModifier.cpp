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

#include "AbstractBoxDomainPdeModifier.hpp"
#include "ReplicatableVector.hpp"
#include "LinearBasisFunction.hpp"

template<unsigned DIM>
AbstractBoxDomainPdeModifier<DIM>::AbstractBoxDomainPdeModifier(boost::shared_ptr<AbstractLinearPde<DIM,DIM> > pPde,
                                                                boost::shared_ptr<AbstractBoundaryCondition<DIM> > pBoundaryCondition,
                                                                bool isNeumannBoundaryCondition,
                                                                boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid,
                                                                double stepSize,
                                                                Vec solution)
    : AbstractPdeModifier<DIM>(pPde,
                               pBoundaryCondition,
                               isNeumannBoundaryCondition,
                               solution),
      mpMeshCuboid(pMeshCuboid),
      mStepSize(stepSize),
      mSetBcsOnBoxBoundary(true)
{
    if (pMeshCuboid)
    {
        // We only need to generate mpFeMesh once, as it does not vary with time
        this->GenerateFeMesh(mpMeshCuboid, mStepSize);
        this->mDeleteFeMesh = true;
    }
}

template<unsigned DIM>
AbstractBoxDomainPdeModifier<DIM>::~AbstractBoxDomainPdeModifier()
{
}

template<unsigned DIM>
double AbstractBoxDomainPdeModifier<DIM>::GetStepSize()
{
     return mStepSize;
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::SetBcsOnBoxBoundary(bool setBcsOnBoxBoundary)
{
    mSetBcsOnBoxBoundary = setBcsOnBoxBoundary;
}

template<unsigned DIM>
bool AbstractBoxDomainPdeModifier<DIM>::AreBcsSetOnBoxBoundary()
{
    return mSetBcsOnBoxBoundary;
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    AbstractPdeModifier<DIM>::SetupSolve(rCellPopulation, outputDirectory);

    InitialiseCellPdeElementMap(rCellPopulation);
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::GenerateFeMesh(boost::shared_ptr<ChasteCuboid<DIM> > pMeshCuboid, double stepSize)
{
    // Create a regular coarse tetrahedral mesh
    this->mpFeMesh = new TetrahedralMesh<DIM,DIM>();
    switch (DIM)
    {
        case 1:
            this->mpFeMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0));
            break;
        case 2:
            this->mpFeMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0), pMeshCuboid->GetWidth(1));
            break;
        case 3:
            this->mpFeMesh->ConstructRegularSlabMesh(stepSize, pMeshCuboid->GetWidth(0), pMeshCuboid->GetWidth(1), pMeshCuboid->GetWidth(2));
            break;
        default:
            NEVER_REACHED;
    }

    // Get centroid of meshCuboid
    ChastePoint<DIM> upper = pMeshCuboid->rGetUpperCorner();
    ChastePoint<DIM> lower = pMeshCuboid->rGetLowerCorner();
    c_vector<double,DIM> centre_of_cuboid = 0.5*(upper.rGetLocation() + lower.rGetLocation());

    // Find the centre of the PDE mesh
    c_vector<double,DIM> centre_of_coarse_mesh = zero_vector<double>(DIM);
    for (unsigned i=0; i<this->mpFeMesh->GetNumNodes(); i++)
    {
        centre_of_coarse_mesh += this->mpFeMesh->GetNode(i)->rGetLocation();
    }
    centre_of_coarse_mesh /= this->mpFeMesh->GetNumNodes();

    // Now move the mesh to the correct location
    this->mpFeMesh->Translate(centre_of_cuboid - centre_of_coarse_mesh);
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Store the PDE solution in an accessible form
    ReplicatableVector solution_repl(this->mSolution);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // The cells are not nodes of the mesh, so we must interpolate
        double solution_at_cell = 0.0;

        // Find the element in the FE mesh that contains this cell. CellElementMap has been updated so use this.
        unsigned elem_index = mCellPdeElementMap[*cell_iter];
        Element<DIM,DIM>* p_element = this->mpFeMesh->GetElement(elem_index);

        const ChastePoint<DIM>& node_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        c_vector<double,DIM+1> weights = p_element->CalculateInterpolationWeights(node_location);

        for (unsigned i=0; i<DIM+1; i++)
        {
            double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(i)];
            solution_at_cell += nodal_value * weights(i);
        }

        cell_iter->GetCellData()->SetItem(this->mDependentVariableName, solution_at_cell);

        if (this->mOutputGradient)
        {
            // Now calculate the gradient of the solution and store this in CellVecData
            c_vector<double, DIM> solution_gradient = zero_vector<double>(DIM);

            // Calculate the basis functions at any point (e.g. zero) in the element
            c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
            double jacobian_det;
            this->mpFeMesh->GetInverseJacobianForElement(elem_index, jacobian, jacobian_det, inverse_jacobian);
            const ChastePoint<DIM> zero_point;
            c_matrix<double, DIM, DIM+1> grad_phi;
            LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(zero_point, inverse_jacobian, grad_phi);

            for (unsigned node_index=0; node_index<DIM+1; node_index++)
            {
                double nodal_value = solution_repl[p_element->GetNodeGlobalIndex(node_index)];

                for (unsigned j=0; j<DIM; j++)
                {
                    solution_gradient(j) += nodal_value* grad_phi(j, node_index);
                }
            }

            switch (DIM)
            {
                case 1:
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_x", solution_gradient(0));
                    break;
                case 2:
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_x", solution_gradient(0));
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_y", solution_gradient(1));
                    break;
                case 3:
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_x", solution_gradient(0));
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_y", solution_gradient(1));
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableName+"_grad_z", solution_gradient(2));
                    break;
                default:
                    NEVER_REACHED;
            }
        }
    }
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::InitialiseCellPdeElementMap(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    mCellPdeElementMap.clear();

    // Find the element of mpFeMesh that contains each cell and populate mCellPdeElementMap
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = this->mpFeMesh->GetContainingElementIndex(r_position_of_cell);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::UpdateCellPdeElementMap(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Find the element of mpCoarsePdeMesh that contains each cell and populate mCellPdeElementMap
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        const ChastePoint<DIM>& r_position_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        unsigned elem_index = this->mpFeMesh->GetContainingElementIndexWithInitialGuess(r_position_of_cell, mCellPdeElementMap[*cell_iter]);
        mCellPdeElementMap[*cell_iter] = elem_index;
    }
}

template<unsigned DIM>
void AbstractBoxDomainPdeModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractPdeModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class AbstractBoxDomainPdeModifier<1>;
template class AbstractBoxDomainPdeModifier<2>;
template class AbstractBoxDomainPdeModifier<3>;
