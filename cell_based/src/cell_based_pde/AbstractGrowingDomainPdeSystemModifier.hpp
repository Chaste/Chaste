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

#ifndef ABSTRACTGROWINGDOMAINPDESYSTEMMODIFIER_HPP_
#define ABSTRACTGROWINGDOMAINPDESYSTEMMODIFIER_HPP_

#include "AbstractPdeSystemModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "ReplicatableVector.hpp"
#include "LinearBasisFunction.hpp"

/**
 * An abstract modifier class containing functionality common to EllipticGrowingDomainPdeSystemModifier
 * and ParabolicGrowingDomainPdeSystemModifier, which both solve a linear elliptic or parabolic PDE
 * coupled to a cell-based simulation on an evolving domain defined by the cell population.
 */
template<unsigned DIM, unsigned PROBLEM_DIM=1>
class AbstractGrowingDomainPdeSystemModifier : public AbstractPdeSystemModifier<DIM,PROBLEM_DIM>
{
    friend class TestEllipticGrowingDomainPdeSystemModifier;
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
    AbstractGrowingDomainPdeSystemModifier(
        boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > pPdeSystem=boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> >(),
        std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions=std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > >(),
        bool isNeumannBoundaryCondition=true,
        Vec solution=nullptr);

    /**
     * Destructor.
     */
    virtual ~AbstractGrowingDomainPdeSystemModifier();

    /**
     * Helper method to generate the mesh from the Cell population.
     *
     * @param rCellPopulation reference to the cell population
     */
    void GenerateFeMesh(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Helper method to copy the PDE solution to CellData
     *
     * Here there is a 1-1 correspondence between cells and nodes in the FE mesh
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);


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
AbstractGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::AbstractGrowingDomainPdeSystemModifier(
    boost::shared_ptr<AbstractLinearPdeSystem<DIM,DIM,PROBLEM_DIM> > pPdeSystem,
    std::vector<boost::shared_ptr<AbstractBoundaryCondition<DIM> > > pBoundaryConditions,
    bool isNeumannBoundaryCondition,
    Vec solution)
    : AbstractPdeSystemModifier<DIM,PROBLEM_DIM>(pPdeSystem,
                                                 pBoundaryConditions,
                                                 isNeumannBoundaryCondition,
                                                 solution)
{
}

template <unsigned DIM, unsigned PROBLEM_DIM>
AbstractGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::~AbstractGrowingDomainPdeSystemModifier()
{
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void AbstractGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::GenerateFeMesh(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (this->mDeleteFeMesh)
    {
        // If a mesh has been created on a previous time step then we need to tidy it up
        assert(this->mpFeMesh != nullptr);
        delete this->mpFeMesh;
    }
    else
    {
        ///\todo We should only set mDeleteFeMesh once, not every time step (#2687, #2863)
        // This placement assumes that if this->mDeleteFeMesh is false it is uninitialised and needs to
        // be checked. If true, it has been checked elsewhere.
        this->mDeleteFeMesh = (dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr);
    }

    // Get the finite element mesh via the cell population. Set to NULL first in case mesh generation fails.
    this->mpFeMesh = nullptr;
    this->mpFeMesh = rCellPopulation.GetTetrahedralMeshForPdeModifier();
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void AbstractGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Store the PDE solution in an accessible form
    ReplicatableVector solution_repl(this->mSolution);

    // Local cell index used by the CA simulation
    unsigned cell_index = 0;

    unsigned index_in_solution_repl = 0;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned tet_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        ///\todo Consider how to remove dynamic_casts here
        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr)
        {
            // Offset to relate elements in vertex mesh to nodes in tetrahedral mesh
            tet_node_index += rCellPopulation.GetNumNodes();
        }
        else if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr)
        {
            // Here local cell index corresponds to tet node
            tet_node_index = cell_index;
            cell_index++;
        }
        else if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation) != nullptr)
        {
            tet_node_index = index_in_solution_repl;
            index_in_solution_repl++;
        }

        double solution_at_node = solution_repl[tet_node_index];

        cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0], solution_at_node);

        if (this->mOutputGradient)
        {
            // Now calculate the gradient of the solution and store this in CellVecData
            c_vector<double, DIM> solution_gradient = zero_vector<double>(DIM);

            Node<DIM>* p_tet_node = this->mpFeMesh->GetNode(tet_node_index);

            // Get the containing elements and average the contribution from each one
            for (typename Node<DIM>::ContainingElementIterator element_iter = p_tet_node->ContainingElementsBegin();
                 element_iter != p_tet_node->ContainingElementsEnd();
                 ++element_iter)
            {
                // Calculate the basis functions at any point (eg zero) in the element
                c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
                double jacobian_det;
                this->mpFeMesh->GetInverseJacobianForElement(*element_iter, jacobian, jacobian_det, inverse_jacobian);
                const ChastePoint<DIM> zero_point;
                c_matrix<double, DIM, DIM+1> grad_phi;
                LinearBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(zero_point, inverse_jacobian, grad_phi);

                // Add the contribution from this element
                for (unsigned node_index=0; node_index<DIM+1; node_index++)
                {
                    double nodal_value = solution_repl[this->mpFeMesh->GetElement(*element_iter)->GetNodeGlobalIndex(node_index)];

                    for (unsigned j=0; j<DIM; j++)
                    {
                        solution_gradient(j) += nodal_value* grad_phi(j, node_index);
                    }
                }
            }

            // Divide by number of containing elements
            solution_gradient /= p_tet_node->GetNumContainingElements();

            switch (DIM)
            {
                case 1:
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_x", solution_gradient(0));
                    break;
                case 2:
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_x", solution_gradient(0));
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_y", solution_gradient(1));
                    break;
                case 3:
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_x", solution_gradient(0));
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_y", solution_gradient(1));
                    cell_iter->GetCellData()->SetItem(this->mDependentVariableNames[0]+"_grad_z", solution_gradient(2));
                    break;
                default:
                    NEVER_REACHED;
            }
        }
    }
}

template <unsigned DIM, unsigned PROBLEM_DIM>
void AbstractGrowingDomainPdeSystemModifier<DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractPdeSystemModifier<DIM,PROBLEM_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

#endif /*ABSTRACTGROWINGDOMAINPDESYSTEMMODIFIER_HPP_*/
