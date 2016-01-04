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

#include "CellBasedPdeHandlerOnCuboid.hpp"
#include "ConstBoundaryCondition.hpp"

template<unsigned DIM>
CellBasedPdeHandlerOnCuboid<DIM>::CellBasedPdeHandlerOnCuboid(AbstractCellPopulation<DIM>* pCellPopulation,
                                              bool deleteMemberPointersInDestructor)
    : CellBasedPdeHandler<DIM>(pCellPopulation, deleteMemberPointersInDestructor)
{
}

template<unsigned DIM>
CellBasedPdeHandlerOnCuboid<DIM>::~CellBasedPdeHandlerOnCuboid()
{
    // Delete all boundary conditions
    for (unsigned i=0; i<mConstBoundaryConditions.size(); i++)
    {
        delete mConstBoundaryConditions[i];
    }
}

template<unsigned DIM>
std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > CellBasedPdeHandlerOnCuboid<DIM>::ConstructBoundaryConditionsContainer(
        PdeAndBoundaryConditions<DIM>* pPdeAndBc,
        TetrahedralMesh<DIM,DIM>* pMesh)
{
    // Not using the inputs as there's only one BC
    assert(DIM==2);

    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

    // Use these two vectors to define what's happening on the top right bottom and left boundaries

    // Specify which sides are Neumann boundaries
    c_vector<bool,4> are_neumann_boundaries;
    are_neumann_boundaries[0] = false; // Top
    are_neumann_boundaries[1] = true;  // Right
    are_neumann_boundaries[2] = false; // Bottom
    are_neumann_boundaries[3] = true;  // Left

    // Specify the value of the boundary condition on each boundary
    c_vector<double,4> boundary_condition_values;
    boundary_condition_values[0] = 1.0; // Top
    boundary_condition_values[1] = 1.0; // Right
    boundary_condition_values[2] = 1.0; // Bottom
    boundary_condition_values[3] = 1.0; // Left

    // Delete all previous boundary conditions
    for (unsigned i=0; i<mConstBoundaryConditions.size(); i++)
    {
        delete mConstBoundaryConditions[i];
    }
    mConstBoundaryConditions.clear();

    mConstBoundaryConditions.push_back(new ConstBoundaryCondition<DIM>(boundary_condition_values[0]));
    mConstBoundaryConditions.push_back(new ConstBoundaryCondition<DIM>(boundary_condition_values[1]));
    mConstBoundaryConditions.push_back(new ConstBoundaryCondition<DIM>(boundary_condition_values[2]));
    mConstBoundaryConditions.push_back(new ConstBoundaryCondition<DIM>(boundary_condition_values[3]));

    ChasteCuboid<DIM> cuboid = pMesh->CalculateBoundingBox();

    double left = cuboid.rGetLowerCorner()[0];
    double bottom = cuboid.rGetLowerCorner()[1];
    double right = cuboid.rGetUpperCorner()[0];
    double top = cuboid.rGetUpperCorner()[1];
    double fudge_factor = 1e-6;

    // Apply Neumann boundary conditions
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = pMesh->GetBoundaryElementIteratorBegin();
         elem_iter != pMesh->GetBoundaryElementIteratorEnd();
         ++elem_iter)
    {
        double x_1 = (*elem_iter)->GetNodeLocation(0)[0];
//        double y_1 = (*elem_iter)->GetNodeLocation(0)[1];
        double x_2 = (*elem_iter)->GetNodeLocation(1)[0];
//        double y_2 = (*elem_iter)->GetNodeLocation(1)[1];

        assert(!are_neumann_boundaries[0]);
//        if (are_neumann_boundaries[0]) // Top is Neumann boundary
//        {
//            if ( (y_1  > (top-fudge_factor) ) &&  (y_2 > (top-fudge_factor) ))
//            {
//                p_bcc->AddNeumannBoundaryCondition(*elem_iter, mConstBoundaryConditions[0]);
//            }
//        }
        if (are_neumann_boundaries[1]) // Right is Neumann boundary
        {
            if ( (x_1  > (right-fudge_factor) ) &&  (x_2 > (right-fudge_factor) ))
            {
                p_bcc->AddNeumannBoundaryCondition(*elem_iter, mConstBoundaryConditions[1]);
            }
        }
        assert(!are_neumann_boundaries[2]);
//        if (are_neumann_boundaries[2]) // Bottom is Neumann boundary
//        {
//            if ( (y_1  < (bottom+fudge_factor) ) &&  (y_2 < (bottom+fudge_factor) ))
//            {
//                p_bcc->AddNeumannBoundaryCondition(*elem_iter, mConstBoundaryConditions[2]);
//            }
//        }
        if (are_neumann_boundaries[3]) // Left is Neumann Boundary
        {
            if ( (x_1  > (left+fudge_factor) ) &&  (x_2 < (left+fudge_factor) ))
            {
                p_bcc->AddNeumannBoundaryCondition(*elem_iter, mConstBoundaryConditions[3]);
            }
        }
    }

    // Apply Dirichlet boundaries
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = pMesh->GetBoundaryNodeIteratorBegin();
         node_iter != pMesh->GetBoundaryNodeIteratorEnd();
         ++node_iter)
    {
//        double x = (*node_iter)->GetPoint()[0];
        double y = (*node_iter)->GetPoint()[1];

        if (!are_neumann_boundaries[0]) // Top is Dirichlet boundary
        {
            if (y > top-fudge_factor)
            {
                p_bcc->AddDirichletBoundaryCondition(*node_iter, mConstBoundaryConditions[0]);
            }
        }
        assert(are_neumann_boundaries[1]);
//        if (!are_neumann_boundaries[1]) // Right is Dirichlet boundary
//        {
//            if (x > right-fudge_factor)
//            {
//                p_bcc->AddDirichletBoundaryCondition(*node_iter, mConstBoundaryConditions[1]);
//            }
//        }
        if (!are_neumann_boundaries[2]) // Bottom is Dirichlet boundary
        {
            if (y < bottom+fudge_factor)
            {
                p_bcc->AddDirichletBoundaryCondition(*node_iter, mConstBoundaryConditions[2]);
            }
        }
        assert(are_neumann_boundaries[3]);
//        if (!are_neumann_boundaries[3]) // Left is Dirichlet boundary
//        {
//            if (x < left+fudge_factor)
//            {
//                p_bcc->AddDirichletBoundaryCondition(*node_iter, mConstBoundaryConditions[3]);
//            }
//        }
    }

    return p_bcc;
}

template<unsigned DIM>
void CellBasedPdeHandlerOnCuboid<DIM>::OutputParameters(out_stream& rParamsFile)
{
    CellBasedPdeHandler<DIM>::OutputParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellBasedPdeHandlerOnCuboid)

///////// Explicit instantiation
template class CellBasedPdeHandlerOnCuboid<1>;
template class CellBasedPdeHandlerOnCuboid<2>;
template class CellBasedPdeHandlerOnCuboid<3>;
