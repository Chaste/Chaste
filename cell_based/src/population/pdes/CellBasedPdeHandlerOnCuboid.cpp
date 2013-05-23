/*

Copyright (c) 2005-2013, University of Oxford.
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
}

template<unsigned DIM>
BoundaryConditionsContainer<DIM,DIM,1> CellBasedPdeHandlerOnCuboid<DIM>::ConstructBoundaryConditionsContiner(PdeAndBoundaryConditions<DIM>* pPdeAndBc,TetrahedralMesh<DIM,DIM>* pMesh)
{
    // Not using the inputs as theres only one BCS
    assert(DIM==2);

    BoundaryConditionsContainer<DIM,DIM,1> bcc(false);


    // Use these 2 vectors to define whats happening on the top right bottom and left boundaries

    // Specify which sides are neuman boundaries
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


    ConstBoundaryCondition<DIM>* p_bc_top = new ConstBoundaryCondition<DIM>(boundary_condition_values[0]);
    ConstBoundaryCondition<DIM>* p_bc_right = new ConstBoundaryCondition<DIM>(boundary_condition_values[1]);
    ConstBoundaryCondition<DIM>* p_bc_bottom = new ConstBoundaryCondition<DIM>(boundary_condition_values[2]);
    ConstBoundaryCondition<DIM>* p_bc_left = new ConstBoundaryCondition<DIM>(boundary_condition_values[3]);

    ChasteCuboid<DIM> cuboid = pMesh->CalculateBoundingBox();

    double left = cuboid.rGetLowerCorner()[0];
    double bottom = cuboid.rGetLowerCorner()[1];
    double right = cuboid.rGetUpperCorner()[0];
    double top = cuboid.rGetUpperCorner()[1];
    double fudge_factor = 1e-6;


    // Apply Neumann Boundaries
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = pMesh->GetBoundaryElementIteratorBegin();
         elem_iter != pMesh->GetBoundaryElementIteratorEnd();
         ++elem_iter)
    {
        double x_1 = (*elem_iter)->GetNodeLocation(0)[0];
        double y_1 = (*elem_iter)->GetNodeLocation(0)[1];
        double x_2 = (*elem_iter)->GetNodeLocation(1)[0];
        double y_2 = (*elem_iter)->GetNodeLocation(1)[1];

        if (are_neumann_boundaries[0]) // Top is Neumann Boundary
        {
            if ( (y_1  > (top-fudge_factor) ) &&  (y_2 > (top-fudge_factor) ))
            {
                bcc.AddNeumannBoundaryCondition(*elem_iter, p_bc_top);
            }
        }
        if (are_neumann_boundaries[1]) // Right is Neumann Boundary
        {
            if ( (x_1  > (right-fudge_factor) ) &&  (x_2 > (right-fudge_factor) ))
            {
                bcc.AddNeumannBoundaryCondition(*elem_iter, p_bc_right);
            }
        }
        if (are_neumann_boundaries[2]) // Bottom is Neumann Boundary
        {
            if ( (y_1  < (bottom+fudge_factor) ) &&  (y_2 < (bottom+fudge_factor) ))
            {
                bcc.AddNeumannBoundaryCondition(*elem_iter, p_bc_bottom);
            }
        }
        if (are_neumann_boundaries[3]) // Left is Neumann Boundary
        {
            if ( (x_1  > (left+fudge_factor) ) &&  (x_2 < (left+fudge_factor) ))
            {
                bcc.AddNeumannBoundaryCondition(*elem_iter, p_bc_left);
            }
        }
    }


    // Apply dirichlet Boundaries
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = pMesh->GetBoundaryNodeIteratorBegin();
         node_iter != pMesh->GetBoundaryNodeIteratorEnd();
         ++node_iter)
    {
        double x = (*node_iter)->GetPoint()[0];
        double y = (*node_iter)->GetPoint()[1];

        if (!are_neumann_boundaries[0]) // Top is Dirichlet Boundary
        {
            if (y > top-fudge_factor)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, p_bc_top);
            }
        }
        if (!are_neumann_boundaries[1]) // Right is Dirichlet Boundary
        {
            if (x > right-fudge_factor)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, p_bc_right);
            }
        }
        if (!are_neumann_boundaries[2]) // Bottom is Dirichlet Boundary
        {
            if (y < bottom+fudge_factor)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, p_bc_bottom);
            }
        }
        if (!are_neumann_boundaries[3]) // Left is Dirichlet Boundary
        {
            if (x < left+fudge_factor)
            {
                bcc.AddDirichletBoundaryCondition(*node_iter, p_bc_left);
            }
        }
    }

    return bcc;
}

template<unsigned DIM>
void CellBasedPdeHandlerOnCuboid<DIM>::OutputParameters(out_stream& rParamsFile)
{
    CellBasedPdeHandler<DIM>::OutputParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellBasedPdeHandlerOnCuboid)

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CellBasedPdeHandlerOnCuboid<1>;
template class CellBasedPdeHandlerOnCuboid<2>;
template class CellBasedPdeHandlerOnCuboid<3>;
