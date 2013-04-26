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
    assert(DIM==2);

    BoundaryConditionsContainer<DIM,DIM,1> bcc(false);


    // Not using the inputs as theres only one BCS
    AbstractBoundaryCondition<DIM>* p_bc = pPdeAndBc->GetBoundaryCondition();

    ConstBoundaryCondition<DIM> bc_left(0.0);
    ConstBoundaryCondition<DIM> bc_right(0.5);
    ConstBoundaryCondition<DIM> bc_top(1.0);
    ConstBoundaryCondition<DIM> bc_bottom(-1.0);

    ChasteCuboid<DIM> cuboid = pMesh->CalculateBoundingBox();

    double left = cuboid.rGetLowerCorner()[0];
    double bottom = cuboid.rGetLowerCorner()[1];
    double right = cuboid.rGetUpperCorner()[0];
    double top = cuboid.rGetUpperCorner()[1];


    // If on right or left apply Neuman BCS
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = pMesh->GetBoundaryElementIteratorBegin();
         elem_iter != pMesh->GetBoundaryElementIteratorEnd();
         ++elem_iter)
    {
        if ( ((*elem_iter)->GetNodeLocation(0)[0] < (left+1e-6) ) &&  ((*elem_iter)->GetNodeLocation(1)[0] < (left+1e-6) ))
        {
            bcc.AddNeumannBoundaryCondition(*elem_iter, p_bc);
        }
        if ( ((*elem_iter)->GetNodeLocation(0)[0] > (right-1e-6) ) &&  ((*elem_iter)->GetNodeLocation(1)[0] > (right-1e-6) ))
        {
            bcc.AddNeumannBoundaryCondition(*elem_iter, p_bc);
        }
    }


    // If on Top or bottom then apply Dirichlet BCS
    for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = pMesh->GetBoundaryNodeIteratorBegin();
         node_iter != pMesh->GetBoundaryNodeIteratorEnd();
         ++node_iter)
    {
        double y = (*node_iter)->GetPoint()[1];

        if (y<bottom+1e-6)
        {
            bcc.AddDirichletBoundaryCondition(*node_iter, p_bc);
        }
        if (y>top-1e-6)
        {
            bcc.AddDirichletBoundaryCondition(*node_iter, p_bc);
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
