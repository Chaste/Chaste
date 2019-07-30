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

#ifndef TESTABSTRACTFESURFACEINTEGRALASSEMBLER_HPP_
#define TESTABSTRACTFESURFACEINTEGRALASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractFeSurfaceIntegralAssembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"

template<unsigned DIM>
class BasicSurfaceAssembler : public AbstractFeSurfaceIntegralAssembler<DIM,DIM,1>
{
private:
    c_vector<double, 1*DIM> ComputeVectorSurfaceTerm(
        const BoundaryElement<DIM-1,DIM>& rSurfaceElement,
        c_vector<double, DIM>& rPhi,
        ChastePoint<DIM>& rX)
    {
        c_vector<double,1*(DIM)> vec;
        for (unsigned i=0; i<DIM; i++)
        {
            vec(i)=1.0;
        }
        return vec;
    }

public:
    BasicSurfaceAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh, BoundaryConditionsContainer<DIM,DIM,1>* pBcc)
        : AbstractFeSurfaceIntegralAssembler<DIM,DIM,1>(pMesh,pBcc)
    {
    }
};

class TestAbstractFeSurfaceIntegralAssembler : public CxxTest::TestSuite
{
public:
    // Test surface element intregral additions in 1d
    void TestSurfaceElementContributions()
    {
        TetrahedralMesh<1,1> mesh;
        double h = 0.1;
        mesh.ConstructRegularSlabMesh(h, 1.0);

        // create a BCC. We just want to tell the assembler that there is a Neumann
        // bc at each end node, doesn't matter what the bc is as the concrete class
        // doesn't use it
        BoundaryConditionsContainer<1,1,1> bcc;
        ConstBoundaryCondition<1>* p_boundary_condition = new ConstBoundaryCondition<1>(1.0);
        TetrahedralMesh<1,1>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
        iter++;
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);
        iter++;
        assert(iter==mesh.GetBoundaryElementIteratorEnd());
        // initialise the vector to (a,a,...,a);
        double offset = 0.342423432;
        Vec vec = PetscTools::CreateAndSetVec(mesh.GetNumNodes(),offset);

        // create the assembler, this time saying don't zero before assembling and
        // apply surface integrals..
        BasicSurfaceAssembler<1> assembler(&mesh,&bcc);
        assembler.SetVectorToAssemble(vec,false); // don't zero before assembling
        assembler.Assemble();

        PetscVecTools::Finalise(vec);

        ReplicatableVector vec_repl(vec);

        TS_ASSERT_DELTA(vec_repl[0], 1.0 + offset, 1e-4);
        TS_ASSERT_DELTA(vec_repl[mesh.GetNumNodes()-1], 1.0 + offset, 1e-4);

        for (unsigned i=1; i<mesh.GetNumNodes()-1; i++)
        {
            TS_ASSERT_DELTA(vec_repl[i], 0.0 + offset, 1e-4);
        }

        PetscTools::Destroy(vec);
    }


    // Test surface element intregral additions in 2d
    void TestSurfaceElementContributions2d()
    {
        // two element mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRegularSlabMesh(1.0, 1.0, 1.0);

        // create a BCC with ONE neumann bc
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_boundary_condition = new ConstBoundaryCondition<2>(1.0);
        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorBegin();
        bcc.AddNeumannBoundaryCondition(*iter, p_boundary_condition);

        //// The following shows the nodes 2 and 3 make up the boundary element
        //std::cout << (*iter)->GetNodeGlobalIndex(0) << " " << (*iter)->GetNodeGlobalIndex(1) << "\n";

        Vec vec = PetscTools::CreateVec(mesh.GetNumNodes());

        BasicSurfaceAssembler<2> assembler(&mesh,&bcc);
        assembler.SetVectorToAssemble(vec,true);
        assembler.Assemble();

        PetscVecTools::Finalise(vec);

        ReplicatableVector vec_repl(vec);

        // nodes 2 and 3 should have 1 (= area-of-surface-elem) added on
        TS_ASSERT_DELTA(vec_repl[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(vec_repl[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(vec_repl[2], 1.0, 1e-4);
        TS_ASSERT_DELTA(vec_repl[3], 1.0, 1e-4);

        PetscTools::Destroy(vec);
    }
};

#endif // TESTABSTRACTFESURFACEINTEGRALASSEMBLER_HPP_
