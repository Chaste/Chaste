/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
    void TestSurfaceElementContributions() throw(Exception)
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

        VecDestroy(vec);
    }


    // Test surface element intregral additions in 2d
    void TestSurfaceElementContributions2d() throw(Exception)
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

        VecDestroy(vec);
    }
};

#endif // TESTABSTRACTFESURFACEINTEGRALASSEMBLER_HPP_
