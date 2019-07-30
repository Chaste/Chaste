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

#ifndef _TESTSIMPLELINEARELLIPTICSOLVER_HPP_
#define _TESTSIMPLELINEARELLIPTICSOLVER_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include "SimplePoissonEquation.hpp"
#include "SimpleLinearEllipticSolver.hpp"
#include <vector>
#include <cmath>
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "TrianglesMeshReader.hpp"

/*
 * These are need for the nD problems in mD space (n!=m), as those
 * particular cases are not explicitly instantiated.
 */
#include "AbstractBoundaryConditionsContainerImplementation.hpp"
#include "BoundaryConditionsContainerImplementation.hpp"

#include "PetscSetupAndFinalize.hpp"

class Test1D3DLinearEllipticSolver : public CxxTest::TestSuite
{
public:

    void TestWithPoissonsEquation1dMeshIn2dSpace()
    {
        const unsigned SPACE_DIM = 2;
        const unsigned ELEMENT_DIM = 1;

        // Create mesh from mesh reader
        TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader("mesh/test/data/trivial_1d_in_2d_mesh");
        TetrahedralMesh<ELEMENT_DIM,SPACE_DIM> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<ELEMENT_DIM,SPACE_DIM> pde;

        // Boundary conditions (u=0 on one end, u'=0 on other end)
        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1> bcc;
        ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition = new ConstBoundaryCondition<SPACE_DIM>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);

        // Solver
        SimpleLinearEllipticSolver<ELEMENT_DIM,SPACE_DIM> assembler(&mesh,&pde,&bcc);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = 0.5*x*(3-x)
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = 0.5*x*(3-x);
            TS_ASSERT_DELTA(result_repl[i], u, 0.001);
        }

        PetscTools::Destroy(result);
    }

    void TestWithPoissonsEquation1dMeshIn3dSpace()
    {
        const unsigned SPACE_DIM = 3;
        const unsigned ELEMENT_DIM = 1;

        // Create mesh from mesh reader
        TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader("mesh/test/data/trivial_1d_in_3d_mesh");
        TetrahedralMesh<ELEMENT_DIM,SPACE_DIM> mesh;

        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<ELEMENT_DIM,SPACE_DIM> pde;

        // Boundary conditions (u=0 on one end, u'=0 on other end)
        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1> bcc;
        ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition = new ConstBoundaryCondition<SPACE_DIM>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);

        // Solver
        SimpleLinearEllipticSolver<ELEMENT_DIM,SPACE_DIM> assembler(&mesh,&pde,&bcc);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = 0.5*x*(3-x)
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            double x = mesh.GetNode(i)->GetPoint()[0];
            double u = 0.5*x*(3-x);
            TS_ASSERT_DELTA(result_repl[i], u, 0.001);
        }

        PetscTools::Destroy(result);
    }

    void TestBranchedPoissonsEquation1dMeshIn3dSpace()
    {
        const unsigned SPACE_DIM = 3;
        const unsigned ELEMENT_DIM = 1;

        // Create mesh from mesh reader
        TrianglesMeshReader<ELEMENT_DIM,SPACE_DIM> mesh_reader("mesh/test/data/branched_1d_in_3d_mesh");
        TetrahedralMesh<ELEMENT_DIM,SPACE_DIM> mesh;

        mesh.ConstructFromMeshReader(mesh_reader);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryElements(), 3u);
        TS_ASSERT_EQUALS(mesh.GetNumBoundaryNodes(), 3u);

        // Instantiate PDE object
        SimplePoissonEquation<ELEMENT_DIM,SPACE_DIM> pde;

        // Boundary conditions (u=0 on one end, u'=0 on other end)
        BoundaryConditionsContainer<ELEMENT_DIM,SPACE_DIM,1> bcc;
        ConstBoundaryCondition<SPACE_DIM>* p_boundary_condition = new ConstBoundaryCondition<SPACE_DIM>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(20), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(30), p_boundary_condition);

        // Solver
        SimpleLinearEllipticSolver<ELEMENT_DIM,SPACE_DIM> assembler(&mesh,&pde,&bcc);

        Vec result = assembler.Solve();
        ReplicatableVector result_repl(result);

        // Solution should be u = -r^2/2 + r, where r is the distance along the branches
        for (unsigned i=0; i<result_repl.GetSize(); i++)
        {
            //This only works because the x,y,z coordinates are aligned to the axes
            //should really calculate distances along the branches properly!
            double x = mesh.GetNode(i)->GetPoint()[0];
            double y = mesh.GetNode(i)->GetPoint()[1];
            double z = mesh.GetNode(i)->GetPoint()[2];
            double r = x+y+z;
            double u = -r*r/2 + r;
            TS_ASSERT_DELTA(result_repl[i], u, 0.001);
        }

        PetscTools::Destroy(result);
    }
};

#endif //_TESTSIMPLELINEARELLIPTICSOLVER_HPP_
