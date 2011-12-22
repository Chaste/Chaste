/*

Copyright (C) University of Oxford, 2005-2011

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

    void TestWithPoissonsEquation1dMeshIn2dSpace() throw (Exception)
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

        VecDestroy(result);
    }

    void TestWithPoissonsEquation1dMeshIn3dSpace() throw (Exception)
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

        VecDestroy(result);
    }

    void TestBranchedPoissonsEquation1dMeshIn3dSpace() throw (Exception)
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

        VecDestroy(result);
    }
};

#endif //_TESTSIMPLELINEARELLIPTICSOLVER_HPP_
