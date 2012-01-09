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

#ifndef _TESTCELLBASEDPDESOLVER_HPP_
#define _TESTCELLBASEDPDESOLVER_HPP_

#include "UblasCustomFunctions.hpp"

#include <cxxtest/TestSuite.h>

#include <petsc.h>
#include <cmath>
#include <pde/test/pdes/SimplePoissonEquation.hpp>

#include "SimpleLinearEllipticSolver.hpp"
#include "CellBasedPdeSolver.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ConstBoundaryCondition.hpp"

class TestCellBasedPdeSolver : public CxxTest::TestSuite
{
public:

    void Test2dHeatEquationOnUnitSquare()
    {
        // Create mesh from mesh reader
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Instantiate PDE object
        SimplePoissonEquation<2,2> pde;

        // Boundary conditions
        BoundaryConditionsContainer<2,2,1> bcc;
        ConstBoundaryCondition<2>* p_boundary_condition = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(1), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(2), p_boundary_condition);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(3), p_boundary_condition);

        // Creat PDE solvers
        SimpleLinearEllipticSolver<2,2> simple_solver(&mesh, &pde, &bcc);
        CellBasedPdeSolver<2> pde_solver(&mesh, &pde, &bcc);

        Vec simple_result = simple_solver.Solve();
        Vec pde_result = pde_solver.Solve();

        ReplicatableVector simple_result_repl(simple_result);
        ReplicatableVector pde_result_repl(pde_result);

        for (unsigned i=0; i<simple_result_repl.GetSize(); i++)
        {
            TS_ASSERT_EQUALS(simple_result_repl[i], pde_result_repl[i]);
        }

        // Tidy up
        VecDestroy(simple_result);
        VecDestroy(pde_result);
    }
};

#endif //_TESTCELLBASEDPDESOLVER_HPP_
