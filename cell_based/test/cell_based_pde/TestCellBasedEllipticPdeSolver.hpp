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

#ifndef _TESTCELLBASEDELLIPTICPDESOLVER_HPP_
#define _TESTCELLBASEDELLIPTICPDESOLVER_HPP_

#include "UblasCustomFunctions.hpp"

#include <cxxtest/TestSuite.h>

#include <petsc.h>
#include <cmath>
#include <pde/test/pdes/SimplePoissonEquation.hpp>

#include "SimpleLinearEllipticSolver.hpp"
#include "CellBasedEllipticPdeSolver.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ConstBoundaryCondition.hpp"

class TestCellBasedEllipticPdeSolver : public CxxTest::TestSuite
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

        // Create PDE solvers
        SimpleLinearEllipticSolver<2,2> simple_solver(&mesh, &pde, &bcc);
        CellBasedEllipticPdeSolver<2> pde_solver(&mesh, &pde, &bcc);

        Vec simple_result = simple_solver.Solve();
        Vec pde_result = pde_solver.Solve();

        ReplicatableVector simple_result_repl(simple_result);
        ReplicatableVector pde_result_repl(pde_result);

        for (unsigned i=0; i<simple_result_repl.GetSize(); i++)
        {
            TS_ASSERT_EQUALS(simple_result_repl[i], pde_result_repl[i]);
        }

        // Tidy up
        PetscTools::Destroy(simple_result);
        PetscTools::Destroy(pde_result);
    }
};

#endif //_TESTCELLBASEDELLIPTICPDESOLVER_HPP_
