/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef _TESTTISSUESIMULATIONWITHNUTRIENTSASSEMBLER_HPP_
#define _TESTTISSUESIMULATIONWITHNUTRIENTSASSEMBLER_HPP_

#include "UblasCustomFunctions.hpp"

#include <cxxtest/TestSuite.h>

#include <petsc.h>
#include <cmath>
#include <pde/test/pdes/SimplePoissonEquation.hpp>

#include "SimpleLinearEllipticAssembler.hpp"
#include "TissueSimulationWithNutrientsAssembler.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ConstBoundaryCondition.hpp"


class TestTissueSimulationWithNutrientsAssembler : public CxxTest::TestSuite
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

        // Assembler
        SimpleLinearEllipticAssembler<2,2> simple_assembler(&mesh,&pde,&bcc);
        TissueSimulationWithNutrientsAssembler<2> nutrients_assembler(&mesh,&pde,&bcc);

        Vec simple_result = simple_assembler.Solve();
        Vec nutrients_result = nutrients_assembler.Solve();

        ReplicatableVector simple_result_repl(simple_result);
        ReplicatableVector nutrients_result_repl(nutrients_result);

        for (unsigned i=0; i<simple_result_repl.GetSize(); i++)
        {
            TS_ASSERT_EQUALS(simple_result_repl[i], nutrients_result_repl[i]);
        }

        VecDestroy(simple_result);
        VecDestroy(nutrients_result);
    }

};

#endif //_TESTTISSUESIMULATIONWITHNUTRIENTSASSEMBLER_HPP_
