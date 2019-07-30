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

#ifndef TESTCONTINUUMMECHANICSNEUMANNBCSASSEMBLER_HPP_
#define TESTCONTINUUMMECHANICSNEUMANNBCSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ContinuumMechanicsNeumannBcsAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestContinuumMechanicsNeumannBcsAssembler : public CxxTest::TestSuite
{
public:
    void TestAssembler2d()
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/canonical_triangle_quadratic", 2, 2, false);
        mesh.ConstructFromMeshReader(mesh_reader);

        ContinuumMechanicsProblemDefinition<2> problem_defn(mesh);

        double t1 = 2.6854233;
        double t2 = 3.2574578;

        // for the boundary element on the y=0 surface, create a traction
        std::vector<BoundaryElement<1,2>*> boundary_elems;
        std::vector<c_vector<double,2> > tractions;
        c_vector<double,2> traction;
        traction(0) = t1;
        traction(1) = t2;
        for (TetrahedralMesh<2,2>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[1])<1e-4)
            {
                BoundaryElement<1,2>* p_element = *iter;
                boundary_elems.push_back(p_element);
                tractions.push_back(traction);
            }
        }
        assert(boundary_elems.size()==1);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);


        ContinuumMechanicsNeumannBcsAssembler<2> assembler(&mesh, &problem_defn);

        TS_ASSERT_THROWS_THIS(assembler.AssembleVector(), "Vector to be assembled has not been set");
        Vec bad_sized_vec = PetscTools::CreateVec(2);
        assembler.SetVectorToAssemble(bad_sized_vec, true);
        TS_ASSERT_THROWS_THIS(assembler.AssembleVector(), "Vector provided to be assembled has size 2, not expected size of 18 ((dim+1)*num_nodes)");


        Vec vec = PetscTools::CreateVec(3*mesh.GetNumNodes());

        assembler.SetVectorToAssemble(vec, true);
        assembler.AssembleVector();

        ReplicatableVector vec_repl(vec);

        // Note: on a 1d boundary element,  intgl phi_i dx = 1/6 for the bases on the vertices
        // and intgl phi_i dx = 4/6 for the basis at the interior node. (Here the integrals
        // are over the canonical 1d element, [0,1], which is also the physical element for this
        // mesh.

        // node 0 is on the surface, and is a vertex
        TS_ASSERT_DELTA(vec_repl[0], t1/6.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[1], t2/6.0, 1e-8);

        // node 1 is on the surface, and is a vertex
        TS_ASSERT_DELTA(vec_repl[3], t1/6.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[4], t2/6.0, 1e-8);

        // nodes 2, 3, 4 are not on the surface
        for (unsigned i=2; i<5; i++)
        {
            TS_ASSERT_DELTA(vec_repl[3*i], 0.0, 1e-8);
            TS_ASSERT_DELTA(vec_repl[3*i], 0.0, 1e-8);
        }

        // node 5 is on the surface, and is an interior node
        TS_ASSERT_DELTA(vec_repl[15], 4*t1/6.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[16], 4*t2/6.0, 1e-8);


        // the rest of the vector is the pressure block and should be zero.
        TS_ASSERT_DELTA(vec_repl[2], 0.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[5], 0.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[8], 0.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[11], 0.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[14], 0.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[17], 0.0, 1e-8);

        PetscTools::Destroy(vec);
        PetscTools::Destroy(bad_sized_vec);
    }


    void TestAssembler3d()
    {
        QuadraticMesh<3> mesh(1.0, 1.0, 1.0, 1.0);

        ContinuumMechanicsProblemDefinition<3> problem_defn(mesh);

        double t1 = 2.6854233;
        double t2 = 3.2574578;
        double t3 = 4.5342308;

        // for the boundary element on the z=0 surface, create a traction
        std::vector<BoundaryElement<2,3>*> boundary_elems;
        std::vector<c_vector<double,3> > tractions;
        c_vector<double,3> traction;
        traction(0) = t1;
        traction(1) = t2;
        traction(2) = t3;

        for (TetrahedralMesh<3,3>::BoundaryElementIterator iter
              = mesh.GetBoundaryElementIteratorBegin();
             iter != mesh.GetBoundaryElementIteratorEnd();
             ++iter)
        {
            if (fabs((*iter)->CalculateCentroid()[2])<1e-4)
            {
                BoundaryElement<2,3>* p_element = *iter;
                boundary_elems.push_back(p_element);
                tractions.push_back(traction);
            }
        }
        assert(boundary_elems.size()==2);
        problem_defn.SetTractionBoundaryConditions(boundary_elems, tractions);

        Vec vec = PetscTools::CreateVec(4*mesh.GetNumNodes());

        ContinuumMechanicsNeumannBcsAssembler<3> assembler(&mesh, &problem_defn);
        assembler.SetVectorToAssemble(vec, true);
        assembler.AssembleVector();

        ReplicatableVector vec_repl(vec);

        // For boundary elements in 3d - ie for 2d elements - and with QUADRATIC basis functions, we have
        // \intgl_{canonical element} \phi_i dV = 0.0  for i=0,1,2 (ie bases on vertices)
        // and
        // \intgl_{canonical element} \phi_i dV = 1/6  for i=3,4,5 (ie bases on mid-nodes)

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double z = mesh.GetNode(i)->rGetLocation()[2];

            if (fabs(z) < 1e-8) // ie if z=0
            {
                if (fabs(x-0.5)<1e-8 || fabs(y-0.5)<1e-8) // if x=0.5 or y=0.5
                {
                    unsigned num_surf_elems_contained_in = 1;
                    if (fabs(x+y-1.0) < 1e-8) // if x=0.5 AND y=0.5
                    {
                        num_surf_elems_contained_in = 2;
                    }

                    // interior node on traction surface
                    TS_ASSERT_DELTA(vec_repl[4*i],   num_surf_elems_contained_in * t1/6.0, 1e-8);
                    TS_ASSERT_DELTA(vec_repl[4*i+1], num_surf_elems_contained_in * t2/6.0, 1e-8);
                    TS_ASSERT_DELTA(vec_repl[4*i+2], num_surf_elems_contained_in * t3/6.0, 1e-8);
                }
                else
                {
                    TS_ASSERT_DELTA(vec_repl[4*i],   0.0, 1e-8);
                    TS_ASSERT_DELTA(vec_repl[4*i+1], 0.0, 1e-8);
                    TS_ASSERT_DELTA(vec_repl[4*i+2], 0.0, 1e-8);
                }
            }
            else
            {
                TS_ASSERT_DELTA(vec_repl[4*i],   0.0, 1e-8);
                TS_ASSERT_DELTA(vec_repl[4*i+1], 0.0, 1e-8);
                TS_ASSERT_DELTA(vec_repl[4*i+2], 0.0, 1e-8);
            }
        }

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(vec_repl[4*i+3],   0.0, 1e-8);
        }

        PetscTools::Destroy(vec);
    }

    void TestAssemblerMeshType()
    {
        TetrahedralMesh<2,2> mesh;
        ContinuumMechanicsProblemDefinition<2> problem_defn(mesh);

        TS_ASSERT_THROWS_CONTAINS(ContinuumMechanicsNeumannBcsAssembler<2>(&mesh, &problem_defn),
                                  "Continuum mechanics solvers require a quadratic mesh");

        TetrahedralMesh<3,3> mesh3d;
        ContinuumMechanicsProblemDefinition<3> problem_defn3d(mesh3d);

        TS_ASSERT_THROWS_CONTAINS(ContinuumMechanicsNeumannBcsAssembler<3>(&mesh3d, &problem_defn3d),
                                  "Continuum mechanics solvers require a quadratic mesh");
    }
};


#endif // TESTCONTINUUMMECHANICSNEUMANNBCSASSEMBLER_HPP_
