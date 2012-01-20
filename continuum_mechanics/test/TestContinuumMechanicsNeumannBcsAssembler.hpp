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

#ifndef TESTCONTINUUMMECHANICSNEUMANNBCSASSEMBLER_HPP_
#define TESTCONTINUUMMECHANICSNEUMANNBCSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "ContinuumMechanicsNeumannBcsAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestContinuumMechanicsNeumannBcsAssembler : public CxxTest::TestSuite
{
public:
    void TestAssembler2d() throw (Exception)
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
        Vec bad_sized_vec = PetscTools::CreateVec(3*mesh.GetNumNodes()+mesh.GetNumVertices());
        assembler.SetVectorToAssemble(bad_sized_vec, true);
        TS_ASSERT_THROWS_THIS(assembler.AssembleVector(), "Vector provided to be assembled has size 21, not expected size of 15 (dim*num_nodes+num_vertices)");


        Vec vec = PetscTools::CreateVec(2*mesh.GetNumNodes()+mesh.GetNumVertices());

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
        TS_ASSERT_DELTA(vec_repl[2], t1/6.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[3], t2/6.0, 1e-8);

        // nodes 2, 3, 4 are not on the surface
        for(unsigned i=2; i<5; i++)
        {
            TS_ASSERT_DELTA(vec_repl[2*i], 0.0, 1e-8);
            TS_ASSERT_DELTA(vec_repl[2*i], 0.0, 1e-8);
        }

        // node 5 is on the surface, and is an interior node
        TS_ASSERT_DELTA(vec_repl[10], 4*t1/6.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[11], 4*t2/6.0, 1e-8);


        // the rest of the vector is the pressure block and should be zero.
        TS_ASSERT_DELTA(vec_repl[12], 0.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[13], 0.0, 1e-8);
        TS_ASSERT_DELTA(vec_repl[14], 0.0, 1e-8);

        PetscTools::Destroy(vec);
        PetscTools::Destroy(bad_sized_vec);
    }


    void TestAssembler3d() throw (Exception)
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

        Vec vec = PetscTools::CreateVec(3*mesh.GetNumNodes()+mesh.GetNumVertices());

        ContinuumMechanicsNeumannBcsAssembler<3> assembler(&mesh, &problem_defn);
        assembler.SetVectorToAssemble(vec, true);
        assembler.AssembleVector();

        ReplicatableVector vec_repl(vec);

        // For boundary elements in 3d - ie for 2d elements - and with QUADRATIC basis functions, we have
        // \intgl_{canonical element} \phi_i dV = 0.0  for i=0,1,2 (ie bases on vertices)
        // and
        // \intgl_{canonical element} \phi_i dV = 1/6  for i=3,4,5 (ie bases on mid-nodes)

        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double z = mesh.GetNode(i)->rGetLocation()[2];

            if(fabs(z)<1e-8) // ie if z=0
            {
                if( fabs(x-0.5)<1e-8 || fabs(y-0.5)<1e-8 ) // if x=0.5 or y=0.5
                {
                    unsigned num_surf_elems_contained_in = 1;
                    if( fabs(x+y-1.0)<1e-8 ) // if x=0.5 AND y=0.5
                    {
                        num_surf_elems_contained_in = 2;
                    }

                    // interior node on traction surface
                    TS_ASSERT_DELTA(vec_repl[3*i],   num_surf_elems_contained_in * t1/6.0, 1e-8);
                    TS_ASSERT_DELTA(vec_repl[3*i+1], num_surf_elems_contained_in * t2/6.0, 1e-8);
                    TS_ASSERT_DELTA(vec_repl[3*i+2], num_surf_elems_contained_in * t3/6.0, 1e-8);
                }
                else
                {
                    TS_ASSERT_DELTA(vec_repl[3*i],   0.0, 1e-8);
                    TS_ASSERT_DELTA(vec_repl[3*i+1], 0.0, 1e-8);
                    TS_ASSERT_DELTA(vec_repl[3*i+2], 0.0, 1e-8);
                }
            }
            else
            {
                TS_ASSERT_DELTA(vec_repl[3*i],   0.0, 1e-8);
                TS_ASSERT_DELTA(vec_repl[3*i+1], 0.0, 1e-8);
                TS_ASSERT_DELTA(vec_repl[3*i+2], 0.0, 1e-8);
            }
        }

        for(unsigned i=0; i<mesh.GetNumVertices(); i++)
        {
            TS_ASSERT_DELTA(vec_repl[3*mesh.GetNumNodes()+i],   0.0, 1e-8);
        }

        PetscTools::Destroy(vec);
    }
};


#endif // TESTCONTINUUMMECHANICSNEUMANNBCSASSEMBLER_HPP_
