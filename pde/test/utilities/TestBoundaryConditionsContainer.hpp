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

#ifndef _TESTBOUNDARYCONDITIONCONTAINER_HPP_
#define _TESTBOUNDARYCONDITIONCONTAINER_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "ReplicatableVector.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestBoundaryConditionsContainer : public CxxTest::TestSuite
{
public:
    void TestSetGet()
    {
        // Test in 1d

        int num_nodes = 10;
        BoundaryConditionsContainer<1,1,1> bcc1;

        TS_ASSERT(!bcc1.HasDirichletBoundaryConditions());

        Node<1>* nodes[num_nodes];
        for (int i=0; i<num_nodes; i++)
        {
            nodes[i] = new Node<1>(i,true,0);
            ConstBoundaryCondition<1>* p_boundary_condition =
                new ConstBoundaryCondition<1>((double)i);
            bcc1.AddDirichletBoundaryCondition(nodes[i], p_boundary_condition);
        }
        bcc1.ResetDirichletCommunication();
        TS_ASSERT(bcc1.HasDirichletBoundaryConditions());

        for (int i=0; i<num_nodes; i++)
        {
            double value = bcc1.GetDirichletBCValue(nodes[i]);
            TS_ASSERT_DELTA( value, i, 1e-12 );
        }

        for (int i=0; i<num_nodes; i++)
        {
            delete nodes[i];
        }

        int num_elem = 10;
        std::vector<BoundaryElement<0,1> > elements;
        for (unsigned element_index=0; element_index< (unsigned) num_elem; element_index++)
        {
            std::vector<Node<1>* > nodes;
            Node<1>* node = new Node<1>(element_index,true,0);
            nodes.push_back(node);

            BoundaryElement<0,1> element(element_index, nodes);
            elements.push_back(element);
        }
        for (int i=0; i<num_elem; i++)
        {
            ConstBoundaryCondition<1>* p_boundary_condition =
                new ConstBoundaryCondition<1>((double)i);
            bcc1.AddNeumannBoundaryCondition(&elements[i], p_boundary_condition);
        }

        for (int i=0; i<num_elem; i++)
        {
            double value = bcc1.GetNeumannBCValue(&elements[i], elements[i].GetNode(0)->GetIndex() );
            TS_ASSERT_DELTA( value, i, 1e-12 );
            delete elements[i].GetNode(0);
        }

        // Test in 2d

        num_nodes = 10;
        BoundaryConditionsContainer<2,2,1> bcc2;

        Node<2>* nodes2[num_nodes];
        for (int i=0; i<num_nodes; i++)
        {
            nodes2[i] = new Node<2>(i,true,0,0);
            ConstBoundaryCondition<2>* p_boundary_condition =
                new ConstBoundaryCondition<2>((double)i);
            bcc2.AddDirichletBoundaryCondition(nodes2[i], p_boundary_condition);
        }

        for (int i=0; i<num_nodes; i++)
        {
            double value = bcc2.GetDirichletBCValue(nodes2[i]);
            TS_ASSERT_DELTA( value, i, 1e-12 );
        }

        for (int i=0; i<num_nodes; i++)
        {
            delete nodes2[i];
        }

        num_elem = 10;
        std::vector<BoundaryElement<1,2> > elements2;
        for (unsigned element_index=0; element_index< (unsigned) num_elem; element_index++)
        {
            std::vector<Node<2>* > nodes;
            Node<2>* node0 = new Node<2>(element_index,true,0,0);
            Node<2>* node1 = new Node<2>(element_index,true,1,1);

            nodes.push_back(node0);
            nodes.push_back(node1);
            BoundaryElement<1,2> element(element_index, nodes);

            elements2.push_back(element);
        }
        for (int i=0; i<num_elem; i++)
        {
            ConstBoundaryCondition<2>* p_boundary_condition =
                new ConstBoundaryCondition<2>((double)i);
            bcc2.AddNeumannBoundaryCondition(&elements2[i], p_boundary_condition);
        }

        for (int i=0; i<num_elem; i++)
        {
            double value = bcc2.GetNeumannBCValue(&elements2[i], elements2[i].GetNode(0)->GetIndex() );
            TS_ASSERT_DELTA( value, i, 1e-12 );
            delete elements2[i].GetNode(0);
            delete elements2[i].GetNode(1);
        }

        // Test in 3d

        num_nodes = 10;
        BoundaryConditionsContainer<3,3,1> bcc3;

        Node<3>* nodes3[num_nodes];
        for (int i=0; i<num_nodes; i++)
        {
            nodes3[i] = new Node<3>(i,true,0,0);
            ConstBoundaryCondition<3>* p_boundary_condition =
                new ConstBoundaryCondition<3>((double)i);
            bcc3.AddDirichletBoundaryCondition(nodes3[i], p_boundary_condition);
        }

        for (int i=0; i<num_nodes; i++)
        {
            double value = bcc3.GetDirichletBCValue(nodes3[i]);
            TS_ASSERT_DELTA( value, i, 1e-12 );
        }
        for (int i=0; i<num_nodes; i++)
        {
            delete nodes3[i];
        }

        num_elem = 10;
        std::vector<BoundaryElement<2,3> > elements3;
        for (int element_index=0; element_index<num_elem; element_index++)
        {
            std::vector<Node<3>* > nodes;
            Node<3>* node0 = new Node<3>(element_index,true,0,0,0);
            Node<3>* node1 = new Node<3>(element_index,true,1,0,0);
            Node<3>* node2 = new Node<3>(element_index,true,0,1,0);
            nodes.push_back(node0);
            nodes.push_back(node1);
            nodes.push_back(node2);
            BoundaryElement<2,3> element(element_index, nodes);

            elements3.push_back(element);
        }
        for (int i=0; i<num_elem; i++)
        {
            ConstBoundaryCondition<3>* p_boundary_condition =
                new ConstBoundaryCondition<3>((double)i);
            bcc3.AddNeumannBoundaryCondition(&elements3[i], p_boundary_condition);
        }

        for (int i=0; i<num_elem; i++)
        {
            double value = bcc3.GetNeumannBCValue(&elements3[i], elements3[i].GetNode(0)->GetIndex() );
            TS_ASSERT_DELTA( value, i, 1e-12 );
            delete elements3[i].GetNode(0);
            delete elements3[i].GetNode(1);
            delete elements3[i].GetNode(2);
        }
    }

    void TestApplyToLinearSystem()
    {
        const int SIZE = 10;
        LinearSystem some_system(SIZE, SIZE);
        for (int i=0; i<SIZE; i++)
        {
            for (int j=0; j<SIZE; j++)
            {
                // LHS matrix is all 1s
                some_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            some_system.SetRhsVectorElement(i,2);
        }

        some_system.AssembleIntermediateLinearSystem();

        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,1> bcc3;

        // Apply dirichlet boundary conditions to all but last node
        for (int i=0; i<SIZE-1; i++)
        {
            nodes_array[i] = new Node<3>(i,true);
            ConstBoundaryCondition<3>* p_boundary_condition =
                new ConstBoundaryCondition<3>(-1);
            bcc3.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition);
        }

        //////////////////////////
        // 2010 AD code from here
        //////////////////////////

        // apply dirichlet bcs to matrix but not rhs vector
        bcc3.ApplyDirichletToLinearProblem(some_system, true, false);
        ReplicatableVector vec_repl(some_system.GetRhsVector());
        for (unsigned i=0; i<(unsigned)SIZE; i++)
        {
            TS_ASSERT_EQUALS(vec_repl[i], 2.0);
        }

        // now apply to the rhs vector
        bcc3.ApplyDirichletToLinearProblem(some_system, false, true);

        some_system.AssembleFinalLinearSystem();

        //////////////////////////
        // 2007 AD code from here
        //////////////////////////

        /*
         *  Based on the original system and the boundary conditions applied in a non-symmetric
         *  manner, the resulting system looks like:
         *
         *      1 0 0 ... 0
         *      0 1 0 ... 0
         *      0 0 1 ... 0
         *      ...
         *      1 1 1 ... 1
         *
         */
        /// \todo: this is very naughty. Must be checked in parallel as well.
        if (PetscTools::IsSequential())
        {
            for (int row=0; row<SIZE-1; row++)
            {
                for (int column=0; column<row; column++)
                {
                    TS_ASSERT_EQUALS(some_system.GetMatrixElement(row,column), 0);
                }

                TS_ASSERT_EQUALS(some_system.GetMatrixElement(row,row), 1);

                for (int column=row+1; column<SIZE; column++)
                {
                    TS_ASSERT_EQUALS(some_system.GetMatrixElement(row,column), 0);
                }
            }

            for (int column=0; column<SIZE; column++)
            {
                TS_ASSERT_EQUALS(some_system.GetMatrixElement(SIZE-1,column), 1);
            }
        }

        Vec solution = some_system.Solve();

        DistributedVectorFactory factory(solution);
        DistributedVector d_solution = factory.CreateDistributedVector( solution );
        for (DistributedVector::Iterator index = d_solution.Begin();
             index != d_solution.End();
             ++index)
        {
            double expected = index.Global < SIZE-1 ? -1.0 : 11.0;
            TS_ASSERT_DELTA(d_solution[index], expected, 1e-6 );
        }

        for (int i=0; i<SIZE-1; i++)
        {
            delete nodes_array[i];
        }
        PetscTools::Destroy(solution);
    }

    void TestApplyToSymmetricLinearSystem()
    {
        const int SIZE = 10;
        LinearSystem some_system(SIZE, SIZE);
        some_system.SetMatrixIsSymmetric(true);

        for (int i=0; i<SIZE; i++)
        {
            for (int j=0; j<SIZE; j++)
            {
                // LHS matrix is all 1s
                some_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            some_system.SetRhsVectorElement(i,2);
        }

        some_system.AssembleIntermediateLinearSystem();

        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,1> bcc3;

        // Apply dirichlet boundary conditions to all but last node
        for (int i=0; i<SIZE-1; i++)
        {
            nodes_array[i] = new Node<3>(i,true);
            ConstBoundaryCondition<3>* p_boundary_condition =
                new ConstBoundaryCondition<3>(-1);
            bcc3.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition);
        }
        bcc3.ApplyDirichletToLinearProblem(some_system);

        some_system.AssembleFinalLinearSystem();

        /*
         * Based on the original system and the boundary conditions applied in a symmetric
         * manner, the resulting system looks like:
         *
         *      1 0 0 ... 0
         *      0 1 0 ... 0
         *      0 0 1 ... 0
         *      ...
         *      0 0 0 ... 1
         */
        int lo, hi;
        some_system.GetOwnershipRange(lo, hi);

        for (int row=lo; row<hi; row++)
        {
            for (int column=0; column<row; column++)
            {
                TS_ASSERT_EQUALS(some_system.GetMatrixElement(row,column), 0);
            }

            TS_ASSERT_EQUALS(some_system.GetMatrixElement(row,row), 1);

            for (int column=row+1; column<SIZE; column++)
            {
                TS_ASSERT_EQUALS(some_system.GetMatrixElement(row,column), 0);
            }
        }

        Vec solution = some_system.Solve();

        DistributedVectorFactory factory(solution);
        DistributedVector d_solution = factory.CreateDistributedVector( solution );
        for (DistributedVector::Iterator index = d_solution.Begin();
             index != d_solution.End();
             ++index)
        {
            double expected = index.Global < SIZE-1 ? -1.0 : 11.0;
            TS_ASSERT_DELTA(d_solution[index], expected, 1e-6 );
        }

        for (int i=0; i<SIZE-1; i++)
        {
            delete nodes_array[i];
        }
        PetscTools::Destroy(solution);
    }

    void TestApplyToNonlinearSystem()
    {
        const int SIZE = 10;
        DistributedVectorFactory factory(SIZE);

        Vec solution = factory.CreateVec();
        DistributedVector d_solution = factory.CreateDistributedVector(solution);

        Vec residual = factory.CreateVec();
        DistributedVector d_residual = factory.CreateDistributedVector(residual);


        for (DistributedVector::Iterator index = d_solution.Begin();
             index != d_solution.End();
             ++index)
        {
            d_solution[index]=index.Global;
            d_residual[index]=SIZE+index.Global;
        }

        d_solution.Restore();
        d_residual.Restore();

        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,1> bcc3;

        for (int i = 0; i < SIZE-1; i++)
        {
            nodes_array[i] = new Node<3>(i, true);
            ConstBoundaryCondition<3>* p_boundary_condition =
                new ConstBoundaryCondition<3>(-1);
            bcc3.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition);
        }

        bcc3.ApplyDirichletToNonlinearResidual(solution, residual, factory);


        for (DistributedVector::Iterator index = d_solution.Begin();
             index != d_solution.End();
             ++index)
        {
            if (index.Global < SIZE-1)
            {
                TS_ASSERT_DELTA(d_solution[index], index.Global,   1e-12);
                TS_ASSERT_DELTA(d_residual[index], index.Global+1, 1e-12);
            }
            else
            {
                TS_ASSERT_DELTA(d_solution[index], 9,   1e-12);
                TS_ASSERT_DELTA(d_residual[index], 19,  1e-12);
            }
        }

        for (int i=0; i < SIZE-1; i++)
        {
            delete nodes_array[i];
        }

        PetscTools::Destroy(solution);
        PetscTools::Destroy(residual);
    }

    void TestDefineZeroDirichletOnMeshBoundary()
    {
        // Load a 2D square mesh with 1 central non-boundary node
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        BoundaryConditionsContainer<2,2,1> bcc;

        bcc.DefineZeroDirichletOnMeshBoundary(&mesh);

        // Check boundary nodes have the right condition
        for (int i=0; i<4; i++)
        {
            double value = bcc.GetDirichletBCValue(mesh.GetNode(i));
            TS_ASSERT_DELTA(value, 0.0, 1e-12);
        }

        // Check non-boundary node has no condition
        TS_ASSERT(!bcc.HasDirichletBoundaryCondition(mesh.GetNode(4)));
    }

    void TestAnyNonZeroNeumannConditionsAndApplyNeumannToMeshBoundary()
    {
        // Load a 2D square mesh with 1 central non-boundary node
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        BoundaryConditionsContainer<2,2,1> bcc;
        BoundaryConditionsContainer<2,2,2> bcc_2unknowns;

        TS_ASSERT_EQUALS(bcc.AnyNonZeroNeumannConditions(), false);

        bcc.DefineZeroNeumannOnMeshBoundary(&mesh);

        TetrahedralMesh<2,2>::BoundaryElementIterator iter;
        iter = mesh.GetBoundaryElementIteratorBegin();
        while (iter != mesh.GetBoundaryElementIteratorEnd())
        {
            TS_ASSERT(bcc.HasNeumannBoundaryCondition(*iter));
            double value = bcc.GetNeumannBCValue(*iter, (*iter)->GetNode(0)->GetPoint());
            TS_ASSERT_DELTA(value, 0.0, 1e-8);

            iter++;
        }
        TS_ASSERT_EQUALS(bcc.AnyNonZeroNeumannConditions(), false);

        iter = mesh.GetBoundaryElementIteratorBegin();

        ConstBoundaryCondition<2>* p_boundary_condition2 = new ConstBoundaryCondition<2>(-1);

        bcc_2unknowns.AddNeumannBoundaryCondition(*iter, p_boundary_condition2);
        TS_ASSERT_EQUALS(bcc_2unknowns.AnyNonZeroNeumannConditions(), true);
    }

    void TestValidate()
    {
        // Load a 2D square mesh with 1 central non-boundary node
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        BoundaryConditionsContainer<2,2,1> bcc;

        // No BCs yet, so shouldn't validate
        TS_ASSERT(!bcc.Validate(&mesh));

        // Add some BCs
        ConstBoundaryCondition<2> *bc = new ConstBoundaryCondition<2>(0.0);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(0), bc);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(1), bc);
        bcc.AddDirichletBoundaryCondition(mesh.GetNode(3), bc);
        TetrahedralMesh<2,2>::BoundaryElementIterator iter
        = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, bc); // 2 to 3
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, bc); // 1 to 2

        TS_ASSERT(bcc.Validate(&mesh));
    }

    void TestAddNeumannBoundaryConditions()
    {
          // Load a 2D square mesh with 1 central non-boundary node
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        BoundaryConditionsContainer<2,2,2> bcc;

        // No BCs yet, so shouldn't validate
        TS_ASSERT(!bcc.Validate(&mesh));

        // Add some BCs
        ConstBoundaryCondition<2> *bc1 = new ConstBoundaryCondition<2>(2.0);
        ConstBoundaryCondition<2> *bc2 = new ConstBoundaryCondition<2>(-3.0);

        TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, bc1, 0);
        bcc.AddNeumannBoundaryCondition(*iter, bc2, 1);
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, bc1, 0);
        iter--;
        bcc.AddNeumannBoundaryCondition(*iter, bc2, 1);

        iter = mesh.GetBoundaryElementIteratorEnd();
        iter--;
        TS_ASSERT_DELTA(bcc.GetNeumannBCValue(*iter, ChastePoint<2>(), 0), 2.0, 1e-9);
        TS_ASSERT_DELTA(bcc.GetNeumannBCValue(*iter, ChastePoint<2>(), 1), -3.0, 1e-9);

        iter--;
        TS_ASSERT_DELTA(bcc.GetNeumannBCValue(*iter, ChastePoint<2>(), 0), 2.0, 1e-9);
        TS_ASSERT_DELTA(bcc.GetNeumannBCValue(*iter, ChastePoint<2>(), 1), 0.0, 1e-9);

        iter--;
        TS_ASSERT_DELTA(bcc.GetNeumannBCValue(*iter, ChastePoint<2>(), 0), 0.0, 1e-9);
        TS_ASSERT_DELTA(bcc.GetNeumannBCValue(*iter, ChastePoint<2>(), 1), -3.0, 1e-9);
    }

    void TestApplyToLinearSystem2Unknowns()
    {
        const int SIZE = 10;

        DistributedVectorFactory factory(SIZE);
        Vec template_vec = factory.CreateVec(2);
        LinearSystem some_system(template_vec, 2*SIZE);
        PetscTools::Destroy(template_vec);

        //LinearSystem some_system(2*SIZE);
        for (int i = 0; i < 2*SIZE; i++)
        {
            for (int j = 0; j < 2*SIZE; j++)
            {
                // LHS matrix is all 1s
                some_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            some_system.SetRhsVectorElement(i,2);
        }

        some_system.AssembleIntermediateLinearSystem();

        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,2> bcc32;

        // Apply dirichlet boundary conditions to all but last node
        for (int i = 0; i < SIZE-1; i++)
        {
            nodes_array[i] = new Node<3>(i,true);

            ConstBoundaryCondition<3>* p_boundary_condition0 =
                new ConstBoundaryCondition<3>(-1);

            ConstBoundaryCondition<3>* p_boundary_condition1 =
                new ConstBoundaryCondition<3>(-2);

            bcc32.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition0, 0);
            bcc32.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition1, 1);
        }
        bcc32.ApplyDirichletToLinearProblem(some_system);

        some_system.AssembleFinalLinearSystem();

        /*
         * Matrix should now look like
         *   A = (1 0 0 0 .. 0)
         *       (0 1 0 0 .. 0)
         *       (     ..     )
         *       (1 1 ..     1)
         *       (1 1 ..     1)
         * and rhs vector looks like b=(-1, -2, -1, -2, ..., -1, -2, 2, 2)
         * so solution of Ax = b is  x=(-1, -2, -1, -2, ..., -1, -2, ?, ?).
         */
        Vec solution = some_system.Solve();
        DistributedVector d_solution = factory.CreateDistributedVector(solution);
        DistributedVector::Stripe solution0(d_solution,0);
        DistributedVector::Stripe solution1(d_solution,1);

        for (DistributedVector::Iterator index = d_solution.Begin();
             index != d_solution.End();
             ++index)
        {
            if (index.Global!=SIZE-1) // last element of each stripe is not tested -- see ? in previous comment
            {
                TS_ASSERT_DELTA(solution0[index], -1.0, 0.000001);
                TS_ASSERT_DELTA(solution1[index], -2.0, 0.000001);
            }

        }

        for (int i = 0; i < SIZE-1; i++)
        {
            delete nodes_array[i];
        }

        PetscTools::Destroy(solution);
    }

    void TestApplyToLinearSystem3Unknowns()
    {
        const int SIZE = 10;

        DistributedVectorFactory factory(SIZE);
        Vec template_vec = factory.CreateVec(3);
        LinearSystem some_system(template_vec, 3*SIZE);
        PetscTools::Destroy(template_vec);

        //LinearSystem some_system(3*SIZE);
        for (int i = 0; i < 3*SIZE; i++)
        {
            for (int j = 0; j < 3*SIZE; j++)
            {
                // LHS matrix is all 1s
                some_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            some_system.SetRhsVectorElement(i,2);
        }

        some_system.AssembleIntermediateLinearSystem();

        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,3> bcc33;

        // Apply dirichlet boundary conditions to all but last node
        for (int i = 0; i < SIZE-1; i++)
        {
            nodes_array[i] = new Node<3>(i,true);

            ConstBoundaryCondition<3>* p_boundary_condition0 =
                new ConstBoundaryCondition<3>(-1);

            ConstBoundaryCondition<3>* p_boundary_condition1 =
                new ConstBoundaryCondition<3>(-2);

            ConstBoundaryCondition<3>* p_boundary_condition2 =
                new ConstBoundaryCondition<3>( 0);

            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition0, 0);
            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition1, 1);
            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition2, 2);
        }
        bcc33.ApplyDirichletToLinearProblem(some_system);

        some_system.AssembleFinalLinearSystem();

        Vec solution = some_system.Solve();

        DistributedVector d_solution = factory.CreateDistributedVector(solution);
        DistributedVector::Stripe solution0(d_solution,0);
        DistributedVector::Stripe solution1(d_solution,1);
        DistributedVector::Stripe solution2(d_solution,2);

        for (DistributedVector::Iterator index = d_solution.Begin();
             index != d_solution.End();
             ++index)
        {
            if (index.Global!=SIZE-1)
            {
                TS_ASSERT_DELTA(solution0[index], -1.0, 0.000001);
                TS_ASSERT_DELTA(solution1[index], -2.0, 0.000001);
                TS_ASSERT_DELTA(solution2[index],  0.0, 0.000001);
            }

        }
        for (int i = 0; i < SIZE-1; i++)
        {
            delete nodes_array[i];
        }
        PetscTools::Destroy(solution);
    }

    void TestApplyToNonlinearSystem3Unknowns()
    {
        const int SIZE = 10;
        DistributedVectorFactory factory(SIZE);

        Vec solution = factory.CreateVec(3);
        Vec residual = factory.CreateVec(3);

        double* p_solution;
        VecGetArray(solution, &p_solution);

        double* p_residual;
        VecGetArray(residual, &p_residual);

        int lo, hi;
        VecGetOwnershipRange(solution, &lo, &hi);

        for (int global_index=lo; global_index<hi; global_index++)
        {
            int local_index = global_index - lo;
            p_solution[local_index] = global_index;
            p_residual[local_index] = 100;
        }

        VecRestoreArray(solution, &p_solution);
        VecRestoreArray(residual, &p_residual);

        Node<3>* nodes_array[SIZE];
        BoundaryConditionsContainer<3,3,3> bcc33;

        for (int i = 0; i < SIZE; i++)
        {
            nodes_array[i] = new Node<3>(i,true);

            ConstBoundaryCondition<3>* p_boundary_condition0 = new ConstBoundaryCondition<3>(-1);
            ConstBoundaryCondition<3>* p_boundary_condition1 = new ConstBoundaryCondition<3>(-2);
            ConstBoundaryCondition<3>* p_boundary_condition2 = new ConstBoundaryCondition<3>(-3);

            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition0, 0);
            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition1, 1);
            bcc33.AddDirichletBoundaryCondition(nodes_array[i], p_boundary_condition2, 2);
        }

        bcc33.ApplyDirichletToNonlinearResidual(solution, residual, factory);

        VecGetArray(solution, &p_solution);
        VecGetArray(residual, &p_residual);

        for (int global_index = lo; global_index < hi; global_index++)
        {
            int local_index = global_index - lo;

            if (global_index%3==0)
            {
                TS_ASSERT_DELTA(p_solution[local_index], global_index,   1e-12);
                TS_ASSERT_DELTA(p_residual[local_index], global_index+1, 1e-12);
            }
            if (global_index%3==1)
            {
                TS_ASSERT_DELTA(p_solution[local_index], global_index, 1e-12);
                TS_ASSERT_DELTA(p_residual[local_index], global_index+2, 1e-12);
            }
            if (global_index%3==2)
            {
                TS_ASSERT_DELTA(p_solution[local_index], global_index,   1e-12);
                TS_ASSERT_DELTA(p_residual[local_index], global_index+3, 1e-12);
            }
        }
        for (int i = 0; i < SIZE; i++)
        {
            delete nodes_array[i];
        }

        VecRestoreArray(solution, &p_solution);
        VecRestoreArray(residual, &p_residual);

        PetscTools::Destroy(solution);
        PetscTools::Destroy(residual);
    }

    void TestArchiving()
    {
        OutputFileHandler handler("bcc_archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("bcc.arch");

        // Load a 2D square mesh with 1 central non-boundary node
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        {
            BoundaryConditionsContainer<2,2,2>* p_bcc = new BoundaryConditionsContainer<2,2,2>;

            // No BCs yet, so shouldn't validate
            TS_ASSERT_EQUALS(p_bcc->Validate(&mesh), false);

            // Add some Neumann BCs
            ConstBoundaryCondition<2> *bc1 = new ConstBoundaryCondition<2>(2.0);
            ConstBoundaryCondition<2> *bc2 = new ConstBoundaryCondition<2>(-3.0);

            TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
            iter--;
            p_bcc->AddNeumannBoundaryCondition(*iter, bc1, 0);
            p_bcc->AddNeumannBoundaryCondition(*iter, bc2, 1);
            iter--;
            p_bcc->AddNeumannBoundaryCondition(*iter, bc1, 0);
            iter--;
            p_bcc->AddNeumannBoundaryCondition(*iter, bc2, 1);

            // Add some Dirichlet BCs
            ConstBoundaryCondition<2> *bc3 = new ConstBoundaryCondition<2>(2.0);
            p_bcc->AddDirichletBoundaryCondition(mesh.GetNode(2), bc3);
            p_bcc->AddDirichletBoundaryCondition(mesh.GetNode(3), bc3);

            // Archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch & p_bcc;
            delete p_bcc;
        }

        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Load container...
            BoundaryConditionsContainer<2,2,2>* p_bcc;
            input_arch >> p_bcc;
            // ...and fill it
            p_bcc->LoadFromArchive( input_arch, &mesh );

            TetrahedralMesh<2,2>::BoundaryElementIterator iter = mesh.GetBoundaryElementIteratorEnd();
            iter--;
            TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter, ChastePoint<2>(), 0), 2.0, 1e-9);
            TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter, ChastePoint<2>(), 1), -3.0, 1e-9);

            iter--;
            TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter, ChastePoint<2>(), 0), 2.0, 1e-9);
            TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter, ChastePoint<2>(), 1), 0.0, 1e-9);

            iter--;
            TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter, ChastePoint<2>(), 0), 0.0, 1e-9);
            TS_ASSERT_DELTA(p_bcc->GetNeumannBCValue(*iter, ChastePoint<2>(), 1), -3.0, 1e-9);

            TS_ASSERT_DELTA(p_bcc->GetDirichletBCValue(mesh.GetNode(2)), 2.0, 1.e-6);
            TS_ASSERT_DELTA(p_bcc->GetDirichletBCValue(mesh.GetNode(3)), 2.0, 1.e-6);

            delete p_bcc;
        }
    }
};

#endif //_TESTBOUNDARYCONDITIONCONTAINER_HPP_
