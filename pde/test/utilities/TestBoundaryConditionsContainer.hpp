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

#ifndef _TESTBOUNDARYCONDITIONCONTAINER_HPP_
#define _TESTBOUNDARYCONDITIONCONTAINER_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include <boost/shared_ptr.hpp>
#include <vector>

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

        {
            std::vector<Node<1>*> nodes(num_nodes);

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
        }
        int num_elem = 10;
        std::vector<BoundaryElement<0,1> > elements;
        for (unsigned element_index=0; element_index< (unsigned) num_elem; element_index++)
        {
            std::vector<Node<1>*> nodes;
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

        std::vector<Node<2>*> nodes2(num_nodes);
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

        std::vector<Node<3>*> nodes3(num_nodes);
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
        LinearSystem linear_system(SIZE, SIZE);
        for (int i=0; i<SIZE; i++)
        {
            for (int j=0; j<SIZE; j++)
            {
                // LHS matrix is all 1s
                linear_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            linear_system.SetRhsVectorElement(i,2);
        }

        linear_system.AssembleIntermediateLinearSystem();

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
        bcc3.ApplyDirichletToLinearProblem(linear_system, true, false);
        ReplicatableVector vec_repl(linear_system.GetRhsVector());
        for (unsigned i=0; i<(unsigned)SIZE; i++)
        {
            TS_ASSERT_EQUALS(vec_repl[i], 2.0);
        }

        // now apply to the rhs vector
        bcc3.ApplyDirichletToLinearProblem(linear_system, false, true);

        linear_system.AssembleFinalLinearSystem();

        //////////////////////////
        // 2007 AD code from here
        //////////////////////////

        /*
         *  Based on the original system and the boundary conditions applied in a non-symmetric
         *  manner, the resulting linear system looks like:
         *
         *      1 0 0 ... 0
         *      0 1 0 ... 0
         *      0 0 1 ... 0
         *      ...
         *      1 1 1 ... 1
         *
         */
        /// \todo: this is very naughty. Must be checked in parallel as well.
        PetscInt lo, hi;
        PetscMatTools::GetOwnershipRange(linear_system.rGetLhsMatrix(), lo, hi);
        for (int row=lo; row<hi; row++)
        {
            if (row<SIZE-1)
            {
                for (int column=0; column<row; column++)
                {
                    TS_ASSERT_EQUALS(linear_system.GetMatrixElement(row,column), 0);
                }

                TS_ASSERT_EQUALS(linear_system.GetMatrixElement(row,row), 1);

                for (int column=row+1; column<SIZE; column++)
                {
                    TS_ASSERT_EQUALS(linear_system.GetMatrixElement(row,column), 0);
                }
            }

            if (row==SIZE-1)
            {
                for (int column=0; column<SIZE; column++)
                {
                    TS_ASSERT_EQUALS(linear_system.GetMatrixElement(row,column), 1);
                }
            }
        }

        Vec solution = linear_system.Solve();

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
        LinearSystem linear_system(SIZE, SIZE);
        linear_system.SetMatrixIsSymmetric(true);

        for (int i=0; i<SIZE; i++)
        {
            for (int j=0; j<SIZE; j++)
            {
                // LHS matrix is all 1s
                linear_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            linear_system.SetRhsVectorElement(i,2);
        }

        linear_system.AssembleIntermediateLinearSystem();

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
        bcc3.ApplyDirichletToLinearProblem(linear_system);

        linear_system.AssembleFinalLinearSystem();

        /*
         * Based on the original system and the boundary conditions applied in a symmetric
         * manner, the resulting linear system looks like:
         *
         *      1 0 0 ... 0
         *      0 1 0 ... 0
         *      0 0 1 ... 0
         *      ...
         *      0 0 0 ... 1
         */
        int lo, hi;
        linear_system.GetOwnershipRange(lo, hi);

        for (int row=lo; row<hi; row++)
        {
            for (int column=0; column<row; column++)
            {
                TS_ASSERT_EQUALS(linear_system.GetMatrixElement(row,column), 0);
            }

            TS_ASSERT_EQUALS(linear_system.GetMatrixElement(row,row), 1);

            for (int column=row+1; column<SIZE; column++)
            {
                TS_ASSERT_EQUALS(linear_system.GetMatrixElement(row,column), 0);
            }
        }

        Vec solution = linear_system.Solve();

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
        LinearSystem linear_system(template_vec, 2*SIZE);
        PetscTools::Destroy(template_vec);

        for (int i = 0; i < 2*SIZE; i++)
        {
            for (int j = 0; j < 2*SIZE; j++)
            {
                // LHS matrix is all 1s
                linear_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            linear_system.SetRhsVectorElement(i,2);
        }

        linear_system.AssembleIntermediateLinearSystem();

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
        bcc32.ApplyDirichletToLinearProblem(linear_system);

        linear_system.AssembleFinalLinearSystem();

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
        Vec solution = linear_system.Solve();
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
        LinearSystem linear_system(template_vec, 3*SIZE);
        PetscTools::Destroy(template_vec);

        for (int i = 0; i < 3*SIZE; i++)
        {
            for (int j = 0; j < 3*SIZE; j++)
            {
                // LHS matrix is all 1s
                linear_system.SetMatrixElement(i,j,1);
            }
            // RHS vector is all 2s
            linear_system.SetRhsVectorElement(i,2);
        }

        linear_system.AssembleIntermediateLinearSystem();

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
        bcc33.ApplyDirichletToLinearProblem(linear_system);

        linear_system.AssembleFinalLinearSystem();

        Vec solution = linear_system.Solve();

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

    void TestPerdiodicBoundaryConditions()
    {
        const int SIZE = 5;

        DistributedVectorFactory factory(SIZE);
        Vec template_vec = factory.CreateVec(2);
        LinearSystem linear_system(template_vec, 2*SIZE);
        PetscTools::Destroy(template_vec);

        for (int i = 0; i < 2*SIZE; i++)
        {
            for (int j = 0; j < 2*SIZE; j++)
            {
                // LHS matrix is all 2s
                linear_system.SetMatrixElement(i,j,2);
            }
            // RHS vector is all 3s
            linear_system.SetRhsVectorElement(i,3);
        }

        linear_system.AssembleIntermediateLinearSystem();

        Node<3>* nodes[SIZE];
        BoundaryConditionsContainer<3,3,2> bcc;

        for (unsigned i=0; i<SIZE-1; i++)
        {
            nodes[i] = new Node<3>(i,true);
        }

        bcc.AddPeriodicBoundaryCondition(nodes[0], nodes[1]);
        bcc.AddPeriodicBoundaryCondition(nodes[2], nodes[3]);

        bcc.ApplyPeriodicBcsToLinearProblem(linear_system, true, true);

        linear_system.AssembleFinalLinearSystem();

        ReplicatableVector rhs_repl(linear_system.rGetRhsVector());
        TS_ASSERT_DELTA(rhs_repl[0], 0.0, 1e-12); // node 0, variable 0
        TS_ASSERT_DELTA(rhs_repl[1], 0.0, 1e-12); // node 0, variable 1
        TS_ASSERT_DELTA(rhs_repl[2], 3.0, 1e-12);
        TS_ASSERT_DELTA(rhs_repl[3], 3.0, 1e-12);
        TS_ASSERT_DELTA(rhs_repl[4], 0.0, 1e-12); // node 2, variable 0
        TS_ASSERT_DELTA(rhs_repl[5], 0.0, 1e-12); // node 2, variable 1
        TS_ASSERT_DELTA(rhs_repl[6], 3.0, 1e-12);
        TS_ASSERT_DELTA(rhs_repl[7], 3.0, 1e-12);
        TS_ASSERT_DELTA(rhs_repl[8], 3.0, 1e-12);
        TS_ASSERT_DELTA(rhs_repl[9], 3.0, 1e-12);


        //
        //  Matrix should have
        //   row 0 altered to be [1, 0 -1, 0, 0, ..., 0]
        //   row 1 altered to be [0, 1, 0 -1, 0, ..., 0]
        //   row 4 altered to be [0, 0, 0, 0, 1, 0, -1, 0, 0, 0]
        //   row 5 altered to be [0, 0, 0, 0, 0, 1, 0, -1, 0, 0]
        //   All other rows just [2, 2, ..., 2]

        Mat& r_mat = linear_system.rGetLhsMatrix();

        PetscInt lo, hi;
        PetscMatTools::GetOwnershipRange(r_mat, lo, hi);
        for (int i=lo; i<hi; i++)
        {
            if (i==0 || i==1 || i==4 || i==5)
            {
                unsigned col_one = i;
                unsigned col_minus_one = i+2;

                for (unsigned j=0; j<2*SIZE; j++)
                {
                    double val = 0.0;
                    if (j==col_one)
                    {
                        val = 1.0;
                    }
                    if (j==col_minus_one)
                    {
                        val = -1.0;
                    }
                    TS_ASSERT_DELTA( PetscMatTools::GetElement(r_mat, i, j), val, 1e-12);
                }
            }
            else
            {
                for (unsigned j=0; j<2*SIZE; j++)
                {
                    TS_ASSERT_DELTA( PetscMatTools::GetElement(r_mat, i, j), 2.0, 1e-12);
                }
            }
        }

        for (unsigned i=0; i<SIZE-1; i++)
        {
            delete nodes[i];
        }
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
