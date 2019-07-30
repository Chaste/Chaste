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

#ifndef TESTFINECOARSEMESHPAIR_HPP_
#define TESTFINECOARSEMESHPAIR_HPP_

#include <cxxtest/TestSuite.h>
#include "FineCoarseMeshPair.hpp"
#include "TetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestFineCoarseMeshPair : public CxxTest::TestSuite
{
public:

    // Simple test where the whole of the coarse mesh is in one fine element
    void TestComputeFineElemsAndWeightsForQuadPointsSimple()
    {
        TetrahedralMesh<2,2> fine_mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        fine_mesh.ConstructFromMeshReader(mesh_reader);

        QuadraticMesh<2> coarse_mesh(0.1, 0.1, 0.1);
        coarse_mesh.Translate(0.5, 0.0); // whole of the coarse mesh in now in fine element with index 1

        FineCoarseMeshPair<2> mesh_pair(fine_mesh, coarse_mesh);

        //test get methods
        TS_ASSERT_EQUALS(mesh_pair.GetFineMesh().GetNumAllElements(), 4u);
        TS_ASSERT_EQUALS(mesh_pair.GetCoarseMesh().GetNumAllElements(), 2u);

        mesh_pair.SetUpBoxesOnFineMesh();
        GaussianQuadratureRule<2> quad_rule(3);
        unsigned num_quads_per_element = quad_rule.GetNumQuadPoints();
        TS_ASSERT_EQUALS(num_quads_per_element, 6u);
        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true);

        // All coarse quadrature points should have been found in the fine mesh
        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 0u);
        TS_ASSERT_EQUALS(mesh_pair.mNotInMeshNearestElementWeights.size(), 0u);

        // Check the elements and weights have been set up correctly
        // Each item corresponds to a quad point in an element in the coarse mesh with 6 quad points per element
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), num_quads_per_element*coarse_mesh.GetNumAllElements());
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), 12u);

        for (unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            // All coarse mesh quad points are in the same fine element
            TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[i].ElementNum, 1u);

            /*
             * All the weights should be between 0 and 1 as no coarse nodes are.
             * Note weights = (1-psi_x-psi_y, psi_x, psi_y), where psi is the
             * position of the point in that element when transformed to the
             * canonical element.
             */
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_LESS_THAN(-1e14, mesh_pair.rGetElementsAndWeights()[i].Weights(j));
                TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].Weights(j), 1.0+1e-14);
            }
        }
        TS_ASSERT_EQUALS(mesh_pair.mStatisticsCounters[0], 12u);
        TS_ASSERT_EQUALS(mesh_pair.mStatisticsCounters[1], 0u);
    }

    void TestWithCoarseContainedInFine()
    {
        // Fine mesh is has h=0.1, on unit cube (so 6000 elements)
        TetrahedralMesh<3,3> fine_mesh;
        fine_mesh.ConstructRegularSlabMesh(0.1, 1.0, 1.0, 1.0);

        // Coarse mesh is has h=1 on unit cube (so 6 elements)
        QuadraticMesh<3> coarse_mesh(1.0, 1.0, 1.0, 1.0);

        FineCoarseMeshPair<3> mesh_pair(fine_mesh,coarse_mesh);

        mesh_pair.SetUpBoxesOnFineMesh(0.3);

        TS_ASSERT_EQUALS(mesh_pair.mpFineMeshBoxCollection->GetNumBoxes(), 4*4*4u);

        // For each node, find containing box. That box should contain any element that node is in.
        for (unsigned i=0; i<fine_mesh.GetNumNodes(); i++)
        {
            unsigned box_index = mesh_pair.mpFineMeshBoxCollection->CalculateContainingBox(fine_mesh.GetNode(i));
            if (mesh_pair.mpFineMeshBoxCollection->IsBoxOwned(box_index))
            {
                assert(fine_mesh.GetNode(i)->rGetContainingElementIndices().size() > 0);

                for (std::set<unsigned>::iterator iter = fine_mesh.GetNode(i)->rGetContainingElementIndices().begin();
                        iter != fine_mesh.GetNode(i)->rGetContainingElementIndices().end();
                        ++iter)
                {
                    Element<3,3>* p_element = fine_mesh.GetElement(*iter);
                    TS_ASSERT_DIFFERS( mesh_pair.mpFineMeshBoxCollection->rGetBox(box_index).rGetElementsContained().find(p_element), mesh_pair.mpFineMeshBoxCollection->rGetBox(box_index).rGetElementsContained().end() )
                }
            }
        }

        GaussianQuadratureRule<3> quad_rule(3);
        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true);

        // All coarse quadrature points should have been found in the fine mesh
        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 0u);
        TS_ASSERT_EQUALS(mesh_pair.mNotInMeshNearestElementWeights.size(), 0u);

        // Check the elements and weights have been set up correctly
        // 6 elements, 8 quad points
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), 6*8u);

        // Some hardcoded values, just to check element_nums not all zero
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[0].ElementNum,  3816u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[10].ElementNum, 217u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[20].ElementNum, 1094u);

        for (unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].ElementNum, fine_mesh.GetNumElements());

            /*
             * As all the quadrature points should have been found in the fine mesh,
             * all the weights should be between 0 and 1.
             * Note weights = (1-psi_x-psi_y-psi_z, psi_x, psi_y, psi_z), where psi
             * is the position of the point in that element when transformed to the
             * canonical element.
             */
            for (unsigned j=0; j<4; j++)
            {
                TS_ASSERT_LESS_THAN(-1e14, mesh_pair.rGetElementsAndWeights()[i].Weights(j));
                TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].Weights(j), 1.0+1e-14);
            }
        }

        TS_ASSERT_EQUALS(mesh_pair.mStatisticsCounters[0], 6*8u);
        TS_ASSERT_EQUALS(mesh_pair.mStatisticsCounters[1], 0u);
        mesh_pair.PrintStatistics();

        mesh_pair.DeleteFineBoxCollection();
        TS_ASSERT(mesh_pair.mpFineMeshBoxCollection==NULL);
    }

    void TestWithCoarseSlightlyOutsideFine()
    {
        // Fine mesh is has h=0.1, on unit cube (so 6000 elements)
        TetrahedralMesh<3,3> fine_mesh;
        fine_mesh.ConstructRegularSlabMesh(0.1, 1.0,1.0,1.0);

        // Coarse mesh is slightly bigger than in previous test
        QuadraticMesh<3> coarse_mesh(1.0, 1.0, 1.0, 1.0); // xmax > 1.0
        coarse_mesh.Scale(1.03, 1.0, 1.0);

        FineCoarseMeshPair<3> mesh_pair(fine_mesh,coarse_mesh);

        GaussianQuadratureRule<3> quad_rule(3);

        // Need to call SetUpBoxesOnFineMesh first
        TS_ASSERT_THROWS_CONTAINS(mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true), "Call");

        mesh_pair.SetUpBoxesOnFineMesh(0.3);
        TS_ASSERT_EQUALS(mesh_pair.mpFineMeshBoxCollection->GetNumBoxes(), 4*4*4u);

        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true);


        ///\todo #2308 These quantities are not shared yet...
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 2u); // hardcoded
            TS_ASSERT_EQUALS(mesh_pair.mNotInMeshNearestElementWeights.size(), 2u);
        }
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), 6*8u);

        for (unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].ElementNum, fine_mesh.GetNumElements());

            // comment out this test as now some of the weights are negative/greater than one
            //for (unsigned j=0; j<4; j++)
            //{
            //    TS_ASSERT_LESS_THAN(-1e14, mesh_pair.rGetElementsAndWeights()[i].Weights(j));
            //    TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].Weights(j), 1.0+1e-14);
            //}
        }

        /*
         * For each quadrature point that was not found in the fine mesh, check
         * that its x-value is greater than one - this is the only way it could
         * be outside the fine mesh.
         */
        QuadraturePointsGroup<3> quad_point_posns(coarse_mesh, quad_rule);
        for (unsigned i=0; i<mesh_pair.mNotInMesh.size(); i++)
        {
            double x = quad_point_posns.rGet(mesh_pair.mNotInMesh[i])(0);
            TS_ASSERT_LESS_THAN(1.0, x);
        }

        mesh_pair.PrintStatistics();

        ///\todo #2308 warnings may or may not be triggered locally
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS( Warnings::Instance()->GetNumWarnings(), 1u);
        }
        Warnings::Instance()->QuietDestroy();
    }

////Bring back this functionality if needed
//    void dontTestWithIdenticalMeshes()
//    {
//        TrianglesMeshReader<1,1> reader1("mesh/test/data/1D_0_to_1_10_elements");
//        TetrahedralMesh<1,1> fine_mesh;
//        fine_mesh.ConstructFromMeshReader(reader1);
//
//        TrianglesMeshReader<1,1> reader2("mesh/test/data/1D_0_to_1_10_elements_quadratic",2);
//        QuadraticMesh<1> coarse_mesh;
//        coarse_mesh.ConstructFromMeshReader(reader2);
//
//        FineCoarseMeshPair<1> mesh_pair(fine_mesh,coarse_mesh);
//        TS_ASSERT_EQUALS(mesh_pair.mIdenticalMeshes, true);
//
//        GaussianQuadratureRule<1> quad_rule(0);
//        mesh_pair.SetUpBoxesOnFineMesh(0.3);
//
//        // Covers the mIdenticalMeshes=true part of this method. Would throw exception if can't find
//        // quad point in first choice of element.
//        TS_ASSERT_THROWS_NOTHING(mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true));
//    }

    void TestWithDefaultBoxWidth()
    {
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRegularSlabMesh(0.1, 1.0, 1.0);

        QuadraticMesh<2> coarse_mesh(1.0, 1.0, 1.0);

        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);

        mesh_pair.SetUpBoxesOnFineMesh();

        /*
         * With this mesh the proposed width - the width that would correspond to 20 boxes
         * in the x-direction is:
         *   Proposed width = 0.0552632
         * but
         *   max_edge_length = 0.141421
         * and we want width > max_edge_length, so end up with
         *   box width = 0.155563
         * (1.1 times max edge length)
         */
        TS_ASSERT_EQUALS(mesh_pair.mpFineMeshBoxCollection->GetNumBoxes(), 64u);

        // Now use a mesh with a smaller edge length
        TetrahedralMesh<2,2> fine_mesh2;
        fine_mesh2.ConstructRegularSlabMesh(0.01, 1.0,1.0);

        /*
         * Can use smaller boxes
         *  Proposed width = 0.0552632
         *  max_edge_length = 0.0141421
         *  box width = 0.0552632
         */
        FineCoarseMeshPair<2> mesh_pair2(fine_mesh2,coarse_mesh);
        mesh_pair2.SetUpBoxesOnFineMesh();
        TS_ASSERT_EQUALS(mesh_pair2.mpFineMeshBoxCollection->GetNumBoxes(), 20*20u);
    }

    /*
     * Test when calling ComputeFineElementsAndWeightsForCoarseQuadPoints()
     * in non-safe mode, but using the default value of box width. It is
     * difficult to get the class to run incorrectly (ie fail without an
     * assertion failing) in non-safe mode (ie we can't just specify boxes
     * that are too small), so we just test we get the same results as in
     * safe mode.
     */
    void TestNonSafeMode()
    {
        // Fine mesh is has h=0.1, on unit cube (so 6000 elements)
        TetrahedralMesh<3,3> fine_mesh;
        fine_mesh.ConstructRegularSlabMesh(0.1, 1.0, 1.0, 1.0);

        QuadraticMesh<3> coarse_mesh(1.0, 1.0, 1.0, 1.0); // xmax > 1.0
        coarse_mesh.Scale(1.03, 1.0, 1.0);

        FineCoarseMeshPair<3> mesh_pair(fine_mesh,coarse_mesh);
        GaussianQuadratureRule<3> quad_rule(3);

        // Call SetUpBoxesOnFineMesh() without providing a width
        mesh_pair.SetUpBoxesOnFineMesh();

        /*
         * Whereas before 4 by 4 by 4 boxes were explicitly chosen, here it
         * has been determined that 6 by 6 by 6 are needed.
         */
        TS_ASSERT_EQUALS(mesh_pair.mpFineMeshBoxCollection->GetNumBoxes(), 6*6*6u);

        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, false /* non-safe mode*/);

        ///\todo #2308 These quantities are not shared yet...
        if (PetscTools::IsSequential())
        {
            TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 2u); // hardcoded
            TS_ASSERT_EQUALS(mesh_pair.mNotInMeshNearestElementWeights.size(), 2u);
        }
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), 6*8u);

        for (unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].ElementNum, fine_mesh.GetNumElements());
        }

        /*
         * For each quadrature point that was not found in the fine mesh, check
         * that its x-value is greater than one - this is the only way it could
         * be outside the fine mesh.
         */
        QuadraturePointsGroup<3> quad_point_posns(coarse_mesh, quad_rule);
        for (unsigned i=0; i<mesh_pair.mNotInMesh.size(); i++)
        {
            double x = quad_point_posns.rGet(mesh_pair.mNotInMesh[i])(0);
            TS_ASSERT_LESS_THAN(1.0, x);
        }

        mesh_pair.PrintStatistics();

        TS_ASSERT_EQUALS( Warnings::Instance()->GetNumWarnings(), 1u);
        Warnings::Instance()->QuietDestroy();
    }

    // Covers some bits that aren't covered in the tests above,
    void TestOtherCoverage()
    {
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRegularSlabMesh(0.1, 1.0, 1.0);

        QuadraticMesh<2> coarse_mesh(1.0, 1.0, 1.0);

        /*
         * Rotate the mesh by 45 degrees, makes it possible (since boxes no
         * longer lined up with elements) for the containing element of a
         * quad point to be in a *local* box, ie not an element contained
         * in the box containing this point.
         */
        c_matrix<double,2,2> rotation_mat;
        rotation_mat(0,0) = 1.0/sqrt(2.0);
        rotation_mat(1,0) = -1.0/sqrt(2.0);
        rotation_mat(0,1) = 1.0/sqrt(2.0);
        rotation_mat(1,1) = 1.0/sqrt(2.0);

        fine_mesh.Rotate(rotation_mat);
        coarse_mesh.Rotate(rotation_mat);

        GaussianQuadratureRule<2> quad_rule(3);

        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);
        mesh_pair.SetUpBoxesOnFineMesh();
        mesh_pair.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true);
        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 0u);

        /*
         * Repeat again with smaller boxes, covers the bit requiring the whole
         * mesh to be searched to find an element for a particular quad point.
         */
        FineCoarseMeshPair<2> mesh_pair2(fine_mesh,coarse_mesh);
        mesh_pair2.SetUpBoxesOnFineMesh(0.01);
        mesh_pair2.ComputeFineElementsAndWeightsForCoarseQuadPoints(quad_rule, true);
        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 0u);
    }

    void TestComputeCoarseElementsForFineNodes()
    {
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRegularSlabMesh(0.2, 1.0, 1.0);

        QuadraticMesh<2> coarse_mesh(1.0, 1.0, 1.0); // 2 triangular elements

        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);
        TS_ASSERT_THROWS_CONTAINS(mesh_pair.ComputeCoarseElementsForFineNodes(true),"Call SetUpBoxesOnCoarseMesh()");

        mesh_pair.SetUpBoxesOnCoarseMesh();
        mesh_pair.ComputeCoarseElementsForFineNodes(true);

        //Check that the indices of the coarse mesh elements are as expected
        unsigned lower_left_element_index=1u;
        ChastePoint<2> lower_left(0.25, 0.25);
        TS_ASSERT_EQUALS(coarse_mesh.GetContainingElementIndex(lower_left), lower_left_element_index);
        ChastePoint<2> lower_left1(0.1, 0.25); //Double check that there is a `backslash`
        ChastePoint<2> lower_left2(0.25, 0.1);
        TS_ASSERT_EQUALS(coarse_mesh.GetContainingElementIndex(lower_left1), lower_left_element_index);
        TS_ASSERT_EQUALS(coarse_mesh.GetContainingElementIndex(lower_left2), lower_left_element_index);

        unsigned upper_right_element_index=0u;
        ChastePoint<2> upper_right(0.75, 0.75);
        TS_ASSERT_EQUALS(coarse_mesh.GetContainingElementIndex(upper_right), upper_right_element_index);

        for (unsigned i=0; i<fine_mesh.GetNumNodes(); i++)
        {
            double x = fine_mesh.GetNode(i)->rGetLocation()[0];
            double y = fine_mesh.GetNode(i)->rGetLocation()[1];

            if (x+y < 1.0 - 1e-5)  // x+y < 1
            {
                TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineNodes()[i], lower_left_element_index);
            }
            else if (x+y > 1.0 + 1e-5)  // x+y > 1
            {
                TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineNodes()[i], upper_right_element_index);
            }
            else // x=1-y, so in both elements, result could be either. However, it should find 0 first
            {
                //TS_ASSERT_LESS_THAN(mesh_pair.rGetCoarseElementsForFineNodes()[i], 2u);
                TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineNodes()[i], 0u);
            }
        }

        /*
         * Translate the fine mesh in the (-1, -1) direction --> all fine nodes
         * nearest to (not contained in) element 0. We have to make the fine mesh
         * tiny and then translate a small amount so that it is still in the box
         * collection for the coarse (normally the two meshes should overlap).
         */
        fine_mesh.Scale(1e-2, 1e-2);
        fine_mesh.Translate(-1.1e-2, -1.1e-2);
        mesh_pair.ComputeCoarseElementsForFineNodes(true);
        for (unsigned i=0; i<fine_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineNodes()[i], lower_left_element_index);
        }
        // A little reset
        TS_ASSERT_DIFFERS(mesh_pair.rGetCoarseElementsForFineNodes()[0], 0u);
        mesh_pair.rGetCoarseElementsForFineNodes()[0] = 189342958;
        // Call again with safeMode=false this time (same results, faster)
        mesh_pair.ComputeCoarseElementsForFineNodes(false);
        for (unsigned i=0; i<fine_mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineNodes()[i], lower_left_element_index);
        }

        // Coverage:
        //  call again
        mesh_pair.SetUpBoxesOnCoarseMesh();
        //  delete
        mesh_pair.DeleteCoarseBoxCollection();
    }

    void TestComputeCoarseElementsForFineElementCentroids()
    {
        TetrahedralMesh<2,2> fine_mesh;
        fine_mesh.ConstructRegularSlabMesh(0.2, 1.0, 1.0);

        QuadraticMesh<2> coarse_mesh(1.0, 1.0, 1.0); // 2 triangular elements

        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);

        TS_ASSERT_THROWS_CONTAINS(mesh_pair.ComputeCoarseElementsForFineElementCentroids(true),"Call SetUpBoxesOnCoarseMesh()");

        mesh_pair.SetUpBoxesOnCoarseMesh();
        mesh_pair.ComputeCoarseElementsForFineElementCentroids(true);

        //Check that the indices of the coarse mesh elements are as expected
        unsigned lower_left_element_index=1u;
        ChastePoint<2> lower_left(0.25, 0.25);
        TS_ASSERT_EQUALS(coarse_mesh.GetContainingElementIndex(lower_left), lower_left_element_index);

        unsigned upper_right_element_index=0u;
        ChastePoint<2> upper_right(0.75, 0.75);
        TS_ASSERT_EQUALS(coarse_mesh.GetContainingElementIndex(upper_right), upper_right_element_index);

        TS_ASSERT_EQUALS( mesh_pair.rGetCoarseElementsForFineElementCentroids().size(), fine_mesh.GetNumElements());
        for (unsigned i=0; i<fine_mesh.GetNumElements(); i++)
        {
            double x = fine_mesh.GetElement(i)->CalculateCentroid()(0);
            double y = fine_mesh.GetElement(i)->CalculateCentroid()(1);
            if (x+y < 1.0)
            {
                TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineElementCentroids()[i], lower_left_element_index);
            }
            else
            {
                TS_ASSERT_EQUALS(mesh_pair.rGetCoarseElementsForFineElementCentroids()[i], upper_right_element_index);
            }
        }

        // Coverage
        mesh_pair.DeleteCoarseBoxCollection();
        mesh_pair.SetUpBoxesOnCoarseMesh(0.8); // force a point to be found in a neighbouring box
        mesh_pair.ComputeCoarseElementsForFineElementCentroids(true);

        mesh_pair.DeleteCoarseBoxCollection();
        mesh_pair.SetUpBoxesOnCoarseMesh(0.1); // force a point to be found in a nonlocal box
        mesh_pair.ComputeCoarseElementsForFineElementCentroids(true);

        mesh_pair.DeleteCoarseBoxCollection();
        mesh_pair.SetUpBoxesOnCoarseMesh(); // back to default

        /*
         * Translate the fine mesh in the (-1, -1) direction --> all fine elements
         * nearest to (not contained in) element 0. We have to make the fine mesh
         * tiny and then translate a small amount so that it is still in the box
         * collection for the coarse (normally the two meshes should overlap).
         */
        fine_mesh.Scale(1e-2, 1e-2);
        fine_mesh.Translate(-1.1e-2, -1.1e-2);
        mesh_pair.ComputeCoarseElementsForFineElementCentroids(true);
        TS_ASSERT_EQUALS( mesh_pair.rGetCoarseElementsForFineElementCentroids().size(), fine_mesh.GetNumElements());
        for (unsigned i=0; i<fine_mesh.GetNumElements(); i++)
        {
            TS_ASSERT_EQUALS( mesh_pair.rGetCoarseElementsForFineElementCentroids()[i], lower_left_element_index);
        }
    }

    void TestComputeFineElemsAndWeightsForCoarseNodes()
    {
        TetrahedralMesh<2,2> fine_mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        fine_mesh.ConstructFromMeshReader(mesh_reader);

        QuadraticMesh<2> coarse_mesh(0.5, 0.5, 0.5);
        coarse_mesh.Translate(0.2,0.1);

        FineCoarseMeshPair<2> mesh_pair(fine_mesh,coarse_mesh);

        // Need to call SetUpBoxesOnFineMesh first
        TS_ASSERT_THROWS_CONTAINS(mesh_pair.ComputeFineElementsAndWeightsForCoarseNodes(true), "Call");

        mesh_pair.SetUpBoxesOnFineMesh();
        mesh_pair.ComputeFineElementsAndWeightsForCoarseNodes(true);

        // All coarse quadrature points should have been found in the fine mesh
        TS_ASSERT_EQUALS(mesh_pair.mNotInMesh.size(), 0u);
        TS_ASSERT_EQUALS(mesh_pair.mNotInMeshNearestElementWeights.size(), 0u);

        // Check the elements and weights have been set up correctly
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights().size(), 9u);

        // Check the first four nodes against what they should be
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[0].ElementNum, 1u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[1].ElementNum, 1u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[2].ElementNum, 0u);
        TS_ASSERT_EQUALS(mesh_pair.rGetElementsAndWeights()[3].ElementNum, 2u);

        for (unsigned i=0; i<mesh_pair.rGetElementsAndWeights().size(); i++)
        {
            TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].ElementNum, fine_mesh.GetNumElements());

            // All the weights should be between 0 and 1 as no coarse nodes are
            // Note weights = (1-psi_x-psi_y-psi_z, psi_x, psi_y, psi_z), where psi is the position of the
            // point in that element when transformed to the canonical element
            for (unsigned j=0; j<3; j++)
            {
                TS_ASSERT_LESS_THAN(-1e14, mesh_pair.rGetElementsAndWeights()[i].Weights(j));
                TS_ASSERT_LESS_THAN(mesh_pair.rGetElementsAndWeights()[i].Weights(j), 1.0+1e-14);
            }
        }

        TS_ASSERT_EQUALS(mesh_pair.mStatisticsCounters[0], 9u);
        TS_ASSERT_EQUALS(mesh_pair.mStatisticsCounters[1], 0u);
    }
};

#endif /*TESTFINECOARSEMESHPAIR_HPP_*/
