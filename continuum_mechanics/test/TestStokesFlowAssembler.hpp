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

#ifndef TESTSTOKESFLOWASSEMBLER_HPP_
#define TESTSTOKESFLOWASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "UblasCustomFunctions.hpp"
#include "StokesFlowAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraticMesh.hpp"
#include "TrianglesMeshReader.hpp"

c_vector<double,2> MyBodyForce(c_vector<double,2>& rX, double t)
{
    // check the point rX has been interpolated correctly. The mesh
    // is the canonical triangle translated by (0.5,0.8).
    assert(rX(0)>0.0 + 0.5);
    assert(rX(1)>0.0 + 0.8);
    assert(rX(0)+rX(1)<1.0 + 0.5 + 0.8);

    c_vector<double,2> body_force;
    body_force(0) = 10.0;
    body_force(1) = 20.0;
    return body_force;
}

class TestStokesFlowAssembler : public CxxTest::TestSuite
{
public:
    /*
     * Test that the matrix is calculated correctly on the cannonical triangle.
     * Tests against the analytical solution calculated by hand.
     */
    void TestAssembler()
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/canonical_triangle_quadratic", 2, 2, false);
        mesh.ConstructFromMeshReader(mesh_reader);

        double mu = 2.0;
        c_vector<double,2> body_force;
        double g1 = 1.34254;
        double g2 = 75.3422;
        body_force(0) = g1;
        body_force(1) = g2;

        StokesFlowProblemDefinition<2> problem_defn(mesh);
        problem_defn.SetViscosity(mu);
        problem_defn.SetBodyForce(body_force);

        StokesFlowAssembler<2> assembler(&mesh, &problem_defn);

        // The tests below test the assembler against hand-calculated variables for
        // an OLD weak form (corresponding to different boundary conditions), not the
        // current Stokes weak form. This factor converts the assembler to use the old
        // weak form. See documentation for this variable for more details.
        assembler.mScaleFactor = 0.0;

        Vec vec = PetscTools::CreateVec(18);
        Mat mat;
        PetscTools::SetupMat(mat, 18, 18, 18);

        assembler.SetVectorToAssemble(vec, true);
        assembler.SetMatrixToAssemble(mat, true);
        assembler.Assemble();
        PetscMatTools::Finalise(mat);

        double A[6][6] = {
                           {      1.0,  1.0/6.0,  1.0/6.0,      0.0, -2.0/3.0, -2.0/3.0},
                           {  1.0/6.0,  1.0/2.0,      0.0,      0.0,      0.0, -2.0/3.0},
                           {  1.0/6.0,      0.0,  1.0/2.0,      0.0, -2.0/3.0,      0.0},
                           {      0.0,      0.0,      0.0,  8.0/3.0, -4.0/3.0, -4.0/3.0},
                           { -2.0/3.0,      0.0, -2.0/3.0, -4.0/3.0,  8.0/3.0,      0.0},
                           { -2.0/3.0, -2.0/3.0,      0.0, -4.0/3.0,      0.0,  8.0/3.0}
                         };

        double Bx[6][3] = {
                            { -1.0/6.0,      0.0,      0.0},
                            {      0.0,  1.0/6.0,      0.0},
                            {      0.0,      0.0,      0.0},
                            {  1.0/6.0,  1.0/6.0,  1.0/3.0},
                            { -1.0/6.0, -1.0/6.0, -1.0/3.0},
                            {  1.0/6.0, -1.0/6.0,      0.0},
                         };

        double By[6][3] = {
                            { -1.0/6.0,      0.0,      0.0},
                            {      0.0,      0.0,      0.0},
                            {      0.0,      0.0,  1.0/6.0},
                            {  1.0/6.0,  1.0/3.0,  1.0/6.0},
                            {  1.0/6.0,      0.0, -1.0/6.0},
                            { -1.0/6.0, -1.0/3.0, -1.0/6.0},
                          };

        c_matrix<double,18,18> exact_A = zero_matrix<double>(18);

        // The diagonal 6x6 blocks
        for (unsigned i=0; i<6; i++)
        {
            for (unsigned j=0; j<6; j++)
            {
                exact_A(3*i,  3*j)   = mu*A[i][j];
                exact_A(3*i+1,3*j+1) = mu*A[i][j];
            }
        }


        // The 6x3 Blocks
        for (unsigned i=0; i<6; i++)
        {
            for (unsigned j=0; j<3; j++)
            {
                exact_A(3*i,3*j+2)   = -Bx[i][j];
                exact_A(3*i+1,3*j+2) = -By[i][j];
                //- as -Div(U)=0
                exact_A(3*j+2,3*i)   = -Bx[i][j];
                exact_A(3*j+2,3*i+1) = -By[i][j];
            }
        }

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<18; j++)
            {
                TS_ASSERT_DELTA(PetscMatTools::GetElement(mat,i,j), exact_A(i,j), 1e-9);
            }
        }
        ReplicatableVector vec_repl(vec);

        // The first 6 entries in the vector correspond to nodes 0, 1, 2, i.e. the vertices.
        // For these nodes, it can be shown that the integral of the corresponding
        // basis function is zero, i.e. \intgl_{canonical element} \phi_i dV = 0.0  for i=0,1,2, phi_i the
        // i-th QUADRATIC basis.
        for (unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA(vec_repl[3*i],   g1*0.0, 1e-8);
            TS_ASSERT_DELTA(vec_repl[3*i+1], g2*0.0, 1e-8);
        }

        // The next 6 entries in the vector correspond to nodes 3, 4, 5, i.e. the internal edges.
        // For these nodes, it can be shown that the integral of the corresponding
        // basis function is 1/6, i.e. \intgl_{canonical element} \phi_i dV = 1/6  for i=3,4,5, phi_i the
        // i-th QUADRATIC basis.
        for (unsigned i=3; i<6; i++)
        {
            TS_ASSERT_DELTA(vec_repl[3*i],   g1/6.0, 1e-8);
            TS_ASSERT_DELTA(vec_repl[3*i+1], g2/6.0, 1e-8);
        }

        // The pressure-block of the RHS vector should be zero.
        TS_ASSERT_DELTA(vec_repl[2], 0.0, 1e-9);
        TS_ASSERT_DELTA(vec_repl[5], 0.0, 1e-9);
        TS_ASSERT_DELTA(vec_repl[8], 0.0, 1e-9);
        TS_ASSERT_DELTA(vec_repl[11], 0.0, 1e-9);
        TS_ASSERT_DELTA(vec_repl[14], 0.0, 1e-9);
        TS_ASSERT_DELTA(vec_repl[17], 0.0, 1e-9);

        // Replace the body force with a functional body force (see MyBodyForce) above, and
        // assemble the vector again. This bit isn't so much to test the vector, but
        // to test the physical location being integrated at is interpolated correctly
        // and passed into the force function - see asserts in MyBodyForce.
        mesh.Translate(0.5, 0.8);

        Vec vec2 = PetscTools::CreateVec(18);
        assembler.SetVectorToAssemble(vec2, true);
        problem_defn.SetBodyForce(MyBodyForce);
        assembler.Assemble();

        ReplicatableVector vec_repl2(vec2);
        for (unsigned i=3; i<6; i++)
        {
            TS_ASSERT_DELTA(vec_repl2[3*i],   10.0/6.0, 1e-8);
            TS_ASSERT_DELTA(vec_repl2[3*i+1], 20.0/6.0, 1e-8);
        }

        PetscTools::Destroy(vec);
        PetscTools::Destroy(vec2);
        PetscTools::Destroy(mat);
    }
};

#endif // TESTSTOKESFLOWASSEMBLER_HPP_
