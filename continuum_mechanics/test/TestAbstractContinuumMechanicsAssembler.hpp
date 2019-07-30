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

#ifndef TESTABSTRACTCONTINUUMMECHANICSASSEMBLER_HPP_
#define TESTABSTRACTCONTINUUMMECHANICSASSEMBLER_HPP_

#include <cxxtest/TestSuite.h>
#include "AbstractContinuumMechanicsAssembler.hpp"
#include "TetrahedralMesh.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscMatTools.hpp"
#include "ReplicatableVector.hpp"
#include "MassMatrixAssembler.hpp"
#include "TrianglesMeshReader.hpp"

///////////////////////////////////////////////////////////////////////////////////
//
//
//   See comments about ordering in AbstractContinuumMechanicsAssembler
//
//
///////////////////////////////////////////////////////////////////////////////////


// NOTE: The tests in TestStokesFlowAssembler also, essentially, test that this class is performing
// correctly.


template<unsigned DIM>
class SimpleAssembler : public AbstractContinuumMechanicsAssembler<DIM,true,true>
{
private:
    static const unsigned NUM_VERTICES_PER_ELEMENT = DIM+1;

    /** Number of nodes per element. */
    static const unsigned NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2; // assuming quadratic

    static const unsigned SPATIAL_BLOCK_SIZE_ELEMENTAL = DIM*NUM_NODES_PER_ELEMENT;
    static const unsigned PRESSURE_BLOCK_SIZE_ELEMENTAL = NUM_VERTICES_PER_ELEMENT;

    double mVal1;
    double mVal2;
    double mVal3;
    double mVal4;
    double mVal5;

public:
    SimpleAssembler(QuadraticMesh<DIM>* pMesh, double val1, double val2, double val3, double val4=0, double val5=0)
        : AbstractContinuumMechanicsAssembler<DIM,true,true>(pMesh),
          mVal1(val1),
          mVal2(val2),
          mVal3(val3),
          mVal4(val4),
          mVal5(val5)
    {
    }

    c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL> ComputeSpatialSpatialMatrixTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,SPATIAL_BLOCK_SIZE_ELEMENTAL> ret;
        for (unsigned i=0; i<SPATIAL_BLOCK_SIZE_ELEMENTAL; i++)
        {
            for (unsigned j=0; j<SPATIAL_BLOCK_SIZE_ELEMENTAL; j++)
            {
                ret(i,j) = mVal1;
            }
        }
        return ret;
    }

    c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputeSpatialPressureMatrixTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
        c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double,SPATIAL_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ret;
        for (unsigned i=0; i<SPATIAL_BLOCK_SIZE_ELEMENTAL; i++)
        {
            for (unsigned j=0; j<PRESSURE_BLOCK_SIZE_ELEMENTAL; j++)
            {
                ret(i,j) = mVal2;
            }
        }
        return ret;
    }

    c_matrix<double,PRESSURE_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputePressurePressureMatrixTerm(
        c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
        c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double,PRESSURE_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ret;
        for (unsigned i=0; i<PRESSURE_BLOCK_SIZE_ELEMENTAL; i++)
        {
            for (unsigned j=0; j<PRESSURE_BLOCK_SIZE_ELEMENTAL; j++)
            {
                ret(i,j) = mVal3;
            }
        }
        return ret;
    }


    c_vector<double,SPATIAL_BLOCK_SIZE_ELEMENTAL> ComputeSpatialVectorTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        c_vector<double,SPATIAL_BLOCK_SIZE_ELEMENTAL> ret;
        for (unsigned i=0; i<SPATIAL_BLOCK_SIZE_ELEMENTAL; i++)
        {
            ret(i) = mVal4;
        }
        return ret;
    }

    c_vector<double,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputePressureVectorTerm(
            c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
            c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
            c_vector<double,DIM>& rX,
            Element<DIM,DIM>* pElement)
    {
        c_vector<double,PRESSURE_BLOCK_SIZE_ELEMENTAL> ret;
        for (unsigned i=0; i<PRESSURE_BLOCK_SIZE_ELEMENTAL; i++)
        {
            ret(i) = mVal5;
        }
        return ret;
    }
};

// Doesn't over-ride any non-compulsory methods, so should return a zero matrix.
class ZeroMatrixAssembler : public AbstractContinuumMechanicsAssembler<1,false,true>
{
public:
    ZeroMatrixAssembler(QuadraticMesh<1>* pMesh)
        : AbstractContinuumMechanicsAssembler<1,false,true>(pMesh)
    {
    }
    c_vector<double,SPATIAL_BLOCK_SIZE_ELEMENTAL> ComputeSpatialVectorTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, 1, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double,1>& rX,
        Element<1,1>* pElement)
    {
        return zero_vector<double>(SPATIAL_BLOCK_SIZE_ELEMENTAL);
    }
};

template<unsigned DIM>
class MyMatrixAssembler : public AbstractContinuumMechanicsAssembler<DIM,false,true>
{
private:
    static const unsigned NUM_VERTICES_PER_ELEMENT = DIM+1;

    /** Number of nodes per element. */
    static const unsigned NUM_NODES_PER_ELEMENT = (DIM+1)*(DIM+2)/2; // assuming quadratic

    static const unsigned SPATIAL_BLOCK_SIZE_ELEMENTAL = DIM*NUM_NODES_PER_ELEMENT;
    static const unsigned PRESSURE_BLOCK_SIZE_ELEMENTAL = NUM_VERTICES_PER_ELEMENT;

public:
    MyMatrixAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh)
        : AbstractContinuumMechanicsAssembler<DIM,false,true>(pMesh)
    {
    }

    c_matrix<double,PRESSURE_BLOCK_SIZE_ELEMENTAL,PRESSURE_BLOCK_SIZE_ELEMENTAL> ComputePressurePressureMatrixTerm(
        c_vector<double, NUM_VERTICES_PER_ELEMENT>& rLinearPhi,
        c_matrix<double, DIM, NUM_VERTICES_PER_ELEMENT>& rGradLinearPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        return outer_prod(rLinearPhi,rLinearPhi);
    }
    c_vector<double,SPATIAL_BLOCK_SIZE_ELEMENTAL> ComputeSpatialVectorTerm(
        c_vector<double, NUM_NODES_PER_ELEMENT>& rQuadPhi,
        c_matrix<double, DIM, NUM_NODES_PER_ELEMENT>& rGradQuadPhi,
        c_vector<double,DIM>& rX,
        Element<DIM,DIM>* pElement)
    {
        return zero_vector<double>(SPATIAL_BLOCK_SIZE_ELEMENTAL);
    }
};

class TestAbstractContinuumMechanicsAssembler : public CxxTest::TestSuite
{
public:
    void TestAssemblers1d()
    {
        double h=0.1;
        QuadraticMesh<1> mesh(h,h); // require a one-element mesh
        unsigned size = 6;

        Vec vec = PetscTools::CreateVec(size);
        Mat mat;
        PetscTools::SetupMat(mat, size, size, size);

        SimpleAssembler<1> assembler(&mesh, 2.0, 3.0, 4.0, 111.0, 222.0);

        // cover exceptions
        TS_ASSERT_THROWS_THIS(assembler.AssembleVector(), "Vector to be assembled has not been set");
        TS_ASSERT_THROWS_THIS(assembler.AssembleMatrix(), "Matrix to be assembled has not been set");


        assembler.SetVectorToAssemble(vec, true);
        assembler.SetMatrixToAssemble(mat, true);
        assembler.Assemble();
        PetscMatTools::Finalise(mat);

        ReplicatableVector vec_repl(vec);
        for (unsigned i=0; i<3; i++)
        {
            // spatial vars
            TS_ASSERT_DELTA( vec_repl[2*i], 111*h, 1e-8 );
        }
        for (unsigned i=0; i<2; i++)
        {
            // non-dummy pressure vars
            TS_ASSERT_DELTA( vec_repl[2*i+1], 222*h, 1e-8 );
        }
        // dummy pressure vars
        TS_ASSERT_DELTA( vec_repl[5], 0.0, 1e-8 );


        double correct_matrix[6][6] = { {2*h, 3*h, 2*h, 3*h, 2*h, 0.0},
                                        {3*h, 4*h, 3*h, 4*h, 3*h, 0.0},
                                        {2*h, 3*h, 2*h, 3*h, 2*h, 0.0},
                                        {3*h, 4*h, 3*h, 4*h, 3*h, 0.0},
                                        {2*h, 3*h, 2*h, 3*h, 2*h, 0.0},
                                        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0} };


        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<6; j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(mat,i,j), correct_matrix[i][j], 1e-8 );
            }
        }

        ZeroMatrixAssembler zero_matrix_assembler(&mesh);
        zero_matrix_assembler.SetMatrixToAssemble(mat, true);
        zero_matrix_assembler.Assemble();
        PetscMatTools::Finalise(mat);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<5; j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(mat,i,j), 0.0, 1e-8 );
            }
        }

        Vec bad_size_vec = PetscTools::CreateVec(67);
        Mat bad_size_mat;
        PetscTools::SetupMat(bad_size_mat, 67, 67, 10);

        assembler.SetVectorToAssemble(bad_size_vec, true);
        TS_ASSERT_THROWS_CONTAINS(assembler.AssembleVector(), "Vector provided to be assembled has size 67, not expected size of 6");

        assembler.SetMatrixToAssemble(bad_size_mat, true);
        TS_ASSERT_THROWS_CONTAINS(assembler.AssembleMatrix(), "Matrix provided to be assembled has size 67, not expected size of 6");

        PetscTools::Destroy(mat);
        PetscTools::Destroy(vec);
        PetscTools::Destroy(bad_size_mat);
        PetscTools::Destroy(bad_size_vec);
    }

    // same to main part of TestAssemblers1d except 2d
    void TestAssemblers2d()
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2d_single_triangular_element_quadratic",2,2,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        Vec vec = PetscTools::CreateVec(3*mesh.GetNumNodes());
        Mat mat;
        PetscTools::SetupMat(mat, 3*mesh.GetNumNodes(),3*mesh.GetNumNodes(), 3*mesh.GetNumNodes());

        SimpleAssembler<2> assembler(&mesh, 2.0, 3.0, 4.0, 111.0, 222.0);
        assembler.SetVectorToAssemble(vec, true);
        assembler.SetMatrixToAssemble(mat, true);
        assembler.Assemble();
        PetscMatTools::Finalise(mat);

        ReplicatableVector vec_repl(vec);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA( vec_repl[3*i  ], 111*0.5, 1e-8 );
            TS_ASSERT_DELTA( vec_repl[3*i+1], 111*0.5, 1e-8 );
        }
        for (unsigned i=0; i<mesh.GetNumVertices(); i++)
        {
            TS_ASSERT_DELTA( vec_repl[3*i+2], 222*0.5, 1e-8 );
        }

        // set up correct matrix: disp-disp entries have value 2h, disp-pressure entries
        // have val 3h, and pressure-pressure entries have value 4h....
        c_matrix<double,18,18> correct_matrix;
        for (unsigned i=0; i<6; i++)
        {
            for (unsigned j=0; j<6; j++)
            {
                correct_matrix(3*i,  3*j)   = 2.0*0.5; // 0.5 is area of triangle
                correct_matrix(3*i+1,3*j)   = 2.0*0.5; // 0.5 is area of triangle
                correct_matrix(3*i,  3*j+1) = 2.0*0.5; // 0.5 is area of triangle
                correct_matrix(3*i+1,3*j+1) = 2.0*0.5; // 0.5 is area of triangle

                correct_matrix(3*i,  3*j+2) = 3.0*0.5; // 0.5 is area of triangle
                correct_matrix(3*i+1,3*j+2) = 3.0*0.5; // 0.5 is area of triangle
                correct_matrix(3*i+2,3*j)   = 3.0*0.5; // 0.5 is area of triangle
                correct_matrix(3*i+2,3*j+1) = 3.0*0.5; // 0.5 is area of triangle

                correct_matrix(3*i+2,3*j+2) = 4.0*0.5;
            }
        }

        //...except for the fact that any entry corresponding to pressure on an non-internal node
        //will have value 0
        for (unsigned i=3; i<6; i++)
        {
            for (unsigned j=0; j<18; j++)
            {
                correct_matrix(3*i+2,j) = 0.0;
                correct_matrix(j,3*i+2) = 0.0;
            }
        }

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<18; j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(mat,i,j), correct_matrix(i,j), 1e-8 );
            }
        }

        PetscTools::Destroy(vec);
        PetscTools::Destroy(mat);
    }

    // same as main part of TestAssemblers1d except 3d
    void TestAssemblers3d()
    {
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        Vec vec = PetscTools::CreateVec(4*mesh.GetNumNodes());
        Mat mat;
        PetscTools::SetupMat(mat, 4*mesh.GetNumNodes(), 4*mesh.GetNumNodes(), 4*mesh.GetNumNodes());

        double vol = mesh.GetVolume(); // volume of element equals volume of mesh

        SimpleAssembler<3> assembler(&mesh, 5.0/vol, 10.0/vol, 15.0/vol, 111.0/vol, 222.0/vol);
        assembler.SetVectorToAssemble(vec, true);
        assembler.SetMatrixToAssemble(mat, true);
        assembler.Assemble();
        PetscMatTools::Finalise(mat);

        ReplicatableVector vec_repl(vec);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA( vec_repl[4*i  ], 111, 1e-8 );
            TS_ASSERT_DELTA( vec_repl[4*i+1], 111, 1e-8 );
            TS_ASSERT_DELTA( vec_repl[4*i+2], 111, 1e-8 );
        }
        for (unsigned i=0; i<mesh.GetNumVertices(); i++)
        {
            TS_ASSERT_DELTA( vec_repl[4*i+3], 222, 1e-8 );
        }

        c_matrix<double,40,40> correct_matrix;
        // set up correct matrix: disp-disp entries have value 5, disp-pressure entries
        // have val 10, and pressure-pressure entries have value 15....
        for (unsigned i=0; i<10; i++)
        {
            for (unsigned j=0; j<10; j++)
            {
                correct_matrix(4*i,  4*j)   = 5;
                correct_matrix(4*i+1,4*j)   = 5;
                correct_matrix(4*i+2,4*j)   = 5;
                correct_matrix(4*i,  4*j+1) = 5;
                correct_matrix(4*i+1,4*j+1) = 5;
                correct_matrix(4*i+2,4*j+1) = 5;
                correct_matrix(4*i,  4*j+2) = 5;
                correct_matrix(4*i+1,4*j+2) = 5;
                correct_matrix(4*i+2,4*j+2) = 5;

                correct_matrix(4*i,  4*j+3) = 10;
                correct_matrix(4*i+1,4*j+3) = 10;
                correct_matrix(4*i+2,4*j+3) = 10;

                correct_matrix(4*i+3,4*j)   = 10;
                correct_matrix(4*i+3,4*j+1) = 10;
                correct_matrix(4*i+3,4*j+2) = 10;

                correct_matrix(4*i+3,4*j+3) = 15;
            }
        }

        //...except for the fact that any entry corresponding to pressure on an non-internal node
        //will have value 0
        for (unsigned i=4; i<10; i++)
        {
            for (unsigned j=0; j<40; j++)
            {
                correct_matrix(4*i+3,j) = 0.0;
                correct_matrix(j,4*i+3) = 0.0;
            }
        }

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<34; j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(mat,i,j), correct_matrix(i,j), 1e-8 );
            }
        }

        PetscTools::Destroy(vec);
        PetscTools::Destroy(mat);
    }

    void TestWithMassMatrixInPressurePressureBlock()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> linear_mesh;
        linear_mesh.ConstructFromMeshReader(reader);

        QuadraticMesh<2> quadratic_mesh;
        quadratic_mesh.ConstructFromLinearMeshReader(reader);

        assert(quadratic_mesh.GetNumVertices()==linear_mesh.GetNumNodes());

        Mat mat1;
        PetscTools::SetupMat(mat1, 3*quadratic_mesh.GetNumNodes(), 3*quadratic_mesh.GetNumNodes(), 27);
        MyMatrixAssembler<2> assembler(&quadratic_mesh);
        assembler.SetMatrixToAssemble(mat1, true);
        assembler.Assemble();
        PetscMatTools::Finalise(mat1);

        Mat mat2;
        PetscTools::SetupMat(mat2, linear_mesh.GetNumNodes(), linear_mesh.GetNumNodes(), 4);
        MassMatrixAssembler<2,2> mass_matrix_assembler_linear(&linear_mesh);
        mass_matrix_assembler_linear.SetMatrixToAssemble(mat2, true);
        mass_matrix_assembler_linear.Assemble();
        PetscMatTools::Finalise(mat2);

        // Checking that mat1(2*num_nodes + i, 2*num_nodes + j) = mat2(i,j) directly would be a pain in parallel.
        // Hardcoded the correct mass matrix here, then check both matrices agree with it.
        double correct_matrix[4][4] = { { 0.0833333, 0.0416667, 0, 0.0416667 },
                                        { 0.0416667, 0.166667,  0.0416667, 0.0833333 },
                                        { 0,  0.0416667,  0.0833333, 0.0416667},
                                        { 0.0416667, 0.0833333, 0.0416667, 0.166667} };

        // verify the above matrix is correctly typed.
        int lo, hi;
        MatGetOwnershipRange(mat2, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<4; j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(mat2,i,j), correct_matrix[i][j], 1e-5 );
            }
        }

        MatGetOwnershipRange(mat1, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<3*quadratic_mesh.GetNumNodes(); j++)
            {
                if (i%3==2 && j%3==2 && (i-2)/3<4 && (j-2)/3<4) // the pressure-pressure block, excl pressure values at non-vertices
                {

                    TS_ASSERT_DELTA( PetscMatTools::GetElement(mat1,i,j), correct_matrix[(i-2)/3][(j-2)/3], 1e-5 );
                }
                else
                {
                    TS_ASSERT_DELTA( PetscMatTools::GetElement(mat1,i,j), 0.0, 1e-8 );
                }
            }
        }

        PetscTools::Destroy(mat1);
        PetscTools::Destroy(mat2);
    }

    void TestAbstractContinuumMechanicsAssemblerMeshType()
    {
        TetrahedralMesh<2,2> mesh;
        TS_ASSERT_THROWS_CONTAINS(MyMatrixAssembler<2> assembler2d(&mesh),
                                  "Continuum mechanics assemblers require a quadratic mesh");

        TetrahedralMesh<3,3> mesh3d;
        TS_ASSERT_THROWS_CONTAINS(MyMatrixAssembler<3> assembler3d(&mesh3d),
                                  "Continuum mechanics assemblers require a quadratic mesh");
    }
};

#endif // TESTABSTRACTCONTINUUMMECHANICSASSEMBLER_HPP_
