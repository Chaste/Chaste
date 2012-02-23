/*

Copyright (c) 2005-2012, University of Oxford.
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
#include "PetscSetupAndFinalize.hpp"
#include "PetscMatTools.hpp"
#include "ReplicatableVector.hpp"
#include "MassMatrixAssembler.hpp"
#include "TrianglesMeshReader.hpp"


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
        for(unsigned i=0; i<SPATIAL_BLOCK_SIZE_ELEMENTAL; i++)
        {
            for(unsigned j=0; j<SPATIAL_BLOCK_SIZE_ELEMENTAL; j++)
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
        for(unsigned i=0; i<SPATIAL_BLOCK_SIZE_ELEMENTAL; i++)
        {
            for(unsigned j=0; j<PRESSURE_BLOCK_SIZE_ELEMENTAL; j++)
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
        for(unsigned i=0; i<PRESSURE_BLOCK_SIZE_ELEMENTAL; i++)
        {
            for(unsigned j=0; j<PRESSURE_BLOCK_SIZE_ELEMENTAL; j++)
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
        for(unsigned i=0; i<SPATIAL_BLOCK_SIZE_ELEMENTAL; i++)
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
        for(unsigned i=0; i<PRESSURE_BLOCK_SIZE_ELEMENTAL; i++)
        {
            ret(i) = mVal5;
        }
        return ret;
    }
};

// Doesn't over-ride any methods, so should return a zero matrix. (Note: can't create
// vectors).
class ZeroMatrixAssembler : public AbstractContinuumMechanicsAssembler<1,false,true>
{
public:
    ZeroMatrixAssembler(QuadraticMesh<1>* pMesh)
        : AbstractContinuumMechanicsAssembler<1,false,true>(pMesh)
    {
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
    MyMatrixAssembler(QuadraticMesh<DIM>* pMesh)
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
};

class TestAbstractContinuumMechanicsAssembler : public CxxTest::TestSuite
{
public:
    void TestAssemblers1d() throw (Exception)
    {
        double h=0.1;
        QuadraticMesh<1> mesh(h,h); // require a one-element mesh
        unsigned size = 5; // 1*num_nodes + num_vertices;

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
        for(unsigned i=0; i<3; i++)
        {
            TS_ASSERT_DELTA( vec_repl[i], 111*h, 1e-8 );
        }
        for(unsigned i=3; i<5; i++)
        {
            TS_ASSERT_DELTA( vec_repl[i], 222*h, 1e-8 );
        }


        double correct_matrix[5][5] = { {2*h, 2*h, 2*h, 3*h, 3*h},
                                        {2*h, 2*h, 2*h, 3*h, 3*h},
                                        {2*h, 2*h, 2*h, 3*h, 3*h},
                                        {3*h, 3*h, 3*h, 4*h, 4*h},
                                        {3*h, 3*h, 3*h, 4*h, 4*h} };

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for(unsigned j=0; j<5; j++)
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
            for(unsigned j=0; j<5; j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(mat,i,j), 0.0, 1e-8 );
            }
        }

        Vec bad_size_vec = PetscTools::CreateVec(67);
        Mat bad_size_mat;
        PetscTools::SetupMat(bad_size_mat, 67, 67, 10);

        assembler.SetVectorToAssemble(bad_size_vec, true);
        TS_ASSERT_THROWS_CONTAINS(assembler.AssembleVector(), "Vector provided to be assembled has size 67, not expected size of 5");

        assembler.SetMatrixToAssemble(bad_size_mat, true);
        TS_ASSERT_THROWS_CONTAINS(assembler.AssembleMatrix(), "Matrix provided to be assembled has size 67, not expected size of 5");

        PetscTools::Destroy(mat);
        PetscTools::Destroy(vec);
        PetscTools::Destroy(bad_size_mat);
        PetscTools::Destroy(bad_size_vec);
    }

    // same to main part of TestAssemblers1d except 2d
    void TestAssemblers2d() throw (Exception)
    {
        QuadraticMesh<2> mesh;
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/2d_single_triangular_element_quadratic",2,2,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        Vec vec = PetscTools::CreateVec(2*mesh.GetNumNodes()+mesh.GetNumVertices());
        Mat mat;
        PetscTools::SetupMat(mat, 2*mesh.GetNumNodes()+mesh.GetNumVertices(), 2*mesh.GetNumNodes()+mesh.GetNumVertices(), 2*mesh.GetNumNodes()+mesh.GetNumVertices());

        SimpleAssembler<2> assembler(&mesh, 2.0, 3.0, 4.0, 111.0, 222.0);
        assembler.SetVectorToAssemble(vec, true);
        assembler.SetMatrixToAssemble(mat, true);
        assembler.Assemble();
        PetscMatTools::Finalise(mat);

        ReplicatableVector vec_repl(vec);
        for(unsigned i=0; i<12; i++)
        {
            TS_ASSERT_DELTA( vec_repl[i], 111*0.5, 1e-8 );
        }
        for(unsigned i=13; i<15; i++)
        {
            TS_ASSERT_DELTA( vec_repl[i], 222*0.5, 1e-8 );
        }

        c_matrix<double,15,15> correct_matrix;
        for(unsigned i=0; i<12; i++)
        {
            for(unsigned j=0; j<12; j++)
            {
                correct_matrix(i,j) = 2.0*0.5; // 0.5 is area of triangle
            }
        }

        for(unsigned i=0; i<12; i++)
        {
            for(unsigned j=12; j<15; j++)
            {
                correct_matrix(i,j) = 3.0*0.5; // 0.5 is area of triangle
                correct_matrix(j,i) = 3.0*0.5; // 0.5 is area of triangle
            }
        }

        for(unsigned i=12; i<15; i++)
        {
            for(unsigned j=12; j<15; j++)
            {
                correct_matrix(i,j) = 4.0*0.5; // 0.5 is area of triangle
            }
        }

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for(unsigned j=0; j<15; j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(mat,i,j), correct_matrix(i,j), 1e-8 );
            }
        }

        PetscTools::Destroy(vec);
        PetscTools::Destroy(mat);
    }

    // same as main part of TestAssemblers1d except 3d
    void TestAssemblers3d() throw (Exception)
    {
        QuadraticMesh<3> mesh;
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element_quadratic",2,1,false);
        mesh.ConstructFromMeshReader(mesh_reader);

        Vec vec = PetscTools::CreateVec(3*mesh.GetNumNodes()+mesh.GetNumVertices());
        Mat mat;
        PetscTools::SetupMat(mat, 3*mesh.GetNumNodes()+mesh.GetNumVertices(), 3*mesh.GetNumNodes()+mesh.GetNumVertices(), 3*mesh.GetNumNodes()+mesh.GetNumVertices());

        double vol = mesh.GetVolume(); // volume of element equals volume of mesh

        SimpleAssembler<3> assembler(&mesh, 5.0/vol, 10.0/vol, 15.0/vol, 111.0/vol, 222.0/vol);
        assembler.SetVectorToAssemble(vec, true);
        assembler.SetMatrixToAssemble(mat, true);
        assembler.Assemble();
        PetscMatTools::Finalise(mat);

        ReplicatableVector vec_repl(vec);
        for(unsigned i=0; i<30; i++)
        {
            TS_ASSERT_DELTA( vec_repl[i], 111.0, 1e-8 );
        }
        for(unsigned i=30; i<34; i++)
        {
            TS_ASSERT_DELTA( vec_repl[i], 222.0, 1e-8 );
        }


        c_matrix<double,34,34> correct_matrix;
        for(unsigned i=0; i<30; i++)
        {
            for(unsigned j=0; j<30; j++)
            {
                correct_matrix(i,j) = 5.0;
            }
        }

        for(unsigned i=0; i<30; i++)
        {
            for(unsigned j=30; j<34; j++)
            {
                correct_matrix(i,j) = 10.0;
                correct_matrix(j,i) = 10.0;
            }
        }

        for(unsigned i=30; i<34; i++)
        {
            for(unsigned j=30; j<34; j++)
            {
                correct_matrix(i,j) = 15.0;
            }
        }

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for(unsigned j=0; j<34; j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(mat,i,j), correct_matrix(i,j), 1e-8 );
            }
        }

        PetscTools::Destroy(vec);
        PetscTools::Destroy(mat);
    }

    void TestWithMassMatrixInPressurePressureBlock() throw(Exception)
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_2_elements");
        TetrahedralMesh<2,2> linear_mesh;
        linear_mesh.ConstructFromMeshReader(reader);

        QuadraticMesh<2> quadratic_mesh;
        quadratic_mesh.ConstructFromLinearMeshReader(reader);

        assert(quadratic_mesh.GetNumVertices()==linear_mesh.GetNumNodes());

        Mat mat1;
        PetscTools::SetupMat(mat1, 2*quadratic_mesh.GetNumNodes()+quadratic_mesh.GetNumVertices(), 2*quadratic_mesh.GetNumNodes()+quadratic_mesh.GetNumVertices(), 22);
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
            for(unsigned j=0; j<4; j++)
            {
                TS_ASSERT_DELTA( PetscMatTools::GetElement(mat2,i,j), correct_matrix[i][j], 1e-5 );
            }
        }

        MatGetOwnershipRange(mat1, &lo, &hi);
        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for(unsigned j=0; j<22; j++)
            {
                if(i>=18 && j>=18)
                {
                    TS_ASSERT_DELTA( PetscMatTools::GetElement(mat1,i,j), correct_matrix[i-18][j-18], 1e-5 );
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
};

#endif // TESTABSTRACTCONTINUUMMECHANICSASSEMBLER_HPP_
