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
#ifndef TESTPAPILLARYFIBRECALCULATOR_HPP_
#define TESTPAPILLARYFIBRECALCULATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "PapillaryFibreCalculator.hpp"

class TestPapillaryFibreCalculator : public CxxTest::TestSuite
{
public:

    void TestGetSingleRadiusVector(void)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/simple_cube");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumElements(),12u);

        PapillaryFibreCalculator calculator(mesh);

        // Call GetRadiusVectors on an element
        unsigned element_index = 0;
        c_vector<double, 3> radius_vector = calculator.GetRadiusVectorForOneElement(element_index);

        // Check they are right
        TS_ASSERT_DELTA(radius_vector[0], -0.275, 1e-9);
        TS_ASSERT_DELTA(radius_vector[1], -0.025, 1e-9);
        TS_ASSERT_DELTA(radius_vector[2], -0.275, 1e-9);
    }

    void TestGetRadiusVectorsAndConstructStructureTensors(void)
    {
        // Test in three parts to use the results of one test in the next...
        //
        //

        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/simple_cube");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_EQUALS(mesh.GetNumElements(),12u);

        PapillaryFibreCalculator calculator(mesh);

        // Call GetRadiusVectors on an element

        calculator.GetRadiusVectors();

        std::vector< c_vector<double,3> >& radius_vectors = calculator.mRadiusVectors;

        // Check they are right
        TS_ASSERT_DELTA(radius_vectors[0][0], -0.275, 1e-9);
        TS_ASSERT_DELTA(radius_vectors[0][1], -0.025, 1e-9);
        TS_ASSERT_DELTA(radius_vectors[0][2], -0.275, 1e-9);

        TS_ASSERT_DELTA(radius_vectors[5][0], 0.475, 1e-9);
        TS_ASSERT_DELTA(radius_vectors[5][1], 0.225, 1e-9);
        TS_ASSERT_DELTA(radius_vectors[5][2], 0.475, 1e-9);

        TS_ASSERT_EQUALS(radius_vectors.size(), mesh.GetNumElements());

        ///////////////////////////////////////////////////////////
        // Test ConstructStructureTensors()
        ///////////////////////////////////////////////////////////
        calculator.ConstructStructureTensors();
        std::vector< c_matrix<double,3,3> >& tensor_i = calculator.mStructureTensors;

        // Worked out by hand...
        TS_ASSERT_DELTA(tensor_i[0](0,0),7.5625e-02,1e-9);
        TS_ASSERT_DELTA(tensor_i[0](0,1),6.8750e-03,1e-9);
        TS_ASSERT_DELTA(tensor_i[0](0,2),7.5625e-02,1e-9);
        TS_ASSERT_DELTA(tensor_i[0](1,0),6.8750e-03,1e-9);
        TS_ASSERT_DELTA(tensor_i[0](1,1),6.2500e-04,1e-9);
        TS_ASSERT_DELTA(tensor_i[0](1,2),6.8750e-03,1e-9);
        TS_ASSERT_DELTA(tensor_i[0](2,0),7.5625e-02,1e-9);
        TS_ASSERT_DELTA(tensor_i[0](2,1),6.8750e-03,1e-9);
        TS_ASSERT_DELTA(tensor_i[0](2,2),7.5625e-02,1e-9);

        TS_ASSERT_DELTA(tensor_i[5](0,0),0.225625,1e-9);
        TS_ASSERT_DELTA(tensor_i[5](0,1),0.106875,1e-9);
        TS_ASSERT_DELTA(tensor_i[5](0,2),0.225625,1e-9);
        TS_ASSERT_DELTA(tensor_i[5](1,0),0.106875,1e-9);
        TS_ASSERT_DELTA(tensor_i[5](1,1),0.050625,1e-9);
        TS_ASSERT_DELTA(tensor_i[5](1,2),0.106875,1e-9);
        TS_ASSERT_DELTA(tensor_i[5](2,0),0.225625,1e-9);
        TS_ASSERT_DELTA(tensor_i[5](2,1),0.106875,1e-9);
        TS_ASSERT_DELTA(tensor_i[5](2,2),0.225625,1e-9);

        //////////////////////////////////////////////////////////////
        // Test SmoothStructureTensor()
        //////////////////////////////////////////////////////////////
        calculator.SmoothStructureTensors();
        std::vector< c_matrix<double,3,3> >& tensor_smooth = calculator.mSmoothedStructureTensors;

        // hard-coded (as difficult to test)
        TS_ASSERT_DELTA(tensor_smooth[0](0,0), 0.075625, 1e-5);
    }

    // see also TestPapillaryFibreCalculatorLong() for bigger test of the main method on a cylinder.
    void TestGetFibreOrientationsOnSimpleCube(void)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/simple_cube");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PapillaryFibreCalculator calculator(mesh);
        std::vector<c_vector<double,3> > fibre_orientations = calculator.CalculateFibreOrientations();

        // Not very well defined for a cube (since it's so well structured that there are zero
        // eigenvectors in the smoothed tensors) but necessary to test coverage.
        // Nightly test TestPapillaryFibreCalculatorLong.hpp is a better one if you want to understand it!

        ///\todo There may still be a sign issue between flapack and MKL
        ///\todo THIS TEST IS KNOWN TO STALL IN LAPACK (dgeev_)
        //       ON THIS CONFIGURATION :
        //  * 32-bit virtual machine on Ubuntu 8.04 LTS
        //  * build=GccOpt
        //  * PETSc: petsc-2.3.2-p10 with f2cblaslapack
        TS_ASSERT_DELTA(fabs(fibre_orientations[0](0)), 0.7056, 1e-4);
        TS_ASSERT_DELTA(fabs(fibre_orientations[0](1)), 0.0641, 1e-4);
        TS_ASSERT_DELTA(fabs(fibre_orientations[0](2)), 0.7056, 1e-4);

        TS_ASSERT_DELTA(fabs(fibre_orientations[4](0)), 0.0455, 1e-4);
        TS_ASSERT_DELTA(fabs(fibre_orientations[4](1)), 0.5005, 1e-4);
        TS_ASSERT_DELTA(fabs(fibre_orientations[4](2)),  0.8645, 1e-4);

        TS_ASSERT_DELTA(fabs(fibre_orientations[5](0)), 0.6704, 1e-4);
        TS_ASSERT_DELTA(fabs(fibre_orientations[5](1)), 0.3176, 1e-4);
        TS_ASSERT_DELTA(fabs(fibre_orientations[5](2)), 0.6704, 1e-4);
    }
};

#endif /*TESTPAPILLARYFIBRECALCULATOR_HPP_*/
