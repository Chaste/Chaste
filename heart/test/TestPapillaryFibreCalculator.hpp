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
#ifndef TESTPAPILLARYFIBRECALCULATOR_HPP_
#define TESTPAPILLARYFIBRECALCULATOR_HPP_

#include <cxxtest/TestSuite.h>
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "PapillaryFibreCalculator.hpp"

class TestPapillaryFibreCalculator : public CxxTest::TestSuite
{
public:

    void TestGetSingleRadiusVector(void) throw(Exception)
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

    void TestGetRadiusVectorsAndConstructStructureTensors(void) throw(Exception)
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
    void TestGetFibreOrientationsOnSimpleCube(void) throw(Exception)
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
