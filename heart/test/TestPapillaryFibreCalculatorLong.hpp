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


#ifndef TESTPAPILLARYFIBRECALCULATORLONG_HPP_
#define TESTPAPILLARYFIBRECALCULATORLONG_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "SimpleDataWriter.hpp"
#include "UblasCustomFunctions.hpp"
#include "PapillaryFibreCalculator.hpp"

class TestPapillaryFibreCalculatorLong : public CxxTest::TestSuite
{
public:
    void TestGetFibreOrientationsOnCylinder(void) throw(Exception)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cylinder_14748_elem");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PapillaryFibreCalculator calculator(mesh);
        std::vector<c_vector<double,3> > fibre_orientations = calculator.CalculateFibreOrientations();

        for(unsigned i=0; i<fibre_orientations.size(); i++)
        {
            TS_ASSERT_DELTA(fibre_orientations[i](0), 0.0, 0.02);
            TS_ASSERT_DELTA(fibre_orientations[i](1), 0.0, 0.02);
            TS_ASSERT_DELTA(fabs(fibre_orientations[i](2)), 1.0, 1e-3);
        }
    }
};

#endif /*TESTPAPILLARYFIBRECALCULATORLONG_HPP_*/
