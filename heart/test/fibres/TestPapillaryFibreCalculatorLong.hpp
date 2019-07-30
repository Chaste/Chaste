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


#ifndef TESTPAPILLARYFIBRECALCULATORLONG_HPP_
#define TESTPAPILLARYFIBRECALCULATORLONG_HPP_

#include <cxxtest/TestSuite.h>

#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "UblasIncludes.hpp"
#include "PapillaryFibreCalculator.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestPapillaryFibreCalculatorLong : public CxxTest::TestSuite
{
public:
    void TestGetFibreOrientationsOnCylinder(void)
    {
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cylinder_14748_elem");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PapillaryFibreCalculator calculator(mesh);
        std::vector<c_vector<double,3> > fibre_orientations = calculator.CalculateFibreOrientations();

        for (unsigned i=0; i<fibre_orientations.size(); i++)
        {
            TS_ASSERT_DELTA(fibre_orientations[i](0), 0.0, 0.02);
            TS_ASSERT_DELTA(fibre_orientations[i](1), 0.0, 0.02);
            TS_ASSERT_DELTA(fabs(fibre_orientations[i](2)), 1.0, 1e-3);
        }
    }
};

#endif /*TESTPAPILLARYFIBRECALCULATORLONG_HPP_*/
