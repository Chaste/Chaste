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


#ifndef TESTPSEUDOECGCALCULATORNIGHTLY_HPP_
#define TESTPSEUDOECGCALCULATORNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "TetrahedralMesh.hpp" //must be first, it gets UblasIncludes from the mesh classes (ChastePoint.hpp)
#include "DistributedTetrahedralMesh.hpp"
#include "PseudoEcgCalculator.hpp"
#include "ReplicatableVector.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "HeartConfig.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "NumericFileComparison.hpp"

class TestPseudoEcgCalculatorNightly : public CxxTest::TestSuite
{
public:

    void TestCalculatorRealistic3D()
    {
        //get the mesh, whole heart mesh
        TrianglesMeshReader<3,3> reader("apps/texttest/weekly/Propagation3d/OxfordRabbitHeart_482um");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);  //The partition/numbering has to match that stored in the HDF5 file
        mesh.ConstructFromMeshReader(reader);

        // Compute the pseudo ECG. We set an electrode at x=2.2, y=6, z=1.85.
        ChastePoint<3> probe_electrode(2.2, 6.0, 1.85);

        // The file 3D.h5 contains the first 5 time steps of a whole heart simulations known to produce
        // a reasonably-looking ECG trace.
        PseudoEcgCalculator<3,3,1> calculator (mesh,
                                               probe_electrode,
                                               FileFinder("heart/test/data/PseudoEcg",RelativeTo::ChasteSourceRoot),
                                               "3D",
                                               "V");

        calculator.SetDiffusionCoefficient(1.0);//the default value

        //where to put results...
        std::string pseudo_ecg_output_dir = "TestRealisticPseudoEcg";
        HeartConfig::Instance()->SetOutputDirectory(pseudo_ecg_output_dir);

        //write out the pseudoECG
        calculator.WritePseudoEcg();

        //now compare it with a valid pseudo-ecg file (see above).
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

        NumericFileComparison comparer(test_output_directory + pseudo_ecg_output_dir + "/output/PseudoEcgFromElectrodeAt_2.2_6_1.85.dat",
                                "heart/test/data/PseudoEcg/ValidPseudoEcg.dat");
        TS_ASSERT(comparer.CompareFiles());
    }

    void TestCalculatorRealistic3DNotDistributed()
    {
         //get the mesh, whole heart mesh
        TrianglesMeshReader<3,3> reader("apps/texttest/weekly/Propagation3d/OxfordRabbitHeart_482um");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Compute the pseudo ECG. We set an electrode at x=2.2, y=6, z=1.85.
        ChastePoint<3> probe_electrode(2.2, 6.0, 1.85);

        // The file 3D.h5 contains the first 5 time steps of a whole heart simulations known to produce
        // a reasonably-looking ECG trace.
        PseudoEcgCalculator<3,3,1> calculator  (mesh,
                                                probe_electrode,
                                                FileFinder("heart/test/data/PseudoEcg",RelativeTo::ChasteSourceRoot),
                                                "3D",
                                                "V");

        calculator.SetDiffusionCoefficient(1.0);//the default value

        //where to put results...
        std::string pseudo_ecg_output_dir = "TestRealisticPseudoEcgNonDist";
        HeartConfig::Instance()->SetOutputDirectory(pseudo_ecg_output_dir);

        //write out the pseudoECG
        calculator.WritePseudoEcg();

        //now compare it with a valid pseudo-ecg file (see above).
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

        NumericFileComparison comparer(test_output_directory + pseudo_ecg_output_dir + "/output/PseudoEcgFromElectrodeAt_2.2_6_1.85.dat",
                                "heart/test/data/PseudoEcg/ValidPseudoEcg.dat");
        TS_ASSERT(comparer.CompareFiles());
    }
};

#endif /*TESTPSEUDOECGCALCULATORNIGHTLY_HPP_*/
