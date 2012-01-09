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


#ifndef TESTPSEUDOECGCALCULATORNIGHTLY_HPP_
#define TESTPSEUDOECGCALCULATORNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "TetrahedralMesh.hpp" //must be first, it gets UblasIncludes from the mesh classes (ChastePoint.hpp)
#include "DistributedTetrahedralMesh.hpp"
#include "PseudoEcgCalculator.hpp"
#include "ReplicatableVector.hpp"
#include "Hdf5DataReader.hpp"
#include "Hdf5DataWriter.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "HeartConfig.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestPseudoEcgCalculatorNightly : public CxxTest::TestSuite
{

public:


    void TestCalculatorRealistic3D() throw (Exception)
    {
        //get the mesh, whole heart mesh
        TrianglesMeshReader<3,3> reader("apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles");
        DistributedTetrahedralMesh<3,3> mesh(DistributedTetrahedralMeshPartitionType::DUMB);  //The partition/numbering has to match that stored in the HDF5 file
        mesh.ConstructFromMeshReader(reader);

        // Compute the pseudo ECG. We set an electrode at x=2.2, y=6, z=1.85.
        ChastePoint<3> probe_electrode(2.2, 6.0, 1.85);

        // The file 3D.h5 contains the first 5 time steps of a whole heart simulations known to produce
        // a reasonably-looking ECG trace.
        PseudoEcgCalculator<3,3,1> calculator (mesh, probe_electrode, "heart/test/data/PseudoEcg", "3D", "V", false);

        calculator.SetDiffusionCoefficient(1.0);//the default value

        //where to put results...
        std::string pseudo_ecg_output_dir = "TestRealisticPseudoEcg";
        HeartConfig::Instance()->SetOutputDirectory(pseudo_ecg_output_dir);

        //write out the pseudoECG
        calculator.WritePseudoEcg();

        //now compare it with a valid pseudo-ecg file (see above).
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

        std::string command_second_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + pseudo_ecg_output_dir + "/output/PseudoEcgFromElectrodeAt_2.2_6_1.85.dat"
                                     + " heart/test/data/PseudoEcg/ValidPseudoEcg.dat";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);

    }
    void TestCalculatorRealistic3DNotDistributed() throw (Exception)
    {
         //get the mesh, whole heart mesh
        TrianglesMeshReader<3,3> reader("apps/texttest/weekly/Propagation3d/heart_chaste2_renum_i_triangles");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Compute the pseudo ECG. We set an electrode at x=2.2, y=6, z=1.85.
        ChastePoint<3> probe_electrode(2.2, 6.0, 1.85);

        // The file 3D.h5 contains the first 5 time steps of a whole heart simulations known to produce
        // a reasonably-looking ECG trace.
        PseudoEcgCalculator<3,3,1> calculator (mesh, probe_electrode, "heart/test/data/PseudoEcg", "3D", "V", false);

        calculator.SetDiffusionCoefficient(1.0);//the default value

        //where to put results...
        std::string pseudo_ecg_output_dir = "TestRealisticPseudoEcgNonDist";
        HeartConfig::Instance()->SetOutputDirectory(pseudo_ecg_output_dir);

        //write out the pseudoECG
        calculator.WritePseudoEcg();

        //now compare it with a valid pseudo-ecg file (see above).
        std::string test_output_directory = OutputFileHandler::GetChasteTestOutputDirectory();

        std::string command_second_time_step = "diff -a -I \"Created by Chaste\" " + test_output_directory + pseudo_ecg_output_dir + "/output/PseudoEcgFromElectrodeAt_2.2_6_1.85.dat"
                                     + " heart/test/data/PseudoEcg/ValidPseudoEcg.dat";
        TS_ASSERT_EQUALS(system(command_second_time_step.c_str()), 0);

    }
};


#endif /*TESTPSEUDOECGCALCULATORNIGHTLY_HPP_*/
