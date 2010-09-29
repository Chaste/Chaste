/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef TESTFORCES2NOTFORRELEASE_HPP_
#define TESTFORCES2NOTFORRELEASE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "WntConcentration.hpp"
#include "CryptProjectionForce.hpp"
#include "HoneycombMutableVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"

class TestForces2NotForRelease : public AbstractCellBasedTestSuite
{
public:
    void TestForceOutputParameters()
    {
        std::string output_directory = "TestNotForReleaseForceOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with CryptProjectionForce
        CryptProjectionForce projection_force;
        projection_force.SetCutOffLength(1.5);
        TS_ASSERT_EQUALS(projection_force.GetIdentifier(), "CryptProjectionForce");

        out_stream projection_force_parameter_file = output_file_handler.OpenOutputFile("projection_results.parameters");
        projection_force.OutputForceParameters(projection_force_parameter_file);
        projection_force_parameter_file->close();

        std::string projection_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + projection_force_results_dir + "projection_results.parameters notforrelease_cell_based/test/data/TestForcesNotForRelease/projection_results.parameters").c_str()), 0);
    }
};

#endif /*TESTFORCES2NOTFORRELEASE_HPP_*/

