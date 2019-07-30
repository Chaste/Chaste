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

#ifndef TESTPHASECALCULATIONS_HPP_
#define TESTPHASECALCULATIONS_HPP_

#include <cxxtest/TestSuite.h>

#include "MonodomainProblem.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "LuoRudySpiralWaveCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "FileFinder.hpp"
#include "NumericFileComparison.hpp"
#include "FileComparison.hpp"
#include "Hdf5ToMeshalyzerConverter.hpp"

class TestSpiralWaveAndPhase : public CxxTest::TestSuite
{
private:
    /**
     * Overridden setUp() method. Initialises mesh to be used by any of the tests.
     */
    void setUp()
    {
        mpMesh = new DistributedTetrahedralMesh<2,2>();
        double node_spacing_in_mesh = 0.015;
        mMeshWidth = 3; // cm
        mpMesh->ConstructRegularSlabMesh(node_spacing_in_mesh, mMeshWidth /*length*/, mMeshWidth /*width*/);
    }

    /**
     * Overridden teardown() method. Clears up mesh.
     */
    void tearDown()
    {
        delete mpMesh;
    }

    DistributedTetrahedralMesh<2,2>* mpMesh; /** The mesh to use for all these tests. */
    double mMeshWidth; /** The width of the mesh */

public:

    /**
     * You just need to run this once to get a spiral wave .h5 file.
     * It is put into
     * <CHASTE_TEST_OUTPUT>/SpiralWaveAndPhaseContinued/results.h5
     *
     * This file is now stored in the repository at
     * heart/test/data/PhasePostprocessing/results.h5
     * and is used by the subsequent tests.
     */
    void xTestSpiralWaveSimulationWithPhaseCalculations()
    {
        // Run a simulation to generate a spiral wave in an h5 file.
        {
            HeartConfig::Instance()->SetSimulationDuration(100); //ms
            HeartConfig::Instance()->SetOutputDirectory("SpiralWaveAndPhase");
            HeartConfig::Instance()->SetOutputFilenamePrefix("results");
            HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 1);

            LuoRudySpiralWaveCellFactory cell_factory(mMeshWidth,mMeshWidth);
            MonodomainProblem<2> monodomain_problem( &cell_factory );
            monodomain_problem.SetMesh(mpMesh);
            monodomain_problem.SetWriteInfo();
            monodomain_problem.Initialise();
            monodomain_problem.Solve();

            // Get a shorter results file that is all 'spiral wave'.
            HeartConfig::Instance()->SetOutputDirectory("SpiralWaveAndPhaseContinued");
            HeartConfig::Instance()->SetSimulationDuration(110); //ms
            monodomain_problem.Solve();
        }

        // Get a plot of V(t+tau) against V(t) for each node of above simulation
        // I just used this to estimate a good v_star and tau...
        {
            // Load up the solution
            FileFinder results_directory("SpiralWaveAndPhase",RelativeTo::ChasteTestOutput);
            Hdf5DataReader* p_data_reader = new Hdf5DataReader(results_directory,"results");

            // Get voltage at lots of nodes, NB hi_idx is one past the end.
            std::vector<std::vector<double> > voltages =
                    p_data_reader->GetVariableOverTimeOverMultipleNodes("V",0u, p_data_reader->GetNumberOfRows());

            // Open output files, don't wipe the results we want to examine!
            OutputFileHandler file_handler("SpiralWaveAndPhase",false);
            out_stream voltage_loops = file_handler.OpenOutputFile("voltage_at_each_node.dat");

            std::cout << "nodes size = " << voltages.size() << "\n";
            std::cout << "times size = " << voltages[0].size() << "\n";
            for (unsigned time_idx=0; time_idx<voltages[0].size(); time_idx++)
            {
                for (unsigned node_idx=0; node_idx<voltages.size(); node_idx++)
                {
                    *voltage_loops << voltages[node_idx][time_idx];

                    if (node_idx==voltages.size()-1u)
                    {
                        *voltage_loops << std::endl;
                    }
                    else
                    {
                        *voltage_loops << "\t";
                    }
                }
            }
            voltage_loops->close();
        }
        assert(0); // Should stop test and copy results across now.
    }

    void TestPhaseCalculation()
    {
        // Copy solution from the repository to the relevant test output directory.
        FileFinder reference_voltage("heart/test/data/PhasePostprocessing/results.h5", RelativeTo::ChasteSourceRoot);

        std::string results_folder_name = "SpiralWaveAndPhase";
        std::string results_file_name = "results";

        // Make a clean folder for the results, and load in the existing voltage traces.
        OutputFileHandler file_handler(results_folder_name);
        file_handler.CopyFileTo(reference_voltage);

        // Load up the solution with HDF5 reader.
        FileFinder results_directory(results_folder_name,RelativeTo::ChasteTestOutput);
        Hdf5DataReader* p_data_reader = new Hdf5DataReader(results_directory,results_file_name);
        std::vector<double> times = p_data_reader->GetUnlimitedDimensionValues();

        // Get voltage at lots of nodes, NB hi_idx is one past the end.
        std::vector<std::vector<double> > voltages = p_data_reader->GetVariableOverTimeOverMultipleNodes("V",
                                                                                                0u,
                                                                                                p_data_reader->GetNumberOfRows());

        ///\todo #2137 Figure out how to automate Vstar selection from the voltage loops output by above test.
        double v_star = -40.0; // This will be model-specific(ish). Might work OK as general postprocessing voltage threshold.
        double tau = 2.0; // Time of delay (ms).

        DistributedVectorFactory factory(mpMesh->GetNumNodes());
        Hdf5DataWriter writer(factory, results_folder_name, results_file_name, false, true, "Phase"); // false to wiping, true to extending
        int phase_id = writer.DefineVariable("Phase","radians");
        writer.DefineFixedDimension(mpMesh->GetNumNodes());
        writer.DefineUnlimitedDimension(p_data_reader->GetUnlimitedDimensionName(), p_data_reader->GetUnlimitedDimensionUnit());
        writer.EndDefineMode();

        // Work out the index of the first time greater than or equal to tau.
        unsigned first_index = 0u;
        for (unsigned i=0; i<times.size(); i++)
        {
            if (times[i] >= tau + times[0])
            {
                first_index = i;
                break;
            }
        }

        for (unsigned time_idx=first_index; time_idx<times.size(); time_idx++)
        {
            Vec phase_vec = factory.CreateVec();
            DistributedVector distributed_phase_vector = factory.CreateDistributedVector(phase_vec);
            for (DistributedVector::Iterator index = distributed_phase_vector.Begin();
                 index!= distributed_phase_vector.End();
                 ++index)
            {
                unsigned node_idx = index.Global;
                // Phase calculation, store ready for HDF5 format
                distributed_phase_vector[index] = atan2(voltages[node_idx][time_idx]    -v_star,
                                                        voltages[node_idx][time_idx-tau]-v_star);
            }
            writer.PutVector(phase_id, phase_vec);
            PetscTools::Destroy(phase_vec);
            writer.PutUnlimitedVariable(times[time_idx]);
            writer.AdvanceAlongUnlimitedDimension();
        }

        delete p_data_reader; // You can get some confusing errors in the next test if the reader memory is left alive!
        writer.Close();
    }

    /**
     * Here we test that we can successfully output meshalyzer files of both the original data and Phase.
     */
    void TestMeshalyzerPhasePlot()
    {
        // Write out voltage (just for fun, to compare with phase plot).
        Hdf5ToMeshalyzerConverter<2,2> converter1(FileFinder("SpiralWaveAndPhase", RelativeTo::ChasteTestOutput),
                                                  "results", mpMesh, true);

        // Check the phase file is written correctly.
        {
            FileFinder meshalyzer_phase_file("SpiralWaveAndPhase/output/Phase.dat", RelativeTo::ChasteTestOutput);
            TS_ASSERT(meshalyzer_phase_file.IsFile());
            FileFinder reference_file("heart/test/data/PhasePostprocessing/results_Phase.dat", RelativeTo::ChasteSourceRoot);

            FileComparison comparer(meshalyzer_phase_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Check the times info file is generated correctly.
        {
            FileFinder generated_file("SpiralWaveAndPhase/output/Phase_times.info", RelativeTo::ChasteTestOutput);
            TS_ASSERT(generated_file.IsFile());
            FileFinder reference_file("heart/test/data/PhasePostprocessing/results_Phase_times.info", RelativeTo::ChasteSourceRoot);

            FileComparison comparer(generated_file, reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTPHASECALCULATIONS_HPP_*/
