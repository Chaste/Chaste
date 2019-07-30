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

#ifndef TESTNIEDERER2011BENCHMARKSIMULATION_HPP_
#define TESTNIEDERER2011BENCHMARKSIMULATION_HPP_


#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "TenTusscher2006Epi.hpp"
#include "TenTusscher2006EpiBackwardEuler.hpp"
#include "NumericFileComparison.hpp"

/* This test provides the code for, and runs one simulation of, the benchmark problem defined in
 * the paper:
 *
 * S. Niederer et al. Verification of cardiac tissue electrophysiology simulators using an N-version
 * benchmark. Royal Society Philosophical Transactions A, 369(1954):4331-4351, 2011
 *
 * This paper set-up an unambiguously-defined problem, which is then run using eleven large-scale
 * cardiac electro-physiology solvers, including Chaste, in order to verify that they give the same or
 * similar results.
 *
 * The benchmark problem uses: the monodomain equations, a cuboid geometry, anisotropic conductivies,
 * a ten Tusscher cell model, stimulation in a small region in one corner of the mesh, and 3 choices
 * of spatial stepsize
 * (h = 0.01, 0.02, 0.05 cm) and 3 choices of timestep. This test sets up the benchmark simulation with
 * one choice of h (=0.05cm) and dt (=0.01ms). Activation times for each node are computed and compared
 * to the results provided for the above paper.
 *
 * Note: the obtain the results in this paper, the simulations are run with 'state-variable interpolation'
 * (SVI) switched on (which is not true by default). For full details on SVI, see the documentation page:
 * ChasteGuides -> Advanced -> StateVariableInterpolation
 *
 * Also, a follow-up paper
 *
 * P. Pathmanathan, M. Bernabeu, S. Niederer, D. Gavaghan, D. Kay. Computational modelling of cardiac
 * electro-physiology: explanation of the variability of results from different numerical solvers.
 * J. Numerical Methods in Biomedical Engineering (accepted for publication - see publications
 * page of Chaste website for pre-print).
 *
 * fully explains the reasons for differences between solvers. It runs various numerical methods on this
 * benchmark problem. Some comments on repeating the experiments in this paper are given below.
 *
 */


class BenchmarkCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    BenchmarkCellFactory()
      : AbstractCardiacCellFactory<3>(),
        mpStimulus(new SimpleStimulus(-50000.0, 2))
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];

        if ((x<0.15+1e-6) && (y<0.15+1e-6) && (z<0.15+1e-6))
        {
            return new CellTenTusscher2006EpiFromCellMLBackwardEuler(mpSolver, mpStimulus);
        }
        else
        {
            return new CellTenTusscher2006EpiFromCellMLBackwardEuler(mpSolver, mpZeroStimulus);
        }
    }
};

class TestNiederer2011BenchmarkSimulation : public CxxTest::TestSuite
{
private:
    void RunBenchMark(double h, double dt, double endTime, bool useSvi)
    {
        TetrahedralMesh<3,3> mesh;

        mesh.ConstructRegularSlabMesh(h, 2.0, 0.7, 0.3);

        std::stringstream output_dir;
        output_dir << "Benchmark" << "_h" << h << "_dt" << dt;

        HeartConfig::Instance()->SetOutputDirectory(output_dir.str());
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetSimulationDuration(endTime); //ms

        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.005, dt, 0.1);
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400); // 1400 1/cm
        HeartConfig::Instance()->SetCapacitance(1); // 1uF/cm^2

        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(false);

        // The Chaste results for the benchmark paper use STATE-VARIABLE INTERPOLATION switched on
        // (see comments above)
        HeartConfig::Instance()->SetUseStateVariableInterpolation(useSvi);

        // Regarding the second paper described above, to run the simulations with ICI, comment out the
        // above line. To run the simulation with operator splitting, or with (full) mass-lumping,
        // comment out the above SVI line and uncomment one of the below. (Note: half-lumping is not
        // available).
        //HeartConfig::Instance()->SetUseMassLumping(true); // what is described as full-lumping in this paper
        //HeartConfig::Instance()->SetUseReactionDiffusionOperatorSplitting(true);

        double long_conductance = 0.17 * 0.62/(0.17+0.62) * 10;    // harmonic mean of 0.17, 0.62 S/m converted to mS/cm
        double trans_conductance = 0.019 * 0.24/(0.019+0.24) * 10; // harmonic mean of 0.019,0.24 S/m converted to mS/cm

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(long_conductance, trans_conductance, trans_conductance));

        BenchmarkCellFactory cell_factory;

        MonodomainProblem<3> problem( &cell_factory );
        problem.SetMesh(&mesh);

        problem.Initialise();
        problem.SetWriteInfo();
        problem.Solve();
    }

    // this method loads the output file from the previous method and computes the activation
    // time (defined as the time V becomes positive) for each node.
    void ConvertToActivationMap(double h, double dt, bool useSvi)
    {
        //TetrahedralMesh<3,3> mesh1;
        //double h1=0.01;    // 0.01, 0.02, 0.05
        //mesh1.ConstructRegularSlabMesh(h1, 2.0, 0.7, 0.3);
        //MeshalyzerMeshWriter<3,3> writer("Mesh0.01", "mesh01");
        //writer.WriteFilesUsingMesh(mesh1);

        TetrahedralMesh<3,3> mesh;
        double printing_dt=0.1;
        mesh.ConstructRegularSlabMesh(h, 2.0, 0.7, 0.3);

        std::stringstream input_dir;
        input_dir << "Benchmark" << "_h" << h << "_dt" << dt;
        Hdf5DataReader reader(input_dir.str(),"results");

        unsigned num_timesteps = reader.GetUnlimitedDimensionValues().size();
        DistributedVectorFactory factory(mesh.GetNumNodes());
        Vec voltage = factory.CreateVec();


        std::vector<double> activation_times(mesh.GetNumNodes(), -1.0);
        std::vector<double> last_negative_voltage(mesh.GetNumNodes(), 1.0);
        for (unsigned timestep=0; timestep<num_timesteps; timestep++)
        {
            reader.GetVariableOverNodes(voltage, "V", timestep);
            ReplicatableVector voltage_repl(voltage);

            for (unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double V = voltage_repl[i];
                if (V > 0 && activation_times[i] < 0.0)
                {
                    double old = last_negative_voltage[i];
                    assert(old < 0);
                    activation_times[i] = (timestep-V/(V-old))*printing_dt;
                }
                else if (V<=0)
                {
                    last_negative_voltage[i]=V;
                }
            }
        }

        OutputFileHandler handler("ActivationMaps", false);
        if (PetscTools::AmMaster() == false)
        {
            return;
        }
        //Only master proceeds to write


        c_vector<double, 3> top_corner;
        top_corner[0] = 2.0;
        top_corner[1] = 0.7;
        top_corner[2] = 0.3;
        c_vector<double, 3> unit_diagonal = top_corner/norm_2(top_corner);

        std::stringstream output_file1;
        output_file1 << "diagonal" << "_h" << h << "_dt" << dt;
        if (useSvi)
        {
            output_file1 << "_svi.dat";
        }
        else
        {
            output_file1 << "_ici.dat";
        }
        out_stream p_diag_file = handler.OpenOutputFile(output_file1.str());

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            c_vector<double, 3> position =  mesh.GetNode(i)->rGetLocation();
            c_vector<double, 3>  projected_diagonal = unit_diagonal*inner_prod(unit_diagonal, position);
            double off_diagonal = norm_2(position - projected_diagonal);

            if (off_diagonal < h/3)
            {
                double distance = norm_2(position);
                (*p_diag_file) << distance<<"\t"<< activation_times[i]<<"\t"<<off_diagonal<<"\n";
                if (fabs(position[0]-2.0) < 1e-8)
                {
                    std::cout << "h, dt = " << h << ", " << dt << "\n\t";
                    std::cout << "activation_times[" << i << "] = " << activation_times[i] << "\n";
                }
            }
        }
        p_diag_file->close();

        std::stringstream output_file;
        output_file << "activation" << "_h" << h << "_dt" << dt << ".dat";
        out_stream p_file = handler.OpenOutputFile(output_file.str());

        for (unsigned i=0; i<activation_times.size(); i++)
        {
            *p_file << activation_times[i] << "\n";
        }
        p_file->close();

        for (unsigned i=0; i<activation_times.size(); i++)
        {
            if (activation_times[i] < 0.0)
            {
                std::cout << "\n\n\n**Some nodes unactivated**\n\n\n";
                output_file << "__error";
                out_stream p_file2 = handler.OpenOutputFile(output_file.str());
                p_file2->close();
                return;
            }
        }
    }

    void Run(double h, double dt, double endTime, bool useSvi)
    {
        RunBenchMark(h, dt, endTime, useSvi);
        ConvertToActivationMap(h, dt, useSvi);
    }

public:
    void TestRunOnCoarsestMesh()
    {
        // run with h=0.05 and dt=0.01.  SVI is turned on
        Run(0.05, 0.01, 40, true);

        // test the activation times produced match those given for the Niederer et al paper
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "ActivationMaps";
        std::string output_file = results_dir + "/activation_h0.05_dt0.01.dat";
        std::string base_file = "heart/test/data/benchmark_data_orig/activation_time_h0.05_dt0.01.dat";

        NumericFileComparison num_comp(output_file, base_file);
        TS_ASSERT(num_comp.CompareFiles(1.5e-3)); //Absolute difference of 1.5 microsecond is tolerated
    }
    void donotTestRunOtherSvi()
    {
        /* To reproduce Code A panel of
         *  Figure 1 Pathmanathan et al. == Figure 2 Niederer et al.
         * we require
         * Run(0.05, 0.005, 40, true);
         * Run(0.02, 0.005, 50, true);
         * Run(0.01, 0.005, 50, true);
         *
         */
        /* This (together with the above) produces Pathmanathan et al. Figure 1 panel A "SVI"
         *
         * NB. Verify time-step
         */

        Run(0.02, 0.01, 50, true);
        Run(0.01, 0.01, 50, true);
    }

    void donotTestRunAllIci()
    {
        /* This produces Pathmanathan et al. Figure 1 panel B "ICI"
         *
         * NB. Verify time-step
         */
        Run(0.05, 0.01, 60, false);
        Run(0.02, 0.01, 50, false);
        Run(0.01, 0.01, 50, false);
    }
};

#endif // TESTNIEDERER2011BENCHMARKSIMULATION_HPP_
