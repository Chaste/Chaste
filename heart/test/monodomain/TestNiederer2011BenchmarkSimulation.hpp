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

    AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    {
        double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
        double y = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[1];
        double z = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[2];

        if ( (x<0.15+1e-6) && (y<0.15+1e-6) && (z<0.15+1e-6) )
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
    void RunBenchMark(double h, double dt, double endTime)
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
        HeartConfig::Instance()->SetUseStateVariableInterpolation(true);

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
    void ConvertToActivationMap(double h, double dt)
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

        unsigned the_node = 0; //corner node
        for(unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
         
            double x = mesh.GetNode(i)->rGetLocation()[0];
            double y = mesh.GetNode(i)->rGetLocation()[1];
            double z = mesh.GetNode(i)->rGetLocation()[2];
            if( fabs(x-2.0) + fabs(y-0.7) + fabs(z-0.3) < 1e-8)
            {
                the_node = i;
            }
        } 

        std::vector<double> activation_times(mesh.GetNumNodes(), -1.0);
        std::vector<double> last_negative_voltage(mesh.GetNumNodes(), 1.0);
        for(unsigned timestep=0; timestep<num_timesteps; timestep++)
        {
            reader.GetVariableOverNodes(voltage, "V", timestep);
            ReplicatableVector voltage_repl(voltage);

            for(unsigned i=0; i<mesh.GetNumNodes(); i++)
            {
                double V = voltage_repl[i];
                if(V > 0 && activation_times[i] < 0.0)
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

        if (PetscTools::AmMaster())
        {
            std::cout << "h, dt = " << h << ", " << dt << "\n\t";
            std::cout << "activation_times[" << the_node << "] = " << activation_times[the_node] << "\n";
        }
        OutputFileHandler handler("ActivationMaps", false);
        std::stringstream output_file;
        output_file << "activation" << "_h" << h << "_dt" << dt << ".dat";
        if (PetscTools::AmMaster())
        {
            out_stream p_file = handler.OpenOutputFile(output_file.str());
    
            for(unsigned i=0; i<activation_times.size(); i++)
            {
                *p_file << activation_times[i] << "\n";
            }
            p_file->close();
        }

        for(unsigned i=0; i<activation_times.size(); i++)
        {
            if(activation_times[i] < 0.0)
            {
                std::cout << "\n\n\n**Some nodes unactivated**\n\n\n";
                output_file << "__error";
                out_stream p_file2 = handler.OpenOutputFile(output_file.str());
                p_file2->close();
                return;
            }
        }
    }


    void Run(double h, double dt, double endTime)
    {
        RunBenchMark(h, dt, endTime);
        ConvertToActivationMap(h, dt);
    }

public:
    void TestRunOnCoarsestMesh() throw(Exception)
    {
        // run with h=0.05 and dt=0.01
        Run(0.05, 0.01, 40);

        // test the activation times produced match those given for the Niederer et al paper
        std::string results_dir = OutputFileHandler::GetChasteTestOutputDirectory() + "ActivationMaps";
        std::string output_file = results_dir + "/activation_h0.05_dt0.01.dat";
        std::string base_file = "heart/test/data/benchmark_data_orig/activation_time_h0.05_dt0.01.dat";
   
        NumericFileComparison num_comp(output_file, base_file);
        TS_ASSERT(num_comp.CompareFiles(1.5e-3)); //Absolute difference of 1.5 microsecond is tolerated
         
    }
};

#endif // TESTNIEDERER2011BENCHMARKSIMULATION_HPP_
