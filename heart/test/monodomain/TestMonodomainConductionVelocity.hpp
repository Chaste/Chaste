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


#ifndef _TESTMONODOMAINCONDUCTIONVELOCITY_HPP_
#define _TESTMONODOMAINCONDUCTIONVELOCITY_HPP_


#include <cxxtest/TestSuite.h>
#include <vector>

#include "TetrahedralMesh.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "Hdf5DataReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "MonodomainProblem.hpp"
#include "CheckMonoLr91Vars.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "LuoRudy1991.hpp"
#include "ActivationOutputModifier.hpp"
#include "NumericFileComparison.hpp"

class TestMonodomainConductionVelocity : public CxxTest::TestSuite
{
public:
    void tearDown()
    {
        HeartConfig::Reset();
    }

    // Solve on a 1D string of cells, 1cm long with a space step of 0.1mm.
    //
    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainDg01DWith100elements()
    {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(30); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoConductionVel");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem(&cell_factory);

        std::vector<unsigned> output_nodes;
        output_nodes.push_back(5);
        output_nodes.push_back(95);
        monodomain_problem.SetOutputNodes(output_nodes);

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        boost::shared_ptr<ActivationOutputModifier> activation_map_0(new ActivationOutputModifier("activation_map_0.0.txt", 0.0));
        boost::shared_ptr<ActivationOutputModifier> activation_map_minus70(new ActivationOutputModifier("activation_map_-70.0.txt", -70.0));
        monodomain_problem.AddOutputModifier(activation_map_0);
        monodomain_problem.AddOutputModifier(activation_map_minus70);

        monodomain_problem.Solve();

        // Extra activation map data
        OutputFileHandler handler("MonoConductionVel", false);
        NumericFileComparison comp_0(handler.GetOutputDirectoryFullPath()+ "activation_map_0.0.txt", "heart/test/data/MonoConductionVel/activation_map_0.0.txt");
        TS_ASSERT(comp_0.CompareFiles());
        NumericFileComparison comp_minus_70(handler.GetOutputDirectoryFullPath()+ "activation_map_-70.0.txt", "heart/test/data/MonoConductionVel/activation_map_-70.0.txt");
        TS_ASSERT(comp_minus_70.CompareFiles());

        // Non-round-robin files
        // We also produce file fragments in case anything happens during post-processing.  These are identical in the sequential case
        if (PetscTools::IsSequential())
        {
            NumericFileComparison comp_0_fragment(handler.GetOutputDirectoryFullPath()+ "activation_map_0.0.txt.0", "heart/test/data/MonoConductionVel/activation_map_0.0.txt");
            TS_ASSERT(comp_0_fragment.CompareFiles());
        }


        // test whether voltages and gating variables are in correct ranges
        CheckMonoLr91Vars<1>(monodomain_problem);

        // Calculate the conduction velocity
        Hdf5DataReader simulation_data=monodomain_problem.GetDataReader();

        PropagationPropertiesCalculator ppc(&simulation_data);
        double velocity=0.0;

        // Check action potential propagated to node 95
        TS_ASSERT_THROWS_NOTHING(velocity=ppc.CalculateConductionVelocity(5,95,0.9));

        // The value should be approximately 50cm/sec
        // i.e. 0.05 cm/msec (which is the units of the simulation)
        TS_ASSERT_DELTA(velocity, 0.05, 0.003);
    }

    // Solve on a 1D string of cells, 1cm long with a space step of 0.5mm.
    //
    // Note that this space step ought to be too big!
    //
    // NOTE: This test uses NON-PHYSIOLOGICAL parameters values (conductivities,
    // surface-area-to-volume ratio, capacitance, stimulus amplitude). Essentially,
    // the equations have been divided through by the surface-area-to-volume ratio.
    // (Historical reasons...)
    void TestMonodomainDg01DWith20elements()
    {
#ifndef NDEBUG //Note that this test relies on the debug VerifyStateVariables() method throwing
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetSimulationDuration(1); //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_20_elements");
        HeartConfig::Instance()->SetOutputDirectory("MonoConductionVelThrows");
        HeartConfig::Instance()->SetOutputFilenamePrefix("MonodomainLR91_1d");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        MonodomainProblem<1> monodomain_problem(&cell_factory);

        monodomain_problem.Initialise();

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1.0);
        HeartConfig::Instance()->SetCapacitance(1.0);

        // the mesh is too coarse, and this simulation will result in cell gating
        // variables going out of range. An exception should be thrown in the
        // EvaluateYDerivatives() method of the cell model

        TS_ASSERT_THROWS_CONTAINS(monodomain_problem.Solve(),
                "State variable fast_sodium_current_m_gate__m has gone out of range. "
                "Check numerical parameters, for example time and space stepsizes");
#endif //NDEBUG //Note that this test relies on the debug VerifyStateVariables() method throwing
    }

};
#endif //_TESTMONODOMAINCONDUCTIONVELOCITY_HPP_
