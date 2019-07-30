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


#ifndef TESTExtendedBidomainPROBLEM_HPP_
#define TESTExtendedBidomainPROBLEM_HPP_

#include "UblasIncludes.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "HeartConfig.hpp"
#include "SimpleStimulus.hpp"
#include "LuoRudy1991.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "BidomainProblem.hpp"
#include "ExtendedBidomainProblem.hpp"
#include "ExtendedBidomainSolver.hpp"
#include "Hdf5DataReader.hpp"

/**
 * This cell factory is used for the extended bidomain case
 */
class StimulatedCellFactory: public AbstractCardiacCellFactory<1>
{
private:
        boost::shared_ptr<SimpleStimulus> mpStimulus;
        //static const double magnitude = -105000.0;//volume stimulus in microA/cm3
public:
    StimulatedCellFactory() : AbstractCardiacCellFactory<1>(),
            mpStimulus ( new SimpleStimulus(-105000.0, 1.0))/*amplitude, duration (ms)*/
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        CellLuoRudy1991FromCellML* first_cell;
        if (x < 0.005)
        {
            first_cell = new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            first_cell = new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }

        return first_cell;
    }
};

class StimulatedCellFactoryBidomain: public AbstractCardiacCellFactory<1>
{
private:
        boost::shared_ptr<SimpleStimulus> mpStimulus;
        //static const double magnitude = -105000.0;
        // this is a volume stimulus in microA/cm3
public:
        StimulatedCellFactoryBidomain() : AbstractCardiacCellFactory<1>(),
            mpStimulus ( new SimpleStimulus(-105000.0, 1.0))/*amplitude, duration (ms)*/
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        CellLuoRudy1991FromCellML* first_cell;
        if (x < 0.005)
        {
            first_cell = new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        }
        else
        {
            first_cell = new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        }

        return first_cell;
    }
};

/**
 * Unstimulated cell factory used by both bidomain and extended bidomain problems
 */
class UnStimulatedCellFactory: public AbstractCardiacCellFactory<1>
{

public:
    UnStimulatedCellFactory() : AbstractCardiacCellFactory<1>()
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<1>* pNode)
    {
        CellLuoRudy1991FromCellML* second_cell = new CellLuoRudy1991FromCellML(mpSolver, mpZeroStimulus);
        return second_cell;
    }
};

class ExtracellularStimulusFactory: public AbstractStimulusFactory<1>
{

public:
    ExtracellularStimulusFactory() : AbstractStimulusFactory<1>()
    {
    }

    boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(Node<1>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        boost::shared_ptr<SimpleStimulus> p_stimulus;
        if (x < 0.005)
        {
            p_stimulus.reset (new SimpleStimulus(-428, 1.0, 0.1));
        }
        else
        {
            p_stimulus.reset (new SimpleStimulus(0.0, 0.5, 0.1));
        }
        return p_stimulus;
    }
};

/**
 * In these tests we compare the extended bidomain problem and the normal bidomain problem in
 * a situation where they are expected to give almost the same answer.
 *
 * The situation is the following:
 * - the two cells don't talk to each other via gap junction (Ggap = 0, the default value - by the way)
 * - Extracellular conductivity much higher than the intracellular (this dampens the changes in Phi_e caused by one of the cells and so the two cells
 *   are really almost isolated from each other).
 * - No extracellular stimulus.
 *
 *  We apply stimulation only to the one cell and compare with bidomain.
 *  The tolerance for the comparison is not very strict although acceptable. There might two reasons for this
 *  The two cells are 'almost' isolated because the extracellular conductivity is high (but not infinite).
 *
 *   We test both the SetNodeForAverageOfPhiZeroed method and the nullbasis method of solution of the linear system
 */

class TestExtendedVsBidomainProblem : public CxxTest::TestSuite
{

public:

    void SetupParameters()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(65.0));//very high compared to intracellular
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400.0);
        HeartConfig::Instance()->SetOdeTimeStep(0.01);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetSimulationDuration(2.0);  //ms. Tried to run for up tp 35 ms and all was fine.
    }

    //////////////////////////////////////////////////////////////////////////////
    // Tests and simulations for the setAveragePhie method
    //////////////////////////////////////////////////////////////////////////////

    void RunExtendedBidomainStimulateFirstCell()
    {
        SetupParameters();
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);

        StimulatedCellFactory stimulated_cell_factory;
        UnStimulatedCellFactory unstimulated_cell_factory;

        HeartConfig::Instance()->SetOutputDirectory("TestExtendedVs_Extended1dStimulateFirstCell");
        HeartConfig::Instance()->SetOutputFilenamePrefix("extended1d");

        //stimulated cell factory passed in as first cell
        ExtendedBidomainProblem<1> extended_problem( &stimulated_cell_factory , &unstimulated_cell_factory);

        extended_problem.SetWriteInfo(false);
        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(0.0005));//shouldn't be necessary if it is the same as the first cell
        extended_problem.SetNodeForAverageOfPhiZeroed(50u);
        extended_problem.Initialise();
        extended_problem.Solve();
    }

    void RunBidomain()
    {
        SetupParameters();
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);

        HeartConfig::Instance()->SetOutputDirectory("TestExtendedVs_Bidomain1d");
        HeartConfig::Instance()->SetOutputFilenamePrefix("normal1d");
        StimulatedCellFactoryBidomain bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        bidomain_problem.SetNodeForAverageOfPhiZeroed(50u);
        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }

    void TestCompareStimulationOfFirstCell()
    {
        //run the two simulations, bidomain and extended bidomain
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        RunBidomain();
        RunExtendedBidomainStimulateFirstCell();
        // Both the simulations use the "averaged phi" method which gives a non-symmetric matrix and hence uses GMRES
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 2u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Code has changed the KSP solver type from cg to gmres");
        Warnings::Instance()->QuietDestroy();

        //pick up the results...
        Hdf5DataReader reader_extended("TestExtendedVs_Extended1dStimulateFirstCell", "extended1d");
        Hdf5DataReader reader_bidomain("TestExtendedVs_Bidomain1d", "normal1d");
        TS_ASSERT_EQUALS(reader_extended.GetNumberOfRows(), reader_bidomain.GetNumberOfRows());

        bool ap_generated =false;
        //... and compare
        for (unsigned node = 0; node < reader_extended.GetNumberOfRows(); node++)
        {
            std::vector<double> voltage_first_cell_extended = reader_extended.GetVariableOverTime("V", node);
            std::vector<double> voltage_second_cell_extended = reader_extended.GetVariableOverTime("V_2", node);
            std::vector<double> phi_e_extended = reader_extended.GetVariableOverTime("Phi_e", node);

            std::vector<double> voltage_bidomain = reader_bidomain.GetVariableOverTime("V", node);
            std::vector<double> phi_e_bidomain = reader_bidomain.GetVariableOverTime("Phi_e", node);

            TS_ASSERT_EQUALS(voltage_first_cell_extended.size(), voltage_second_cell_extended.size());
            TS_ASSERT_EQUALS(voltage_second_cell_extended.size(), phi_e_extended.size());
            TS_ASSERT_EQUALS(phi_e_extended.size(), voltage_bidomain.size());
            TS_ASSERT_EQUALS(voltage_bidomain.size(), phi_e_bidomain.size());

            for ( unsigned index = 0; index < voltage_bidomain.size(); index ++)
            {
                //check that an AP was generated somewhere at some time
                if (voltage_bidomain[index] >0.0)
                {
                    ap_generated = true;
                }
                TS_ASSERT_DELTA(voltage_first_cell_extended[index] , voltage_bidomain[index], 2e-3);
                TS_ASSERT_LESS_THAN(voltage_second_cell_extended[index] , -83.5);//second unstimulated cell should be at rest
                TS_ASSERT_DELTA(phi_e_extended[index] , phi_e_bidomain[index], 2e-3);
            }
        }
        TS_ASSERT_EQUALS(ap_generated, true);
    }


    void RunExtendedBidomainStimulateSecondCell()
    {
        SetupParameters();
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);

        StimulatedCellFactory stimulated_cell_factory;
        UnStimulatedCellFactory unstimulated_cell_factory;

        HeartConfig::Instance()->SetOutputDirectory("TestExtendedVs_Extended1dStimulateSecondCell");
        HeartConfig::Instance()->SetOutputFilenamePrefix("extended1d");

        //stimulated cell factory passed in as second cell
        ExtendedBidomainProblem<1> extended_problem( &unstimulated_cell_factory,  &stimulated_cell_factory);

        //covering the case where the user sets extended parameters explicitly
        extended_problem.SetExtendedBidomainParameters(1400, 1400, 1400, 1.0 , 1.0, 0.0);

        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(0.0005));
        extended_problem.SetNodeForAverageOfPhiZeroed(50u);
        extended_problem.Initialise();
        extended_problem.Solve();
    }

    /**
     * This test is the same as TestCompareStimulationOfFirstCell above, but with the stimulus applied to the second cell instead of the first.
     */
    void TestCompareStimulationOfSecondCell()
    {
        //first, run the extended bidomain simulation stimulating the second cell
        RunExtendedBidomainStimulateSecondCell();
        // Check that the solver switched to GMRES for the simulation again
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Code has changed the KSP solver type from cg to gmres");
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "gmres")==0);
        Warnings::Instance()->QuietDestroy();

        //pick up the results...(note the bidomain simulatioon is the same as the test above for the first cell)
        Hdf5DataReader reader_extended("TestExtendedVs_Extended1dStimulateSecondCell", "extended1d");
        Hdf5DataReader reader_bidomain("TestExtendedVs_Bidomain1d", "normal1d");
        TS_ASSERT_EQUALS(reader_extended.GetNumberOfRows(), reader_bidomain.GetNumberOfRows());

        bool ap_generated = false;
        //... and compare
        for (unsigned node = 0; node < reader_extended.GetNumberOfRows(); node++)
        {
            std::vector<double> voltage_first_cell_extended = reader_extended.GetVariableOverTime("V", node);
            std::vector<double> voltage_second_cell_extended = reader_extended.GetVariableOverTime("V_2", node);
            std::vector<double> phi_e_extended = reader_extended.GetVariableOverTime("Phi_e", node);

            std::vector<double> voltage_bidomain = reader_bidomain.GetVariableOverTime("V", node);
            std::vector<double> phi_e_bidomain = reader_bidomain.GetVariableOverTime("Phi_e", node);

            TS_ASSERT_EQUALS(voltage_first_cell_extended.size(), voltage_second_cell_extended.size());
            TS_ASSERT_EQUALS(voltage_second_cell_extended.size(), phi_e_extended.size());
            TS_ASSERT_EQUALS(phi_e_extended.size(), voltage_bidomain.size());
            TS_ASSERT_EQUALS(voltage_bidomain.size(), phi_e_bidomain.size());

            for ( unsigned index = 0; index < voltage_bidomain.size(); index ++)
            {
                //check that an AP was generated somewhere at some time
                if (voltage_bidomain[index] >0.0)
                {
                    ap_generated = true;
                }
                TS_ASSERT_DELTA(voltage_second_cell_extended[index] , voltage_bidomain[index], 2e-3);
                TS_ASSERT_LESS_THAN(voltage_first_cell_extended[index] , -83.5);//first unstimulated cell should be at rest
                TS_ASSERT_DELTA(phi_e_extended[index] , phi_e_bidomain[index], 2e-3);
            }
        }
        TS_ASSERT_EQUALS(ap_generated, true);
    }

    //////////////////////////////////////////////////////////////////////////////
    // Tests and simulations for the null space solution
    //////////////////////////////////////////////////////////////////////////////

    void RunExtendedSimulationWithNullBasis()
    {
        SetupParameters();

        StimulatedCellFactory stimulated_cell_factory;
        UnStimulatedCellFactory unstimulated_cell_factory;

        HeartConfig::Instance()->SetOutputDirectory("TestExtendedVs_Extended1dStimulateSecondCellNullBasis");
        HeartConfig::Instance()->SetOutputFilenamePrefix("extended1d");

        //stimulated cell factory passed in as second cell
        ExtendedBidomainProblem<1> extended_problem( &unstimulated_cell_factory,  &stimulated_cell_factory);

        //covering the case where the user sets extended parameters explicitly
        extended_problem.SetExtendedBidomainParameters(1400, 1400, 1400, 1.0 , 1.0, 0.0);

        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(0.0005));
        extended_problem.Initialise();
        extended_problem.Solve();
    }

    void RunBidomainNullBasis()
    {
        SetupParameters();

        HeartConfig::Instance()->SetOutputDirectory("TestExtendedVs_Bidomain1dNullBasis");
        HeartConfig::Instance()->SetOutputFilenamePrefix("normal1d");
        StimulatedCellFactoryBidomain bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );


        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }

    void TestCompareNullBasis()
    {
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);

        //first, run the extended bidomain simulation stimulating the second cell
        RunExtendedSimulationWithNullBasis();
        RunBidomainNullBasis();
        TS_ASSERT(strcmp(HeartConfig::Instance()->GetKSPSolver(), "cg")==0);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);

        //pick up the results...(note the bidomain simulatioon is the same as the test above for the first cell)
        Hdf5DataReader reader_extended("TestExtendedVs_Extended1dStimulateSecondCellNullBasis", "extended1d");
        Hdf5DataReader reader_bidomain("TestExtendedVs_Bidomain1dNullBasis", "normal1d");
        TS_ASSERT_EQUALS(reader_extended.GetNumberOfRows(), reader_bidomain.GetNumberOfRows());

        bool ap_generated = false;
        //... and compare
        for (unsigned node = 0; node < reader_extended.GetNumberOfRows(); node++)
        {
            std::vector<double> voltage_first_cell_extended = reader_extended.GetVariableOverTime("V", node);
            std::vector<double> voltage_second_cell_extended = reader_extended.GetVariableOverTime("V_2", node);
            std::vector<double> phi_e_extended = reader_extended.GetVariableOverTime("Phi_e", node);

            std::vector<double> voltage_bidomain = reader_bidomain.GetVariableOverTime("V", node);
            std::vector<double> phi_e_bidomain = reader_bidomain.GetVariableOverTime("Phi_e", node);

            TS_ASSERT_EQUALS(voltage_first_cell_extended.size(), voltage_second_cell_extended.size());
            TS_ASSERT_EQUALS(voltage_second_cell_extended.size(), phi_e_extended.size());
            TS_ASSERT_EQUALS(phi_e_extended.size(), voltage_bidomain.size());
            TS_ASSERT_EQUALS(voltage_bidomain.size(), phi_e_bidomain.size());

            for ( unsigned index = 0; index < voltage_bidomain.size(); index ++)
            {
                //check that an AP was generated somewhere at some time
                if (voltage_bidomain[index] >0.0)
                {
                    ap_generated = true;
                }
                TS_ASSERT_DELTA(voltage_second_cell_extended[index] , voltage_bidomain[index], 2e-3);
                TS_ASSERT_LESS_THAN(voltage_first_cell_extended[index] , -83.5);//first unstimulated cell should be at rest
                TS_ASSERT_DELTA(phi_e_extended[index] , phi_e_bidomain[index], 2e-3);
            }
        }
        TS_ASSERT_EQUALS(ap_generated, true);
    }

    //////////////////////////////////////////////////////////////////////////////
    // Tests and simulations for pinning a node
    //////////////////////////////////////////////////////////////////////////////

    void RunExtendedSimulationPinnedNode()
    {
        SetupParameters();

        StimulatedCellFactory stimulated_cell_factory;
        UnStimulatedCellFactory unstimulated_cell_factory;

        HeartConfig::Instance()->SetOutputDirectory("TestExtendedVs_Extended1dPinnedNode");
        HeartConfig::Instance()->SetOutputFilenamePrefix("extended1d");
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-7);
        //stimulated cell factory passed in as second cell
        ExtendedBidomainProblem<1> extended_problem( &unstimulated_cell_factory,  &stimulated_cell_factory);

        //covering the case where the user sets extended parameters explicitly
        extended_problem.SetExtendedBidomainParameters(1400, 1400, 1400, 1.0 , 1.0, 0.0);

        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(0.0005));
        //////////////////////////
        ///Pinning a node in the extended bidomain problem
        //////////////////////////
        std::vector<unsigned> pinned_nodes;
        pinned_nodes.push_back(100);
        extended_problem.SetFixedExtracellularPotentialNodes(pinned_nodes);

        extended_problem.Initialise();
        extended_problem.Solve();
    }

    void RunBidomainPinnedNode()
    {
        SetupParameters();

        HeartConfig::Instance()->SetOutputDirectory("TestExtendedVs_Bidomain1dPinnedNode");
        HeartConfig::Instance()->SetOutputFilenamePrefix("normal1d");
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-7);
        StimulatedCellFactoryBidomain bidomain_cell_factory;
        BidomainProblem<1> bidomain_problem( &bidomain_cell_factory );

        //////////////////////////
        ///Pinning a node in the bidomain problem (same as above)
        //////////////////////////
        std::vector<unsigned> pinned_nodes;
        pinned_nodes.push_back(100);
        bidomain_problem.SetFixedExtracellularPotentialNodes(pinned_nodes);

        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }

    void TestComparePinnedNode()
    {
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        //first, run the two simulations
        RunBidomainPinnedNode();
        RunExtendedSimulationPinnedNode();

        //pick up the results...
        Hdf5DataReader reader_extended("TestExtendedVs_Extended1dPinnedNode", "extended1d");
        Hdf5DataReader reader_bidomain("TestExtendedVs_Bidomain1dPinnedNode", "normal1d");
        TS_ASSERT_EQUALS(reader_extended.GetNumberOfRows(), reader_bidomain.GetNumberOfRows());

        bool ap_generated = false;
        //... and compare
        for (unsigned node = 0; node < reader_extended.GetNumberOfRows(); node++)
        {
            std::vector<double> voltage_first_cell_extended = reader_extended.GetVariableOverTime("V", node);
            std::vector<double> voltage_second_cell_extended = reader_extended.GetVariableOverTime("V_2", node);
            std::vector<double> phi_e_extended = reader_extended.GetVariableOverTime("Phi_e", node);

            std::vector<double> voltage_bidomain = reader_bidomain.GetVariableOverTime("V", node);
            std::vector<double> phi_e_bidomain = reader_bidomain.GetVariableOverTime("Phi_e", node);

            TS_ASSERT_EQUALS(voltage_first_cell_extended.size(), voltage_second_cell_extended.size());
            TS_ASSERT_EQUALS(voltage_second_cell_extended.size(), phi_e_extended.size());
            TS_ASSERT_EQUALS(phi_e_extended.size(), voltage_bidomain.size());
            TS_ASSERT_EQUALS(voltage_bidomain.size(), phi_e_bidomain.size());

            for ( unsigned index = 0; index < voltage_bidomain.size(); index ++)
            {
                //check that an AP was generated somewhere at some time
                if (voltage_bidomain[index] >0.0)
                {
                    ap_generated = true;
                }
                TS_ASSERT_DELTA(voltage_second_cell_extended[index] , voltage_bidomain[index], 2e-3);
                TS_ASSERT_LESS_THAN(voltage_first_cell_extended[index] , -83.5);//first unstimulated cell should be at rest
                TS_ASSERT_DELTA(phi_e_extended[index] , phi_e_bidomain[index], 2e-3);
            }
        }
        TS_ASSERT_EQUALS(ap_generated, true);
    }

    //////////////////////////////////////////////////////////////////////////////
    // Tests for heterogeneous Ggap
    //////////////////////////////////////////////////////////////////////////////

    void TestHeterogeneousGgap()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(65.0));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400.0);
        HeartConfig::Instance()->SetOdeTimeStep(0.01);
        HeartConfig::Instance()->SetCapacitance(1.0);
        //HeartConfig::Instance()->SetKSPSolver("gmres");
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        HeartConfig::Instance()->SetSimulationDuration(5.0);

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> unstimulated_cell_factory(0.0);
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> stimulated_cell_factory(-105000.0);

        ///////////////
        // Putting the second half of the cable with Ggap=0
        ///////////////
        std::vector<boost::shared_ptr<AbstractChasteRegion<1> > > heterogeneity_areas;
        std::vector<double> Ggap_values;
        //second half of the cable
        ChastePoint<1> cornerA(-0.1);
        ChastePoint<1> cornerB(0.1);
        boost::shared_ptr<ChasteCuboid<1> > p_cuboid_1(new ChasteCuboid<1>(cornerA, cornerB));
        heterogeneity_areas.push_back(p_cuboid_1);
        //...where cells are coupled
        Ggap_values.push_back(2.0);

        HeartConfig::Instance()->SetOutputDirectory("TestExtendedVs_HeterogeneousGgaps");
        HeartConfig::Instance()->SetOutputFilenamePrefix("extended1d");
        //stimulated cell factory passed in as second cell
        ExtendedBidomainProblem<1> extended_problem( &unstimulated_cell_factory,  &stimulated_cell_factory);
        extended_problem.SetMesh(&mesh);
        //Ggap is zero everywhere else
        extended_problem.SetExtendedBidomainParameters(1400, 1400, 1400, 1.0 , 1.0, 0.0);
        extended_problem.SetGgapHeterogeneities(heterogeneity_areas, Ggap_values);

        extended_problem.Initialise();
        extended_problem.Solve();

        /**
         * We test the output.  What we have in the test was that Ggap was 2.0 at the beginning
         * and 0 everywhere else.
         *
         * If it was 0 everywhere (i.e., the Ggap methods failed to apply the proper values)
         * an AP would be triggered in the second (stimulated cell), while the second would stay at rest (this was checked).
         * Here, therefore, we check that the first cell peaks at about -68 and (and it is not at rest)
         * the second (stimulated) at about -57.8 (much lower than it would be if an AP was triggered).
         */
        Hdf5DataReader reader_extended("TestExtendedVs_HeterogeneousGgaps", "extended1d");
        unsigned probe_node = 0u; //we probe the stimulated one

        const std::vector<unsigned>& r_permutation = mesh.rGetNodePermutation();
        if (!r_permutation.empty())
        {
            probe_node = r_permutation[probe_node];
        }

        std::vector<double> voltage_first_cell_extended = reader_extended.GetVariableOverTime("V", probe_node);
        std::vector<double> voltage_second_cell_extended = reader_extended.GetVariableOverTime("V_2", probe_node);
        TS_ASSERT_EQUALS(voltage_first_cell_extended.size(), voltage_second_cell_extended.size());

        double first_cell_peak = -1000000.0;
        double second_cell_peak = -1000000.0;
        for (unsigned i = 0; i < voltage_first_cell_extended.size(); i++)
        {
            if (voltage_first_cell_extended[i] >=first_cell_peak)
            {
                first_cell_peak = voltage_first_cell_extended[i];
            }
            if (voltage_second_cell_extended[i] >= second_cell_peak)
            {
                second_cell_peak = voltage_second_cell_extended[i];
            }
        }

        TS_ASSERT_DELTA(first_cell_peak, -67.9974, 2e-4);
        TS_ASSERT_DELTA(second_cell_peak, -57.7988, 2e-4);


    }

    /**
     * This test is just to try out (and cover) the WriteInfo method and the set  method of the extracellular stimulus
     * Also tests method for accessing the flag indicating whether the user suipplied an extracellular stimulus or not
     * */
    void TestSomeOtherMethods()
    {
        HeartEventHandler::Reset();
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.0005));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(65.0));//very high compared to intracellular
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        HeartConfig::Instance()->SetSimulationDuration(0.03);  //ms.

        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        UnStimulatedCellFactory first_cell_factory;
        UnStimulatedCellFactory second_cell_factory;
        ExtracellularStimulusFactory extra_factory;

        HeartConfig::Instance()->SetOutputDirectory("TestExtendedVs_Extended1dOtherMethods");
        HeartConfig::Instance()->SetOutputFilenamePrefix("extended1d");

        //stimulated cell factory passed in as first cell
        ExtendedBidomainProblem<1> extended_problem( &first_cell_factory , &second_cell_factory);
        extended_problem.SetMesh(&mesh);
        extended_problem.SetWriteInfo(true);
        extended_problem.SetNodeForAverageOfPhiZeroed(0u);
        extended_problem.SetExtracellularStimulusFactory(&extra_factory);
        extended_problem.Initialise();
        TS_ASSERT_EQUALS(extended_problem.GetExtendedBidomainTissue()->HasTheUserSuppliedExtracellularStimulus(), true);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 0u);
        extended_problem.Solve();
        // Averaged phi_e forces the solver to use GMRES
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNumWarnings(), 1u);
        TS_ASSERT_EQUALS(Warnings::Instance()->GetNextWarningMessage(),"Code has changed the KSP solver type from cg to gmres");
        Warnings::Instance()->QuietDestroy();
    }

    void TestExceptions()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetSimulationDuration(1.0);  //ms
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/1D_0_to_1_100_elements");
        HeartConfig::Instance()->SetOutputDirectory("TestExtendedVs_ExtendedExceptions");
        HeartConfig::Instance()->SetOutputFilenamePrefix("exceptions");

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory;
        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> bidomain_cell_factory_2;
        ExtendedBidomainProblem<1> extended_problem( &bidomain_cell_factory , &bidomain_cell_factory_2);

        extended_problem.Initialise();

        // Exception in PreSolveChecks()
        HeartConfig::Instance()->SetUseRelativeTolerance(1e-4);
        TS_ASSERT_THROWS_THIS( extended_problem.Solve()
                , "Bidomain external voltage is not bounded in this simulation - use KSP *absolute* tolerance");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-4);

        //Exception in SetupLinearSystem within MatrixBasedAssembler
        extended_problem.SetHasBath(true);
        TS_ASSERT_THROWS_THIS( extended_problem.Solve()
                ,"Bath simulations are not yet supported for extended bidomain problems");
        extended_problem.SetHasBath(false);
        TS_ASSERT_EQUALS(extended_problem.GetHasBath(),false);

        HeartEventHandler::Reset();
        // check throws if the fixed node num isn't valid
        std::vector<unsigned> pinned_nodes;
        pinned_nodes.push_back(1000);
        extended_problem.SetFixedExtracellularPotentialNodes(pinned_nodes);

        // Covering exception thrown in ExtendedBidomainProblem<DIM>::CreateAssembler()
        TS_ASSERT_THROWS_THIS(extended_problem.Solve(),"Fixed node number must be less than total number nodes");

        //cover an exception in SetGgapHeterogeneities
        std::vector<boost::shared_ptr<AbstractChasteRegion<1> > > heterogeneity_areas;
        std::vector<double> Ggap_values;

        ChastePoint<1> cornerA(-1);
        ChastePoint<1> cornerB(0.001);
        boost::shared_ptr<ChasteCuboid<1> > p_cuboid_1(new ChasteCuboid<1>(cornerA, cornerB));
        heterogeneity_areas.push_back(p_cuboid_1);
        Ggap_values.push_back(143.0);
        Ggap_values.push_back(91.0);//now Ggap_values has 2 members while heterogeneity_areas only 1.
        TS_ASSERT_THROWS_THIS(extended_problem.SetGgapHeterogeneities(heterogeneity_areas, Ggap_values),
            "Gap junction heterogeneity areas must be of the same number as the heterogeneity values");

        // Coverage of the exception in the solver itself
        BoundaryConditionsContainer<1,1,3> container;
        ExtendedBidomainSolver<1,1> bidomain_solver(false,
                                                       &extended_problem.rGetMesh(),
                                                       extended_problem.GetExtendedBidomainTissue(),
                                                       &container);

        TS_ASSERT_THROWS_THIS(bidomain_solver.SetRowForAverageOfPhiZeroed(4),
                "Row for applying the constraint 'Average of phi_e = zero' should be every 3 rows");
    }
};

#endif /*TESTExtendedBidomainPROBLEM_HPP_*/
