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


#ifndef TESTARCHIVINGEXTENDEDBIDOMAIN_HPP_
#define TESTARCHIVINGEXTENDEDBIDOMAIN_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <vector>

#include "UblasIncludes.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "HeartConfig.hpp"
#include "SimpleStimulus.hpp"
#include "ArchiveOpener.hpp"
#include "LuoRudy1991.hpp"
#include "ExtendedBidomainProblem.hpp"
#include "CompareHdf5ResultsFiles.hpp"

/**
 * Cell factories for the archiving tests.
 */
class StimulatedCellFactory: public AbstractCardiacCellFactory<2>
{
private:
        boost::shared_ptr<SimpleStimulus> mpStimulus;
        //static const double magnitude = -445000.0;
public:
        StimulatedCellFactory() : AbstractCardiacCellFactory<2>(),
            mpStimulus ( new SimpleStimulus(-445000.0, 1.0))/*amplitude, duration (ms)*/
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        CellLuoRudy1991FromCellML* first_cell;
        if ((x < 0.005) )
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

class UnStimulatedCellFactory: public AbstractCardiacCellFactory<2>
{
private:
        boost::shared_ptr<SimpleStimulus> mpStimulus;
        //static const double magnitude = 0.0;
public:
        UnStimulatedCellFactory() : AbstractCardiacCellFactory<2>(),
                    mpStimulus ( new SimpleStimulus(0.0, 1.0)) ///\todo Shouldn't this be a ZeroStimulus?
    {
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        CellLuoRudy1991FromCellML* second_cell = new CellLuoRudy1991FromCellML(mpSolver, mpStimulus);
        return second_cell;
    }
};

/**
* This should set a compatible stimulus with the mesh mesh/test/data/2D_0_to_1mm_400_elements
* It is intended to start after the 'save' (3 ms) so that we check it is unarchived properly
*/
class ExtracellularStimulusFactory: public AbstractStimulusFactory<2>
{

public:
    ExtracellularStimulusFactory() : AbstractStimulusFactory<2>()
    {
    }

    boost::shared_ptr<AbstractStimulusFunction> CreateStimulusForNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];

        boost::shared_ptr<SimpleStimulus> p_stimulus;
        if ((x <= 0.01))
        {
            p_stimulus.reset( new SimpleStimulus(25000, 0.5, 3.5) );
        }
        else if (x >= 0.09)
        {
            p_stimulus.reset( new SimpleStimulus(-25000, 0.5, 3.5) );
        }
        else
        {
            p_stimulus.reset( new SimpleStimulus(0.0, 1.0, 0.1) );
        }
        return p_stimulus;
    }
};

class TestArchivingExtendedBidomain : public CxxTest::TestSuite
{

public:

    /**
     * Test of archiving an extended bidomain problem with intracellular stimulus.
     *
     * The test does the following:
     * - Set up and run a simulation for 3 ms and save everything to an archive.
     * - Load the simulation, check that all the member variables are in the correct state and run additional 2 ms.
     * - Run a new full 5 ms (3 + 2) simulation with the same setup.
     * - Compare hdf5 the hdf5 output from the full and the checkpointed simulations and check they are the same.
     *
     * in this test, an heterogeneous pattern pf Ggap is set up to check for proper archiving
     */
    void TestArchivingProblemIntraStim()
    {
        FileFinder archive_dir("extended_bidomain_problem_archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "extended.arch";

        mProbeNode = 15u;//for comparing node traces

        // Run 3 ms and save
        Run2DSimulationSaveAfterThreemilliSecondsIntraStim(archive_dir, archive_file);
        // Resume and run remaining 2 ms
        LoadAndRunRemainingTwoMilliSecondsIntraStim(archive_dir, archive_file);

        // Now run, from scratch the same full extended bidomain simulation (5 ms)
        RunFull2DSimulationIntraStim();

        // check that full and checkpointed simulations generate the same hdf5 result
        TS_ASSERT(CompareFilesViaHdf5DataReader("Extended2DFull", "extended2d", true,
                                                "Extended2DArchived", "extended2d", true,
                                                5e-9));
    }


    /**
     * Test of archiving an extended bidomain problem with extracellular stimulus.
     * This is the same as above but we check that the extracellular stimulus is archived properly
     */
    void TestArchivingProblemExtraStim()
    {
        FileFinder archive_dir("extended_bidomain_problem_archive_extrastim", RelativeTo::ChasteTestOutput);
        std::string archive_file = "extended.arch";

        mProbeNode = 15u;//for comparing node traces

        // Run 3 ms and save
        Run2DSimulationSaveAfterThreemilliSecondsExtraStim(archive_dir, archive_file);
        // Resume and run remaining 2 ms
        LoadAndRunRemainingTwoMilliSecondsExtraStim(archive_dir, archive_file);

        // Now run, from scratch the same full extended bidomain simulation (5 ms)
        RunFull2DSimulationExtraStim();

        // check that full and checkpointed simulations generate the same hdf5 result
        TS_ASSERT(CompareFilesViaHdf5DataReader("Extended2DFullExtraStim", "extended2d", true,
                                                "Extended2DArchivedExtraStim", "extended2d", true,
                                                5e-9));
    }
private:

    unsigned mProbeNode; //probe node in the original mesh numbering

    void SetupParameters()
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.5, 0.5));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(5.0, 5.0));
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(1400.0);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.01, 0.01, 0.1);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetKSPSolver("gmres");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-5);
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        HeartConfig::Instance()->SetMeshFileName("mesh/test/data/2D_0_to_1mm_400_elements");
        //Dumb partitioning is needed because we unarchive only with dumb partitioning (see #1199).
        //This does not affect the hdf5 file, but it does affect the numbering of the node for the traces (which we test here).
        HeartConfig::Instance()->SetMeshPartitioning("dumb");
    }

////////////////////////////////////////////////////////////////////////////////////////////
//// Functions for archiving problems with intracellular stimulus
////////////////////////////////////////////////////////////////////////////////////////////

    void RunFull2DSimulationIntraStim()
    {
        SetupParameters();

        HeartConfig::Instance()->SetSimulationDuration(5.0);

        StimulatedCellFactory stimulated_cell_factory;
        UnStimulatedCellFactory unstimulated_cell_factory;

        HeartConfig::Instance()->SetOutputDirectory("Extended2DFull");
        HeartConfig::Instance()->SetOutputFilenamePrefix("extended2d");

        //stimulated cell factory passed in as second cell
        ExtendedBidomainProblem<2> extended_problem( &unstimulated_cell_factory,  &stimulated_cell_factory);

        //change some values of parameters
        extended_problem.SetExtendedBidomainParameters(1400, 1500, 1600, 1.0 , 1.0, 0.0);
        //and the case where second cell has different intracellular conductivities.
        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(1.5,1.5));

        //also put in non-default Ggap heterogeneity pattern to check for archiving
        std::vector<boost::shared_ptr<AbstractChasteRegion<2> > > heterogeneity_areas;
        std::vector<double> Ggap_values;
        ChastePoint<2> cornerA(-1, -1);
        ChastePoint<2> cornerB(0.05, 0.05);
        boost::shared_ptr<ChasteCuboid<2> > p_cuboid_1(new ChasteCuboid<2>(cornerA, cornerB));
        heterogeneity_areas.push_back(p_cuboid_1);
        Ggap_values.push_back(0.1);
        extended_problem.SetGgapHeterogeneities(heterogeneity_areas, Ggap_values);

        extended_problem.Initialise();
        extended_problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }

    void Run2DSimulationSaveAfterThreemilliSecondsIntraStim(FileFinder archive_dir, std::string archive_file)
    {
        SetupParameters();

        HeartConfig::Instance()->SetSimulationDuration(3.0);

        StimulatedCellFactory stimulated_cell_factory;
        UnStimulatedCellFactory unstimulated_cell_factory;

        HeartConfig::Instance()->SetOutputDirectory("Extended2DArchived");
        HeartConfig::Instance()->SetOutputFilenamePrefix("extended2d");

        //stimulated cell factory passed in as second cell
        ExtendedBidomainProblem<2> extended_problem( &unstimulated_cell_factory,  &stimulated_cell_factory);

        //change some values of parameters
        extended_problem.SetExtendedBidomainParameters(1400, 1500, 1600, 1.0 , 1.0, 0.0);
        //and the case where second cell has different intracellular conductivities.
        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(1.5,1.5));

        //also put in non-default Ggap heterogeneity pattern to check for archiving
        std::vector<boost::shared_ptr<AbstractChasteRegion<2> > > heterogeneity_areas;
        std::vector<double> Ggap_values;
        ChastePoint<2> cornerA(-1, -1);
        ChastePoint<2> cornerB(0.05, 0.05);
        boost::shared_ptr<ChasteCuboid<2> > p_cuboid_1(new ChasteCuboid<2>(cornerA, cornerB));
        heterogeneity_areas.push_back(p_cuboid_1);
        Ggap_values.push_back(0.1);
        extended_problem.SetGgapHeterogeneities(heterogeneity_areas, Ggap_values);

        extended_problem.Initialise();
        extended_problem.Solve();

        ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
        boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

        AbstractCardiacProblem<2,2,3>* const p_extended_problem = &extended_problem;
        (*p_arch) & p_extended_problem; //archive
    }

    void LoadAndRunRemainingTwoMilliSecondsIntraStim(FileFinder archive_dir, std::string archive_file)
    {
        ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
        boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

        AbstractCardiacProblem<2,2,3> *p_problem;
        (*p_arch) >> p_problem;//load

        //see if the heartconfig is saved properly
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 3.0);

        //check specific variables of the extended problem (need to do dynamic cast)
        ExtendedBidomainProblem<2>* p_extended_problem =  dynamic_cast<ExtendedBidomainProblem<2>*>(p_problem);
        TS_ASSERT_EQUALS(p_extended_problem->GetHasBath(), false);
        TS_ASSERT_EQUALS(p_extended_problem->mUserSpecifiedSecondCellConductivities, true);
        TS_ASSERT_EQUALS(p_extended_problem->mAmFirstCell, 1400);
        TS_ASSERT_EQUALS(p_extended_problem->mAmSecondCell, 1500);
        TS_ASSERT_EQUALS(p_extended_problem->mAmGap, 1600);
        TS_ASSERT_EQUALS(p_extended_problem->mCmFirstCell, 1.0);
        TS_ASSERT_EQUALS(p_extended_problem->mCmSecondCell, 1.0);
        TS_ASSERT_EQUALS(p_extended_problem->mUserHasSetBidomainValuesExplicitly, true);
        TS_ASSERT_EQUALS(p_extended_problem->mFixedExtracellularPotentialNodes.size(), 0u);//default value here
        TS_ASSERT_EQUALS(p_extended_problem->mGgapHeterogeneityRegions.size(), 1u);
        ChastePoint<2> contained_point(0.0, 0.0);
        ChastePoint<2> not_contained_point(25.0, 0.0);
        TS_ASSERT_EQUALS(p_extended_problem->mGgapHeterogeneityRegions[0]->DoesContain(contained_point), true);
        TS_ASSERT_EQUALS(p_extended_problem->mGgapHeterogeneityRegions[0]->DoesContain(not_contained_point), false);
        TS_ASSERT_EQUALS(p_extended_problem->mGgapHeterogenousValues.size(), 1u);
        TS_ASSERT_EQUALS(p_extended_problem->mGgapHeterogenousValues[0], 0.1);
        TS_ASSERT_EQUALS(p_extended_problem->mIntracellularConductivitiesSecondCell(0), 1.5);
        TS_ASSERT_EQUALS(p_extended_problem->mIntracellularConductivitiesSecondCell(1), 1.5);
        TS_ASSERT_EQUALS(p_extended_problem->mVariablesIDs.size(), 3u);
        TS_ASSERT_EQUALS(p_extended_problem->mRowForAverageOfPhiZeroed, INT_MAX);//default value in this case
        TS_ASSERT_EQUALS(p_extended_problem->mApplyAveragePhieZeroConstraintAfterSolving, false);
        TS_ASSERT_EQUALS(p_extended_problem->mUserSuppliedExtracellularStimulus, false);

        //check that the tissue object has distributed vectors of the right size
        unsigned size_first_cell, size_second_cell, size_extrastim, global_size_first_cell, global_size_second_cell, global_size_extrastim;
        for (unsigned proc = 0; proc < PetscTools::GetNumProcs(); proc++)
        {
            size_first_cell = p_problem->GetTissue()->rGetCellsDistributed().size();
            size_second_cell = p_extended_problem->GetExtendedBidomainTissue()->rGetSecondCellsDistributed().size();
            size_extrastim = p_extended_problem->GetExtendedBidomainTissue()->rGetExtracellularStimulusDistributed().size();
        }
        int mpi_ret = MPI_Allreduce(&size_first_cell, &global_size_first_cell, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        UNUSED_OPT(mpi_ret);
        assert(mpi_ret == MPI_SUCCESS);
        mpi_ret = MPI_Allreduce(&size_second_cell, &global_size_second_cell, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        assert(mpi_ret == MPI_SUCCESS);
        mpi_ret = MPI_Allreduce(&size_extrastim, &global_size_extrastim, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        assert(mpi_ret == MPI_SUCCESS);

        unsigned number_of_nodes = p_problem->GetTissue()->pGetMesh()->GetNumNodes();
        TS_ASSERT_EQUALS(global_size_first_cell, number_of_nodes);
        TS_ASSERT_EQUALS(global_size_second_cell, number_of_nodes);
        TS_ASSERT_EQUALS(global_size_extrastim, number_of_nodes);

        //now solve for the remining 2 ms
        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        p_problem->Solve();
        delete p_problem;
    }


////////////////////////////////////////////////////////////////////////////////////////////
//// Functions for archiving problems with extracellular stimulus
////////////////////////////////////////////////////////////////////////////////////////////


    void RunFull2DSimulationExtraStim()
    {
        SetupParameters();

        HeartConfig::Instance()->SetSimulationDuration(5.0);

        UnStimulatedCellFactory unstimulated_cell_factory_1;
        UnStimulatedCellFactory unstimulated_cell_factory_2;
        ExtracellularStimulusFactory extra_stim;

        HeartConfig::Instance()->SetOutputDirectory("Extended2DFullExtraStim");
        HeartConfig::Instance()->SetOutputFilenamePrefix("extended2d");

        //stimulated cell factory passed in as second cell
        ExtendedBidomainProblem<2> extended_problem( &unstimulated_cell_factory_1,  &unstimulated_cell_factory_2);
        extended_problem.SetExtracellularStimulusFactory(&extra_stim);

        //change some values of parameters
        extended_problem.SetExtendedBidomainParameters(1800, 1900, 2100, 1.0 , 1.0, 0.0);
        //and the case where second cell has different intracellular conductivities.
        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(1.5,1.5));

        extended_problem.Initialise();
        extended_problem.Solve();

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
    void Run2DSimulationSaveAfterThreemilliSecondsExtraStim(FileFinder archive_dir, std::string archive_file)
    {
        SetupParameters();

        HeartConfig::Instance()->SetSimulationDuration(3.0);

        UnStimulatedCellFactory unstimulated_cell_factory_1;
        UnStimulatedCellFactory unstimulated_cell_factory_2;
        ExtracellularStimulusFactory extra_stim;

        HeartConfig::Instance()->SetOutputDirectory("Extended2DArchivedExtraStim");
        HeartConfig::Instance()->SetOutputFilenamePrefix("extended2d");

        ExtendedBidomainProblem<2> extended_problem( &unstimulated_cell_factory_1,  &unstimulated_cell_factory_2);
        extended_problem.SetExtracellularStimulusFactory(&extra_stim);

        //change some values of parameters
        extended_problem.SetExtendedBidomainParameters(1800, 1900, 2100, 1.0 , 1.0, 0.0);
        //and the case where second cell has different intracellular conductivities.
        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(1.5,1.5));

        extended_problem.Initialise();
        extended_problem.Solve();

        ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
        boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

        AbstractCardiacProblem<2,2,3>* const p_extended_problem = &extended_problem;
        (*p_arch) & p_extended_problem; //archive
    }

    void LoadAndRunRemainingTwoMilliSecondsExtraStim(FileFinder archive_dir, std::string archive_file)
    {
        ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
        boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

        AbstractCardiacProblem<2,2,3> *p_problem;
        (*p_arch) >> p_problem;//load

        //see if the heartconfig is saved properly
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(), 3.0);

        //check specific variables of the extended problem (need to do dynamic cast)
        ExtendedBidomainProblem<2>* p_extended_problem =  dynamic_cast<ExtendedBidomainProblem<2>*>(p_problem);
        TS_ASSERT_EQUALS(p_extended_problem->GetHasBath(), false);
        TS_ASSERT_EQUALS(p_extended_problem->mUserSpecifiedSecondCellConductivities, true);
        TS_ASSERT_EQUALS(p_extended_problem->mAmFirstCell, 1800);
        TS_ASSERT_EQUALS(p_extended_problem->mAmSecondCell, 1900);
        TS_ASSERT_EQUALS(p_extended_problem->mAmGap, 2100);
        TS_ASSERT_EQUALS(p_extended_problem->mCmFirstCell, 1.0);
        TS_ASSERT_EQUALS(p_extended_problem->mCmSecondCell, 1.0);
        TS_ASSERT_EQUALS(p_extended_problem->mUserHasSetBidomainValuesExplicitly, true);
        TS_ASSERT_EQUALS(p_extended_problem->mFixedExtracellularPotentialNodes.size(), 0u);//default value here
        TS_ASSERT_EQUALS(p_extended_problem->mIntracellularConductivitiesSecondCell(0), 1.5);
        TS_ASSERT_EQUALS(p_extended_problem->mIntracellularConductivitiesSecondCell(1), 1.5);
        TS_ASSERT_EQUALS(p_extended_problem->mVariablesIDs.size(), 3u);
        TS_ASSERT_EQUALS(p_extended_problem->mRowForAverageOfPhiZeroed, INT_MAX);//default value in this case
        TS_ASSERT_EQUALS(p_extended_problem->mApplyAveragePhieZeroConstraintAfterSolving, false);
        TS_ASSERT_EQUALS(p_extended_problem->mUserSuppliedExtracellularStimulus, true);

        //chec that the tissue object has distributed vectors of the right size
        unsigned size_first_cell, size_second_cell, size_extrastim, global_size_first_cell, global_size_second_cell, global_size_extrastim;
        for (unsigned proc = 0; proc < PetscTools::GetNumProcs(); proc++)
        {
            size_first_cell = p_problem->GetTissue()->rGetCellsDistributed().size();
            size_second_cell = p_extended_problem->GetExtendedBidomainTissue()->rGetSecondCellsDistributed().size();
            size_extrastim = p_extended_problem->GetExtendedBidomainTissue()->rGetExtracellularStimulusDistributed().size();
        }
        int mpi_ret = MPI_Allreduce(&size_first_cell, &global_size_first_cell, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        UNUSED_OPT(mpi_ret);
        assert(mpi_ret == MPI_SUCCESS);
        mpi_ret = MPI_Allreduce(&size_second_cell, &global_size_second_cell, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        assert(mpi_ret == MPI_SUCCESS);
        mpi_ret = MPI_Allreduce(&size_extrastim, &global_size_extrastim, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        assert(mpi_ret == MPI_SUCCESS);

        unsigned number_of_nodes = p_problem->GetTissue()->pGetMesh()->GetNumNodes();
        TS_ASSERT_EQUALS(global_size_first_cell, number_of_nodes);
        TS_ASSERT_EQUALS(global_size_second_cell, number_of_nodes);
        TS_ASSERT_EQUALS(global_size_extrastim, number_of_nodes);

        //now solve for the remining 2 ms
        HeartConfig::Instance()->SetSimulationDuration(5.0); //ms
        p_problem->Solve();
        delete p_problem;
    }


};
#endif //TESTARCHIVINGEXTENDEDBIDOMAIN_HPP_
