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

#ifndef TESTODEBASEDCELLCYCLEMODELS_HPP_
#define TESTODEBASEDCELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <boost/shared_ptr.hpp>

#include "Alarcon2004OxygenBasedCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * This class contains tests for methods on classes
 * inheriting from AbstractOdeBasedCellCycleModel.
 */
class TestOdeBasedCellCycleModels : public AbstractCellBasedTestSuite
{
public:

    void TestTysonNovakCellCycleModel()
    {
        // Set up
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_timesteps = 50;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, num_timesteps);

        double standard_divide_time = 75.19/60.0;

        // Test TysonNovakCellCycleModel methods for a healthy cell
        TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel;
        p_cell_model->SetBirthTime(p_simulation_time->GetTime());

        TS_ASSERT_EQUALS(p_cell_model->CanCellTerminallyDifferentiate(), false);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();
        p_cell_model->SetDt(0.1/60.0);

        /*
         * For coverage, we create another cell-cycle model that is identical except that we
         * manually pass in an ODE solver. In this case, our ODE solver (BackwardEulerIvpOdeSolver)
         * is the same type as the solver used by the cell-cycle model if no solver is provided
         * (unless CVODE is used), so our results should be identical.
         */
        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver> >
            p_solver(CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver>::Instance());
        p_solver->SetSizeOfOdeSystem(6);
        p_solver->Initialise();

        TysonNovakCellCycleModel* p_other_cell_model = new TysonNovakCellCycleModel(p_solver);

        // Coverage of GetOdeSolver()
        boost::shared_ptr<AbstractCellCycleModelOdeSolver> p_solver_from_model = p_other_cell_model->GetOdeSolver();
        TS_ASSERT_EQUALS(p_solver_from_model->GetSizeOfOdeSystem(), 6u);

        p_other_cell_model->SetBirthTime(p_simulation_time->GetTime());

        // Timestep for non-adaptive solvers defaults to 0.0001
        TS_ASSERT_EQUALS(p_other_cell_model->GetDt(), 0.0001);
        p_other_cell_model->SetDt(0.1/60.0);

        CellPtr p_other_cell(new Cell(p_healthy_state, p_other_cell_model));
        p_other_cell->SetCellProliferativeType(p_stem_type);
        p_other_cell->InitialiseCellCycleModel();

        // Test the cell is ready to divide at the right time
        for (unsigned i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool result = p_cell_model->ReadyToDivide();
            bool other_result = p_other_cell_model->ReadyToDivide();

            if (time > standard_divide_time)
            {
                TS_ASSERT_EQUALS(result, true);
                TS_ASSERT_EQUALS(other_result, true);
            }
            else
            {
                TS_ASSERT_EQUALS(result, false);
                TS_ASSERT_EQUALS(other_result, false);
            }
        }

        // Test ODE solution
        std::vector<double> proteins = p_cell_model->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins.size(), 6u);
        TS_ASSERT_DELTA(proteins[0], 0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(proteins[1], 0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(proteins[2], 1.54216806705641, 1e-1);
        TS_ASSERT_DELTA(proteins[3], 1.40562614481544, 2e-2);
        TS_ASSERT_DELTA(proteins[4], 0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(proteins[5], 0.95328206604519, 2e-2);

        std::vector<double> other_proteins = p_other_cell_model->GetProteinConcentrations();
        TS_ASSERT_EQUALS(other_proteins.size(), 6u);
        TS_ASSERT_DELTA(other_proteins[0], 0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(other_proteins[1], 0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(other_proteins[2], 1.54216806705641, 1e-1);
        TS_ASSERT_DELTA(other_proteins[3], 1.40562614481544, 2e-2);
        TS_ASSERT_DELTA(other_proteins[4], 0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(other_proteins[5], 0.95328206604519, 2e-2);

        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),true);

        // For coverage, we also test TysonNovakCellCycleModel methods for a mutant cell
        p_cell_model->ResetForDivision();

        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),false);

        TysonNovakCellCycleModel* p_cell_model2 = static_cast<TysonNovakCellCycleModel*> (p_cell_model->CreateCellCycleModel());

        MAKE_PTR(ApcOneHitCellMutationState, p_mutation);
        CellPtr p_stem_cell_2(new Cell(p_mutation, p_cell_model2));

        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(),false);

        p_stem_cell_2->SetCellProliferativeType(p_stem_type);

        TS_ASSERT_EQUALS(p_cell_model->ReadyToDivide(),false);
        TS_ASSERT_EQUALS(p_cell_model2->ReadyToDivide(),false);

        // Test the cell is ready to divide at the right time
        for (unsigned i=0; i<num_timesteps/2; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            double time = p_simulation_time->GetTime();

            bool result = p_cell_model->ReadyToDivide();
            bool result2 = p_cell_model2->ReadyToDivide();

            if (time > 2.0* standard_divide_time)
            {
                TS_ASSERT_EQUALS(result, true);
                TS_ASSERT_EQUALS(result2, true);
            }
            else
            {
                TS_ASSERT_EQUALS(result, false);
                TS_ASSERT_EQUALS(result2, false);
            }
        }

        // Test ODE solution
        proteins = p_cell_model->GetProteinConcentrations();
        TS_ASSERT_EQUALS(proteins.size(), 6u);
        TS_ASSERT_DELTA(proteins[0], 0.10000000000000, 1e-2);
        TS_ASSERT_DELTA(proteins[1], 0.98913684535843, 1e-2);
        TS_ASSERT_DELTA(proteins[2], 1.54216806705641, 1e-1);
        TS_ASSERT_DELTA(proteins[3], 1.40562614481544, 1e-1);
        TS_ASSERT_DELTA(proteins[4], 0.67083371879876, 1e-2);
        TS_ASSERT_DELTA(proteins[5], 0.9662, 1e-2);

        // Coverage of AbstractOdeBasedCellCycleModel::SetProteinConcentrationsForTestsOnly()
        std::vector<double> test_results(6);
        for (unsigned i=0; i<6; i++)
        {
            test_results[i] = (double)i;
        }
        p_cell_model->SetProteinConcentrationsForTestsOnly(1.0, test_results);
        proteins = p_cell_model->GetProteinConcentrations();

        for (unsigned i=0; i<6; i++)
        {
            TS_ASSERT_DELTA(proteins[i], test_results[i], 1e-6);
        }
    }

    /**
     * Test for Tyson & Novak self-cycling cells without having their
     * initial conditions reset. When using CVODE, the cell-cycle model
     * resets itself by halving the mass of the cell.
     * When not using CVODE, the cell-cycle model resets its initial
     * conditions, since the oscillatory solution computed using the Chaste
     * ODE solver is not stable.
     */
    void TestTysonNovakCellCycleModelSolver()
    {
        // Set up simulation time
        unsigned num_timesteps = 100000;
        double standard_divide_time = 75.19/60.0;
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.1*standard_divide_time, num_timesteps);

        // Create cell-cycle model and associated cell
        TysonNovakCellCycleModel* p_repeating_cell_model = new TysonNovakCellCycleModel;

        MAKE_PTR(ApcOneHitCellMutationState, p_mutation);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        CellPtr p_tyson_novak_cell(new Cell(p_mutation, p_repeating_cell_model));
        p_tyson_novak_cell->SetCellProliferativeType(p_stem_type);
        p_tyson_novak_cell->InitialiseCellCycleModel();

        // Run through the cell-cycle model for a certain duration
        // and test how many times it has stopped for division
        unsigned num_divisions = 0;
        for (unsigned i=0; i<num_timesteps; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            bool result = p_repeating_cell_model->ReadyToDivide();
//                std::vector<double> proteins = p_repeating_cell_model->GetProteinConcentrations();
//                out << SimulationTime::Instance()->GetTime() << "\t";
//                for (unsigned j=0; j<proteins.size(); j++)
//                {
//                    out << proteins[j] << "\t";
//                }
//                out << "\n" << std::flush;

            if (result)
            {
                p_repeating_cell_model->ResetForDivision();
                p_repeating_cell_model->SetBirthTime(SimulationTime::Instance()->GetTime());
                num_divisions++;
            }
        }
        TS_ASSERT_LESS_THAN(num_divisions, 102u);
        TS_ASSERT_LESS_THAN(99u, num_divisions);
//            out.close();
        /*
         * Matlab code for plotting the output commented above:
         * cdchaste
         * data = load('TN_output.txt');
         * figure
         * for i=1:6
         *   subplot(3,2,i)
         *   plot(data(:,1),data(:,1+i))
         * end
         */
    }

    void TestAlarcon2004OxygenBasedCellCycleModel()
    {
        // Set up SimulationTime
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30.0, 3);

        // Set up oxygen_concentration
        double oxygen_concentration = 1.0;

         // Create cell-cycle models
        Alarcon2004OxygenBasedCellCycleModel* p_model_1d = new Alarcon2004OxygenBasedCellCycleModel();
        p_model_1d->SetDimension(1);

        Alarcon2004OxygenBasedCellCycleModel* p_model_2d = new Alarcon2004OxygenBasedCellCycleModel();
        p_model_2d->SetDimension(2);

        Alarcon2004OxygenBasedCellCycleModel* p_model_3d = new Alarcon2004OxygenBasedCellCycleModel();
        p_model_3d->SetDimension(3);

        // Create cells
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        CellPtr p_cell_1d(new Cell(p_state, p_model_1d));
        p_cell_1d->SetCellProliferativeType(p_stem_type);
        p_cell_1d->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell_1d->InitialiseCellCycleModel();

        CellPtr p_cell_2d(new Cell(p_state, p_model_2d));
        p_cell_2d->SetCellProliferativeType(p_stem_type);
        p_cell_2d->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell_2d->InitialiseCellCycleModel();

        CellPtr p_cell_3d(new Cell(p_state, p_model_3d));
        p_cell_3d->SetCellProliferativeType(p_stem_type);
        p_cell_3d->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell_3d->InitialiseCellCycleModel();

        // For coverage, we create another cell-cycle model that is identical to p_model_2d except for the ODE solver
        boost::shared_ptr<CellCycleModelOdeSolver<Alarcon2004OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver> >
            p_solver(CellCycleModelOdeSolver<Alarcon2004OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver>::Instance());
        p_solver->Initialise();

        Alarcon2004OxygenBasedCellCycleModel* p_other_model_2d = new Alarcon2004OxygenBasedCellCycleModel(p_solver);
        p_other_model_2d->SetDimension(2);

        CellPtr p_other_cell_2d(new Cell(p_state, p_other_model_2d));
        p_other_cell_2d->SetCellProliferativeType(p_stem_type);
        p_other_cell_2d->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_other_cell_2d->InitialiseCellCycleModel();

        // Check oxygen concentration is correct in cell-cycle model
        TS_ASSERT_DELTA(p_model_2d->GetProteinConcentrations()[5], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);

        TS_ASSERT_DELTA(p_other_model_2d->GetProteinConcentrations()[5], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_other_model_2d->ReadyToDivide(), false);

        // Divide the cells
        Alarcon2004OxygenBasedCellCycleModel* p_model_1d_2 = static_cast<Alarcon2004OxygenBasedCellCycleModel*> (p_model_1d->CreateCellCycleModel());
        CellPtr p_cell_1d_2(new Cell(p_state, p_model_1d_2));
        p_cell_1d_2->SetCellProliferativeType(p_stem_type);
        p_cell_1d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);

        Alarcon2004OxygenBasedCellCycleModel* p_model_2d_2 = static_cast<Alarcon2004OxygenBasedCellCycleModel*> (p_model_2d->CreateCellCycleModel());
        CellPtr p_cell_2d_2(new Cell(p_state, p_model_2d_2));
        p_cell_2d_2->SetCellProliferativeType(p_stem_type);
        p_cell_2d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);

        Alarcon2004OxygenBasedCellCycleModel* p_model_3d_2 = static_cast<Alarcon2004OxygenBasedCellCycleModel*> (p_model_3d->CreateCellCycleModel());
        CellPtr p_cell_3d_2(new Cell(p_state, p_model_3d_2));
        p_cell_3d_2->SetCellProliferativeType(p_stem_type);
        p_cell_3d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_model_1d->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_other_model_2d->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_model_3d->ReadyToDivide(), false);

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_model_1d->ReadyToDivide(), true)
        TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), true)
        TS_ASSERT_EQUALS(p_other_model_2d->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_model_3d->ReadyToDivide(), true);

        TS_ASSERT_THROWS_NOTHING(p_model_2d->ResetForDivision());

        // For coverage, create a 1D model
        Alarcon2004OxygenBasedCellCycleModel* p_cell_model3 = new Alarcon2004OxygenBasedCellCycleModel();
        p_cell_model3->SetDimension(1);

        CellPtr p_cell3(new Cell(p_state, p_cell_model3));
        p_cell3->SetCellProliferativeType(p_stem_type);
        p_cell3->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell3->InitialiseCellCycleModel();

        TS_ASSERT_DELTA(p_cell_model3->GetProteinConcentrations()[5], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
    }

    void TestArchiveTysonNovakCellCycleModels()
    {
        // Set up
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "TysonNovakCellCycleModel.arch";

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new TysonNovakCellCycleModel;

            p_model->SetDimension(3);
            p_model->SetBirthTime(-1.5);
            static_cast<TysonNovakCellCycleModel*>(p_model)->SetDt(0.085);

            // We must create a cell to be able to initialise the cell cycle model's ODE system
            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            CellPtr p_cell(new Cell(p_healthy_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->InitialiseCellCycleModel();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            // Note that here, deletion of the cell-cycle model is handled by the cell destructor
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_EQUALS(p_model2->GetDimension(), 3u);
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.5, 1e-12);
            TS_ASSERT_DELTA(static_cast<TysonNovakCellCycleModel*>(p_model2)->GetDt(), 0.085, 1e-3);

            TysonNovakCellCycleModel* p_static_cast_model =
                static_cast<TysonNovakCellCycleModel*>(p_model2);

            TysonNovak2001OdeSystem* p_ode_system =
                static_cast<TysonNovak2001OdeSystem*>(p_static_cast_model->GetOdeSystem());

            TS_ASSERT(p_ode_system != NULL);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveAlarcon2004OxygenBasedCellCycleModels()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "Alarcon2004OxygenBasedCellCycleModel.arch";

        // Set up oxygen_concentration
        double oxygen_concentration = 1.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new Alarcon2004OxygenBasedCellCycleModel;

            p_model->SetDimension(1);
            p_model->SetBirthTime(-1.5);
            static_cast<Alarcon2004OxygenBasedCellCycleModel*>(p_model)->SetDt(0.085);

            // We must create a cell to be able to initialise the cell cycle model's ODE system
            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            CellPtr p_cell(new Cell(p_healthy_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", oxygen_concentration);
            p_cell->InitialiseCellCycleModel();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            // Note that here, deletion of the cell-cycle model is handled by the cell destructor
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_EQUALS(p_model2->GetDimension(), 1u);
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.5, 1e-12);
            TS_ASSERT_DELTA(static_cast<Alarcon2004OxygenBasedCellCycleModel*>(p_model2)->GetDt(), 0.085, 1e-3);

            Alarcon2004OxygenBasedCellCycleModel* p_static_cast_model =
                static_cast<Alarcon2004OxygenBasedCellCycleModel*>(p_model2);

            Alarcon2004OxygenBasedCellCycleOdeSystem* p_ode_system =
                static_cast<Alarcon2004OxygenBasedCellCycleOdeSystem*>(p_static_cast_model->GetOdeSystem());

            TS_ASSERT(p_ode_system != NULL);
            TS_ASSERT_EQUALS(p_ode_system->IsLabelled(), false);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestCellCycleModelOutputParameters()
    {
        std::string output_directory = "TestCellCycleModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with Alarcon2004OxygenBasedCellCycleModel
        Alarcon2004OxygenBasedCellCycleModel alarcon_oxygen_based_cell_cycle_model;
        TS_ASSERT_EQUALS(alarcon_oxygen_based_cell_cycle_model.GetIdentifier(), "Alarcon2004OxygenBasedCellCycleModel");

        out_stream alarcon_oxygen_based_parameter_file = output_file_handler.OpenOutputFile("alarcon_oxygen_based_results.parameters");
        alarcon_oxygen_based_cell_cycle_model.OutputCellCycleModelParameters(alarcon_oxygen_based_parameter_file);
        alarcon_oxygen_based_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("alarcon_oxygen_based_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/alarcon_oxygen_based_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }

        // Test with TysonNovakCellCycleModel
        TysonNovakCellCycleModel tyson_novak_based_cell_cycle_model;
        TS_ASSERT_EQUALS(tyson_novak_based_cell_cycle_model.GetIdentifier(), "TysonNovakCellCycleModel");

        out_stream tyson_novak_based_parameter_file = output_file_handler.OpenOutputFile("tyson_novak_based_results.parameters");
        tyson_novak_based_cell_cycle_model.OutputCellCycleModelParameters(tyson_novak_based_parameter_file);
        tyson_novak_based_parameter_file->close();

        {
            FileFinder generated_file = output_file_handler.FindFile("tyson_novak_based_results.parameters");
            FileFinder reference_file("cell_based/test/data/TestCellCycleModels/tyson_novak_based_results.parameters",
                                      RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated_file,reference_file);
            TS_ASSERT(comparer.CompareFiles());
        }
    }
};

#endif /*TESTODEBASEDCELLCYCLEMODELS_HPP_*/
