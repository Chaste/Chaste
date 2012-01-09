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

#include "CellwiseData.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

/**
 * This class contains tests for methods on classes
 * inheriting from AbstractOdeBasedCellCycleModel.
 */
class TestOdeBasedCellCycleModels : public AbstractCellBasedTestSuite
{
public:

    void TestTysonNovakCellCycleModel() throw(Exception)
    {
        // Set up
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_timesteps = 50;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, num_timesteps);

        double standard_divide_time = 75.19/60.0;

        // Test TysonNovakCellCycleModel methods for a healthy cell
        TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel;
        p_cell_model->SetBirthTime(p_simulation_time->GetTime());
        p_cell_model->SetCellProliferativeType(STEM);

        TS_ASSERT_EQUALS(p_cell_model->CanCellTerminallyDifferentiate(), false);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->InitialiseCellCycleModel();

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
        p_other_cell_model->SetBirthTime(p_simulation_time->GetTime());
        p_other_cell_model->SetCellProliferativeType(STEM);
        // Timestep for non-adaptive solvers defaults to 0.0001
        TS_ASSERT_EQUALS(p_other_cell_model->GetDt(), 0.0001);
        p_other_cell_model->SetDt(0.1/60.0);

        CellPtr p_other_cell(new Cell(p_healthy_state, p_other_cell_model));
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

        // For coverage, we also test TysonNovakCellCycleModel methods for a mutant cell
        p_cell_model->ResetForDivision();

        TysonNovakCellCycleModel* p_cell_model2 = static_cast<TysonNovakCellCycleModel*> (p_cell_model->CreateCellCycleModel());
        p_cell_model2->SetCellProliferativeType(STEM);

        MAKE_PTR(ApcOneHitCellMutationState, p_mutation);

        CellPtr p_stem_cell_2(new Cell(p_mutation, p_cell_model2));

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
    }

    /**
     * Test for Tyson & Novak self-cycling cells without having their
     * initial conditions reset. When using CVODE, the cell-cycle model
     * resets itself by halving the mass of the cell.
     * When not using CVODE, the cell-cycle model resets its initial
     * conditions, since the oscillatory solution computed using the Chaste
     * ODE solver is not stable.
     */
    void TestTysonNovakCellCycleModelSolver() throw(Exception)
    {
        // Set up simulation time
        unsigned num_timesteps = 100000;
        double standard_divide_time = 75.19/60.0;
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.1*standard_divide_time, num_timesteps);

        // Create cell-cycle model and associated cell
        TysonNovakCellCycleModel* p_repeating_cell_model = new TysonNovakCellCycleModel;
        p_repeating_cell_model->SetCellProliferativeType(STEM);

        MAKE_PTR(ApcOneHitCellMutationState, p_mutation);

        CellPtr p_tyson_novak_cell(new Cell(p_mutation, p_repeating_cell_model));
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

    void TestAlarcon2004OxygenBasedCellCycleModel() throw(Exception)
    {
        // Set up SimulationTime
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(30.0, 3);

        // Set up oxygen_concentration
        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);

        // For coverage, we create 1D, 2D and 3D instances
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);
        CellwiseData<2>::Instance()->SetConstantDataForTesting(oxygen_concentration);
        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        // Create cell-cycle models
        Alarcon2004OxygenBasedCellCycleModel* p_model_1d = new Alarcon2004OxygenBasedCellCycleModel();
        p_model_1d->SetDimension(1);
        p_model_1d->SetCellProliferativeType(STEM);

        Alarcon2004OxygenBasedCellCycleModel* p_model_2d = new Alarcon2004OxygenBasedCellCycleModel();
        p_model_2d->SetDimension(2);
        p_model_2d->SetCellProliferativeType(STEM);

        Alarcon2004OxygenBasedCellCycleModel* p_model_3d = new Alarcon2004OxygenBasedCellCycleModel();
        p_model_3d->SetDimension(3);
        p_model_3d->SetCellProliferativeType(STEM);

        // Create cells
        MAKE_PTR(WildTypeCellMutationState, p_state);

        CellPtr p_cell_1d(new Cell(p_state, p_model_1d));
        p_cell_1d->InitialiseCellCycleModel();

        CellPtr p_cell_2d(new Cell(p_state, p_model_2d));
        p_cell_2d->InitialiseCellCycleModel();

        CellPtr p_cell_3d(new Cell(p_state, p_model_3d));
        p_cell_3d->InitialiseCellCycleModel();

        // For coverage, we create another cell-cycle model that is identical to p_model_2d except for the ODE solver
        boost::shared_ptr<CellCycleModelOdeSolver<Alarcon2004OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver> >
            p_solver(CellCycleModelOdeSolver<Alarcon2004OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver>::Instance());
        p_solver->Initialise();

        Alarcon2004OxygenBasedCellCycleModel* p_other_model_2d = new Alarcon2004OxygenBasedCellCycleModel(p_solver);
        p_other_model_2d->SetDimension(2);
        p_other_model_2d->SetCellProliferativeType(STEM);

        CellPtr p_other_cell_2d(new Cell(p_state, p_other_model_2d));
        p_other_cell_2d->InitialiseCellCycleModel();

        // Check oxygen concentration is correct in cell-cycle model
        TS_ASSERT_DELTA(p_model_2d->GetProteinConcentrations()[5], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);

        TS_ASSERT_DELTA(p_other_model_2d->GetProteinConcentrations()[5], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_other_model_2d->ReadyToDivide(), false);

        // Divide the cells
        Alarcon2004OxygenBasedCellCycleModel* p_model_1d_2 = static_cast<Alarcon2004OxygenBasedCellCycleModel*> (p_model_1d->CreateCellCycleModel());
        p_model_1d_2->SetCellProliferativeType(STEM);
        CellPtr p_cell_1d_2(new Cell(p_state, p_model_1d_2));

        Alarcon2004OxygenBasedCellCycleModel* p_model_2d_2 = static_cast<Alarcon2004OxygenBasedCellCycleModel*> (p_model_2d->CreateCellCycleModel());
        p_model_2d_2->SetCellProliferativeType(STEM);
        CellPtr p_cell_2d_2(new Cell(p_state, p_model_2d_2));

        Alarcon2004OxygenBasedCellCycleModel* p_model_3d_2 = static_cast<Alarcon2004OxygenBasedCellCycleModel*> (p_model_3d->CreateCellCycleModel());
        p_model_3d_2->SetCellProliferativeType(STEM);
        CellPtr p_cell_3d_2(new Cell(p_state, p_model_3d_2));

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

        // Tidy up
        CellwiseData<2>::Destroy();

        // For coverage, create a 1D model
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);
        Alarcon2004OxygenBasedCellCycleModel* p_cell_model3 = new Alarcon2004OxygenBasedCellCycleModel();
        p_cell_model3->SetDimension(1);
        p_cell_model3->SetCellProliferativeType(STEM);

        CellPtr p_cell3(new Cell(p_state, p_cell_model3));
        p_cell3->InitialiseCellCycleModel();

        TS_ASSERT_DELTA(p_cell_model3->GetProteinConcentrations()[5], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);

        // Tidy up
        CellwiseData<1>::Destroy();
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
            p_model->SetCellProliferativeType(STEM);
            p_model->SetBirthTime(-1.5);
            static_cast<TysonNovakCellCycleModel*>(p_model)->SetDt(0.085);

            // We must create a cell to be able to initialise the cell cycle model's ODE system
            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            CellPtr p_cell(new Cell(p_healthy_state, p_model));
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
            TS_ASSERT_EQUALS(p_model2->GetCellProliferativeType(), STEM);
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

        std::vector<double> oxygen_concentration;
        oxygen_concentration.push_back(1.0);
        CellwiseData<1>::Instance()->SetConstantDataForTesting(oxygen_concentration);

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new Alarcon2004OxygenBasedCellCycleModel;

            p_model->SetDimension(1);
            p_model->SetCellProliferativeType(STEM);
            p_model->SetBirthTime(-1.5);
            static_cast<Alarcon2004OxygenBasedCellCycleModel*>(p_model)->SetDt(0.085);

            // We must create a cell to be able to initialise the cell cycle model's ODE system
            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            CellPtr p_cell(new Cell(p_healthy_state, p_model));
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
            TS_ASSERT_EQUALS(p_model2->GetCellProliferativeType(), STEM);
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
            CellwiseData<1>::Destroy();
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

        std::string alarcon_oxygen_based_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + alarcon_oxygen_based_results_dir + "alarcon_oxygen_based_results.parameters cell_based/test/data/TestCellCycleModels/alarcon_oxygen_based_results.parameters").c_str()), 0);

        // Test with TysonNovakCellCycleModel
        TysonNovakCellCycleModel tyson_novak_based_cell_cycle_model;
        TS_ASSERT_EQUALS(tyson_novak_based_cell_cycle_model.GetIdentifier(), "TysonNovakCellCycleModel");

        out_stream tyson_novak_based_parameter_file = output_file_handler.OpenOutputFile("tyson_novak_based_results.parameters");
        tyson_novak_based_cell_cycle_model.OutputCellCycleModelParameters(tyson_novak_based_parameter_file);
        tyson_novak_based_parameter_file->close();

        std::string tyson_novak_based_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + tyson_novak_based_results_dir + "tyson_novak_based_results.parameters cell_based/test/data/TestCellCycleModels/tyson_novak_based_results.parameters").c_str()), 0);
    }
};

#endif /*TESTODEBASEDCELLCYCLEMODELS_HPP_*/
