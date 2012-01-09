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


#ifndef _TESTPYCML_HPP_
#define _TESTPYCML_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <string>
#include <ctime>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "RunAndCheckIonicModels.hpp"

#include "Exception.hpp"
#include "SimpleStimulus.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "ArchiveLocationInfo.hpp"
#include "VectorHelperFunctions.hpp"

#include "LuoRudy1991.hpp"
#include "LuoRudy1991Opt.hpp"
#include "LuoRudy1991BackwardEuler.hpp"
#include "LuoRudy1991.hpp"

#include "TenTusscher2006Epi.hpp"
#include "TenTusscher2006EpiOpt.hpp"
#include "TenTusscher2006EpiBackwardEuler.hpp"

// Note: only using the optimised model, to test linking with chaste_libs=0!
#include "NobleVargheseKohlNoble1998aOpt.hpp"

#ifdef CHASTE_CVODE
#include "LuoRudy1991Cvode.hpp"
#include "LuoRudy1991CvodeOpt.hpp"
#endif // CHASTE_CVODE

class TestPyCml : public CxxTest::TestSuite
{
    template<typename VECTOR_TYPE>
    void CheckDerivedQuantities(AbstractParameterisedSystem<VECTOR_TYPE>& rCell,
                                const VECTOR_TYPE& rStateVec)
    {
        TS_ASSERT_EQUALS(rCell.GetNumberOfDerivedQuantities(), 2u);
        TS_ASSERT_EQUALS(rCell.GetDerivedQuantityIndex("FonRT"), 0u);
        TS_ASSERT_EQUALS(rCell.GetDerivedQuantityIndex("potassium_currents"), 1u);
        TS_ASSERT_EQUALS(rCell.GetDerivedQuantityUnits(0u), "per_millivolt");
        TS_ASSERT_EQUALS(rCell.GetDerivedQuantityUnits(1u), "microA_per_cm2");
        VECTOR_TYPE derived = rCell.ComputeDerivedQuantitiesFromCurrentState(0.0);
        const double FonRT = 0.037435728309031795;
        const double i_K_total = 1.0007;
        TS_ASSERT_EQUALS(GetVectorSize(derived), 2u);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 0), FonRT, 1e-12);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 1), i_K_total, 1e-4);
        DeleteVector(derived);
        derived = rCell.ComputeDerivedQuantities(0.0, rStateVec);
        TS_ASSERT_EQUALS(GetVectorSize(derived), 2u);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 0), FonRT, 1e-12);
        TS_ASSERT_DELTA(GetVectorComponent(derived, 1), i_K_total, 1e-4);
        DeleteVector(derived);
    }

    template<typename VECTOR_TYPE>
    void CheckParameter(AbstractParameterisedSystem<VECTOR_TYPE>& rCell)
    {
        TS_ASSERT_EQUALS(rCell.GetNumberOfParameters(), 2u);
        TS_ASSERT_EQUALS(rCell.GetParameterIndex("membrane_fast_sodium_current_conductance"), 0u);
        TS_ASSERT_EQUALS(rCell.GetParameterUnits(0u), "milliS_per_cm2");
        TS_ASSERT_EQUALS(rCell.GetParameter(0u), 23.0);
        rCell.SetParameter(0u, 0.1);
        TS_ASSERT_EQUALS(rCell.GetParameter(0u), 0.1);
        rCell.SetParameter(0u, 23.0);

        // and the system name...
        TS_ASSERT_EQUALS(rCell.GetSystemName(), "luo_rudy_1991");
    }

    template<typename VECTOR_TYPE>
    void CheckAttributes(AbstractParameterisedSystem<VECTOR_TYPE>& rCell)
    {
        TS_ASSERT_EQUALS(rCell.GetNumberOfAttributes(), 2u);
        TS_ASSERT(rCell.HasAttribute("SuggestedCycleLength"));
        TS_ASSERT_DELTA(rCell.GetAttribute("SuggestedCycleLength"), 750, 1e-12);
        TS_ASSERT(rCell.HasAttribute("SuggestedForwardEulerTimestep"));
        TS_ASSERT_DELTA(rCell.GetAttribute("SuggestedForwardEulerTimestep"), 0.005, 1e-12);
    }

    void CheckCai(AbstractCardiacCell& rCell, bool hasCai, double value=0.0)
    {
        if (hasCai)
        {
            TS_ASSERT_DELTA(rCell.GetIntracellularCalciumConcentration(), value, 1e-6);
        }
        else
        {
            TS_ASSERT_THROWS_THIS(rCell.GetIntracellularCalciumConcentration(),
                                  "AbstractCardiacCell::GetIntracellularCalciumConcentration() called. Either model has no [Ca_i] or method has not been implemented yet");
        }
    }

public:
    /**
     * This test is designed to quickly check that PyCml-generated code matches the Chaste interfaces,
     * and gives expected results.
     */
    void TestPyCmlCodeGeneration() throw(Exception)
    {
        clock_t ck_start, ck_end;

        //
        // Set up cells
        //
        double magnitude = -25.5;
        double duration  = 2.0 ;  // ms
        double when = 50.0; // ms
        boost::shared_ptr<SimpleStimulus> p_stimulus(new SimpleStimulus(magnitude, duration, when));
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        double end_time = 1000.0; //One second in milliseconds
        double i_ionic_end_time = 60.0; // ms
        double i_ionic = 1.9411; // test value

        // Normal model
        CellLuoRudy1991FromCellML normal(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(normal.GetVoltageIndex(), 0u);
        CheckCai(normal, true, 0.0002);
        double normal_initial_i_ionic = normal.GetIIonic();
        // Coverage
        normal.SetTimestep(HeartConfig::Instance()->GetOdeTimeStep());

        // Optimised model
        AbstractLookupTableCollection::EventHandler::Enable();
        CellLuoRudy1991FromCellMLOpt opt(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(opt.GetVoltageIndex(), 0u);
        CheckCai(opt, true, 0.0002);

        // Backward Euler optimised model
        CellLuoRudy1991FromCellMLBackwardEuler be(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(be.GetVoltageIndex(), 0u);
        CheckCai(be, true, 0.0002);

        // Check tables using AbstractLookupTableCollection interface
        TS_ASSERT(!normal.GetLookupTableCollection());
        AbstractLookupTableCollection* p_tables = opt.GetLookupTableCollection();
        TS_ASSERT(p_tables);
        TS_ASSERT_EQUALS(p_tables->GetKeyingVariableNames().size(), 2u);
        TS_ASSERT_EQUALS(p_tables->GetKeyingVariableNames()[0], "membrane_voltage");
        TS_ASSERT_EQUALS(p_tables->GetKeyingVariableNames()[1], "cytosolic_calcium_concentration");
        TS_ASSERT_EQUALS(p_tables->GetNumberOfTables("membrane_voltage"), 19u);
        TS_ASSERT_EQUALS(p_tables->GetNumberOfTables("cytosolic_calcium_concentration"), 1u);
        TS_ASSERT_THROWS_THIS(p_tables->GetNumberOfTables("non-var"), "Lookup table keying variable 'non-var' does not exist.");
        double min, max, step;
        p_tables->GetTableProperties("membrane_voltage", min, step, max);
        TS_ASSERT_DELTA(min, -150.0001, 1e-12);
        TS_ASSERT_DELTA(step, 0.01, 1e-12);
        TS_ASSERT_DELTA(max, 199.9999, 1e-12);
        p_tables->GetTableProperties("cytosolic_calcium_concentration", min, step, max);
        TS_ASSERT_DELTA(min, 0.00001, 1e-12);
        TS_ASSERT_DELTA(step, 0.0001, 1e-12);
        TS_ASSERT_DELTA(max, 30.00001, 1e-12);

        // Check set methods for coverage
        AbstractLookupTableCollection::EventHandler::Headings();
        AbstractLookupTableCollection::EventHandler::Report();
        p_tables->SetTimestep(0.1);
        p_tables->SetTableProperties("membrane_voltage", -100.0001, 0.01, 60.9999);
        p_tables->RegenerateTables();
        AbstractLookupTableCollection::EventHandler::Report();
        TS_ASSERT_THROWS_THIS(p_tables->SetTableProperties("membrane_voltage", -1, 0.03, 1),
                              "Table step size does not divide range between table limits.");
        p_tables->SetTimestep(HeartConfig::Instance()->GetOdeTimeStep());
        p_tables->SetTableProperties("membrane_voltage", -150.0001, 0.01, 199.9999);
        p_tables->RegenerateTables();
        AbstractLookupTableCollection::EventHandler::Report();

        // Check that the tables really exist!
        double v = opt.GetVoltage();
        opt.SetVoltage(-100000);
        TS_ASSERT_THROWS_CONTAINS(opt.GetIIonic(), "membrane_voltage outside lookup table range");
        opt.SetVoltage(v);

        be.SetVoltage(-100000);
        TS_ASSERT_THROWS_CONTAINS(be.GetIIonic(), "membrane_voltage outside lookup table range");
        be.SetVoltage(v);

        unsigned cai_index = opt.GetStateVariableIndex("cytosolic_calcium_concentration");
        double cai = opt.GetStateVariable(cai_index);
        opt.SetStateVariable(cai_index, -1.0);
        TS_ASSERT_THROWS_CONTAINS(opt.GetIIonic(), "cytosolic_calcium_concentration outside lookup table range");
        opt.SetStateVariable(cai_index, cai);

        be.SetStateVariable(cai_index, -1.0);
        TS_ASSERT_THROWS_CONTAINS(be.GetIIonic(), "cytosolic_calcium_concentration outside lookup table range");
        be.SetStateVariable(cai_index, cai);

        // Single parameter
        CheckParameter(normal);
        CheckParameter(opt);
        CheckParameter(be);

        // Derived variables
        CheckDerivedQuantities(normal, normal.GetInitialConditions());
        CheckDerivedQuantities(opt, opt.GetInitialConditions());
        CheckDerivedQuantities(be, be.GetInitialConditions());

        // Attributes
        CheckAttributes(normal);
        CheckAttributes(opt);
        CheckAttributes(be);

#ifdef CHASTE_CVODE
        // CVODE version
        CellLuoRudy1991FromCellMLCvode cvode_cell(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(cvode_cell.GetVoltageIndex(), 0u);
        // Optimised CVODE version
        CellLuoRudy1991FromCellMLCvodeOpt cvode_opt(p_solver, p_stimulus);
        TS_ASSERT_EQUALS(cvode_opt.GetVoltageIndex(), 0u);

        // Check tables using AbstractLookupTableCollection interface
        TS_ASSERT(!cvode_cell.GetLookupTableCollection());
        p_tables = cvode_opt.GetLookupTableCollection();
        TS_ASSERT(p_tables);
        TS_ASSERT_EQUALS(p_tables->GetKeyingVariableNames().size(), 2u);
        TS_ASSERT_EQUALS(p_tables->GetKeyingVariableNames()[0], "membrane_voltage");
        TS_ASSERT_EQUALS(p_tables->GetKeyingVariableNames()[1], "cytosolic_calcium_concentration");
        TS_ASSERT_EQUALS(p_tables->GetNumberOfTables("membrane_voltage"), 19u);
        TS_ASSERT_EQUALS(p_tables->GetNumberOfTables("cytosolic_calcium_concentration"), 1u);
        TS_ASSERT_THROWS_THIS(p_tables->GetNumberOfTables("non-var"), "Lookup table keying variable 'non-var' does not exist.");
        p_tables->GetTableProperties("membrane_voltage", min, step, max);
        TS_ASSERT_DELTA(min, -150.0001, 1e-12);
        TS_ASSERT_DELTA(step, 0.01, 1e-12);
        TS_ASSERT_DELTA(max, 199.9999, 1e-12);
        p_tables->GetTableProperties("cytosolic_calcium_concentration", min, step, max);
        TS_ASSERT_DELTA(min, 0.00001, 1e-12);
        TS_ASSERT_DELTA(step, 0.0001, 1e-12);
        TS_ASSERT_DELTA(max, 30.00001, 1e-12);

        // Check that the tables really exist!
        cvode_opt.SetVoltage(-100000);
        TS_ASSERT_THROWS_CONTAINS(cvode_opt.GetIIonic(), "membrane_voltage outside lookup table range");
        cvode_opt.SetVoltage(v);

        cvode_opt.SetStateVariable(cai_index, -1.0);
        TS_ASSERT_THROWS_CONTAINS(cvode_opt.GetIIonic(), "cytosolic_calcium_concentration outside lookup table range");
        cvode_opt.SetStateVariable(cai_index, cai);

        // Single parameter
        CheckParameter(cvode_cell);
        CheckParameter(cvode_opt);

        // Derived variables
        N_Vector inits = cvode_cell.GetInitialConditions();
        CheckDerivedQuantities(cvode_cell, inits);
        CheckDerivedQuantities(cvode_opt, inits);
        DeleteVector(inits);

        // Attributes
        CheckAttributes(cvode_cell);
        CheckAttributes(cvode_opt);
#endif // CHASTE_CVODE

        // Test the archiving code too
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("lr91-pycml.arch");

        // Save all (non-CVODE) cells at initial state
        {
            AbstractCardiacCell* const p_normal_cell = &normal;
            AbstractCardiacCell* const p_opt_cell = &opt;
            AbstractCardiacCell* const p_be_cell = &be;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_normal_cell;
            output_arch << p_opt_cell;
            output_arch << p_be_cell;
        }

        //
        // Solve and write to file
        //

        // Normal
        ck_start = clock();
        RunOdeSolverWithIonicModel(&normal,
                                   end_time,
                                   "Lr91FromPyCml");
        ck_end = clock();
        double normal_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tNormal: " << normal_time << std::endl;

        CheckCellModelResults("Lr91FromPyCml", "Lr91DelayedStim");

        RunOdeSolverWithIonicModel(&normal,
                                   i_ionic_end_time,
                                   "Lr91GetIIonic", 1000, false);
        TS_ASSERT_DELTA( normal.GetIIonic(), i_ionic, 1e-3);

        // Variant form of GetIIonic
        {
            std::vector<double> inits = normal.GetInitialConditions();
            TS_ASSERT_DELTA(normal.GetIIonic(&inits), normal_initial_i_ionic, 1e-12);
        }

        // With zero g_Na
        normal.SetParameter(0u, 0.0);
        normal.ResetToInitialConditions();
        RunOdeSolverWithIonicModel(&normal,
                                   end_time,
                                   "Lr91FromPyCmlZeroGna");
        CheckCellModelResults("Lr91FromPyCmlZeroGna");

        // Optimised
        ck_start = clock();
        RunOdeSolverWithIonicModel(&opt,
                                   end_time,
                                   "Lr91FromPyCmlOpt");
        ck_end = clock();
        double opt_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tOptimised: " << opt_time << std::endl;

        CompareCellModelResults("Lr91FromPyCml", "Lr91FromPyCmlOpt", 2e-3, true);

        RunOdeSolverWithIonicModel(&opt,
                                   i_ionic_end_time,
                                   "Lr91GetIIonicOpt", 1000, false);
        TS_ASSERT_DELTA( opt.GetIIonic(), i_ionic, 1e-3);

        // No stimulus at end time
        TS_ASSERT_DELTA(opt.GetIntracellularAreaStimulus(i_ionic_end_time), 0.0, 1e-12);

        // Backward Euler
        ck_start = clock();
        RunOdeSolverWithIonicModel(&be,
                                   end_time,
                                   "Lr91FromPyCmlBackwardEuler");
        ck_end = clock();
        double be_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tBackward Euler: " << be_time << std::endl;

        CompareCellModelResults("Lr91FromPyCml", "Lr91FromPyCmlBackwardEuler", 2e-2, true);

        RunOdeSolverWithIonicModel(&be,
                                   i_ionic_end_time,
                                   "Lr91GetIIonicBackwardEuler", 1000, false);
        TS_ASSERT_DELTA( be.GetIIonic(), i_ionic, 1e-3);

        // With zero g_Na
        be.SetParameter(0u, 0.0);
        be.ResetToInitialConditions();
        RunOdeSolverWithIonicModel(&be,
                                   end_time,
                                   "Lr91BEFromPyCmlZeroGna");
        CheckCellModelResults("Lr91BEFromPyCmlZeroGna", "Lr91FromPyCmlZeroGna", 2e-2);

#ifdef CHASTE_CVODE
        // CVODE
        double max_dt = 1.0; //ms
        ck_start = clock();
        OdeSolution cvode_solution = cvode_cell.Solve(0.0, end_time, max_dt, max_dt);
        cvode_solution.WriteToFile("TestIonicModels","Lr91FromPyCmlCvode","ms",1,false);
        ck_end = clock();
        double cvode_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tCVODE: " << cvode_time << std::endl;
        CompareCellModelResults("Lr91FromPyCml", "Lr91FromPyCmlCvode", 1e-1, true);
        // Coverage
        cvode_cell.SetVoltageDerivativeToZero();
        cvode_cell.ResetToInitialConditions();
        cvode_cell.Solve(0.0, i_ionic_end_time, max_dt);
        TS_ASSERT_DELTA(cvode_cell.GetIIonic(), 0.0, 1e-1); // Cell should be at rest
        cvode_cell.SetVoltageDerivativeToZero(false);
        // Check GetIIonic
        cvode_cell.ResetToInitialConditions();
        cvode_cell.Solve(0.0, i_ionic_end_time, max_dt);
        TS_ASSERT_DELTA(cvode_cell.GetIIonic(), i_ionic, 1e-1);

        // CVODE Optimised
        ck_start = clock();
        cvode_solution = cvode_opt.Solve(0.0, end_time, max_dt, max_dt);
        cvode_solution.WriteToFile("TestIonicModels","Lr91FromPyCmlCvodeOpt","ms",1,false);
        ck_end = clock();
        double cvode_opt_time = (double)(ck_end - ck_start)/CLOCKS_PER_SEC;
        std::cout << "\n\tCVODE Optimised: " << cvode_opt_time << std::endl;
        CompareCellModelResults("Lr91FromPyCmlCvode", "Lr91FromPyCmlCvodeOpt", 1e-1, true);
        // Coverage
        cvode_opt.SetVoltageDerivativeToZero();
        cvode_opt.ResetToInitialConditions();
        cvode_opt.Solve(0.0, i_ionic_end_time, max_dt);
        TS_ASSERT_DELTA(cvode_opt.GetIIonic(), 0.0, 1e-1); // Cell should be at rest
        cvode_opt.SetVoltageDerivativeToZero(false);
        // Check GetIIonic
        cvode_opt.ResetToInitialConditions();
        cvode_opt.Solve(0.0, i_ionic_end_time, max_dt);
        TS_ASSERT_DELTA(cvode_opt.GetIIonic(), i_ionic, 1e-1);

        // No stimulus at end time
        TS_ASSERT_DELTA(cvode_opt.GetIntracellularAreaStimulus(i_ionic_end_time), 0.0, 1e-12);
#endif // CHASTE_CVODE

        // Load and check simulation results still match
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCardiacCell* p_normal_cell;
            AbstractCardiacCell* p_opt_cell;
            AbstractCardiacCell* p_be_cell;
            input_arch >> p_normal_cell;
            input_arch >> p_opt_cell;
            input_arch >> p_be_cell;

            TS_ASSERT_EQUALS( p_normal_cell->GetNumberOfStateVariables(), 8u );
            TS_ASSERT_EQUALS( p_opt_cell->GetNumberOfStateVariables(), 8u );
            TS_ASSERT_EQUALS( p_be_cell->GetNumberOfStateVariables(), 8u );

            RunOdeSolverWithIonicModel(p_normal_cell,
                                       end_time,
                                       "Lr91FromPyCmlAfterArchive");
            CheckCellModelResults("Lr91FromPyCmlAfterArchive", "Lr91DelayedStim");

            RunOdeSolverWithIonicModel(p_opt_cell,
                                       end_time,
                                       "Lr91FromPyCmlOptAfterArchive");
            CompareCellModelResults("Lr91FromPyCml", "Lr91FromPyCmlOptAfterArchive", 2e-3, true);

            RunOdeSolverWithIonicModel(p_be_cell,
                                       end_time,
                                       "Lr91FromPyCmlBackwardEulerAfterArchive");
            CompareCellModelResults("Lr91FromPyCml", "Lr91FromPyCmlBackwardEulerAfterArchive", 1e-2, true);

            delete p_normal_cell;
            delete p_opt_cell;
            delete p_be_cell;
        }
    }

    void TestModelWithNoIntracellularCalcium() throw(Exception)
    {
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        // Normal model
        CellTenTusscher2006EpiFromCellML normal(p_solver, p_stimulus);
        normal.UseCellMLDefaultStimulus();
        CheckCai(normal, false);

        // Optimised model
        CellTenTusscher2006EpiFromCellMLOpt opt(p_solver, p_stimulus);
        opt.UseCellMLDefaultStimulus();
        CheckCai(opt, false);

        // Backward Euler model
        CellTenTusscher2006EpiFromCellMLBackwardEuler be(p_solver, p_stimulus);
        be.UseCellMLDefaultStimulus();
        CheckCai(be, false);

        // N98
        CellNobleVargheseKohlNoble1998aFromCellMLOpt n98opt(p_solver, p_stimulus);
        CheckCai(n98opt, false);
    }
};


#endif //_TESTPYCML_HPP_
