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


#ifndef _TESTPYCML_HPP_
#define _TESTPYCML_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <string>

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

#include "DiFrancescoNoble1985.hpp"
#include "DiFrancescoNoble1985Opt.hpp"
#include "DiFrancescoNoble1985BackwardEuler.hpp"

// Note: only using the optimised model, to test linking with chaste_libs=0!
#include "NobleVargheseKohlNoble1998aOpt.hpp"

#ifdef CHASTE_CVODE
#include "LuoRudy1991Cvode.hpp"
#include "LuoRudy1991CvodeOpt.hpp"
#endif // CHASTE_CVODE

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

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
        TS_ASSERT_EQUALS(rCell.GetNumberOfParameters(), 3u);
        TS_ASSERT_EQUALS(rCell.GetParameterIndex("membrane_fast_sodium_current_conductance"), 1u);
        TS_ASSERT_EQUALS(rCell.GetParameterUnits(1u), "milliS_per_cm2");
        TS_ASSERT_EQUALS(rCell.GetParameter(1u), 23.0);
        rCell.SetParameter(1u, 0.1);
        TS_ASSERT_EQUALS(rCell.GetParameter(1u), 0.1);
        TS_ASSERT_EQUALS(rCell.GetParameter("membrane_fast_sodium_current_conductance"), 0.1);
        rCell.SetParameter(1u, 23.0);
    }

    template<typename VECTOR_TYPE>
    void CheckAttributes(AbstractParameterisedSystem<VECTOR_TYPE>& rCell)
    {
        TS_ASSERT_EQUALS(rCell.GetNumberOfAttributes(), 2u);
        TS_ASSERT(rCell.HasAttribute("SuggestedCycleLength"));
        TS_ASSERT_DELTA(rCell.GetAttribute("SuggestedCycleLength"), 750, 1e-12);
        TS_ASSERT(rCell.HasAttribute("SuggestedForwardEulerTimestep"));
        TS_ASSERT_DELTA(rCell.GetAttribute("SuggestedForwardEulerTimestep"), 0.005, 1e-12);

        // The system name
        TS_ASSERT_EQUALS(rCell.GetSystemName(), "luo_rudy_1991");
        // Free variable
        TS_ASSERT_EQUALS(rCell.GetSystemInformation()->GetFreeVariableName(), "environment_time");
        TS_ASSERT_EQUALS(rCell.GetSystemInformation()->GetFreeVariableUnits(), "millisecond");
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
                                  "AbstractCardiacCellInterface::GetIntracellularCalciumConcentration() called. Either model has no [Ca_i] or method has not been implemented yet");
        }
    }

public:
    /**
     * This test is designed to quickly check that PyCml-generated code matches the Chaste interfaces,
     * and gives expected results.
     */
    void TestPyCmlCodeGeneration()
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
        TS_ASSERT_DELTA(step, 0.001, 1e-12);
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
        TS_ASSERT_EQUALS(p_tables->GetNumberOfTables("membrane_voltage"), 54u);
        TS_ASSERT_EQUALS(p_tables->GetNumberOfTables("cytosolic_calcium_concentration"), 2u);
        TS_ASSERT_THROWS_THIS(p_tables->GetNumberOfTables("non-var"), "Lookup table keying variable 'non-var' does not exist.");
        p_tables->GetTableProperties("membrane_voltage", min, step, max);
        TS_ASSERT_DELTA(min, -150.0001, 1e-12);
        TS_ASSERT_DELTA(step, 0.001, 1e-12);
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
        N_Vector vec_inits = cvode_cell.GetInitialConditions();
        CheckDerivedQuantities(cvode_cell, vec_inits);
        CheckDerivedQuantities(cvode_opt, vec_inits);
        DeleteVector(vec_inits);

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
        normal.SetParameter("membrane_fast_sodium_current_conductance", 0.0);
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
        be.SetParameter("membrane_fast_sodium_current_conductance", 0.0);
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
        // Upon a reset to initial conditions, when SetVoltageDerivativeToZero has been called we need to do the following line:
        cvode_cell.SetFixedVoltage(cvode_cell.GetVoltage());
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
        cvode_opt.ResetToInitialConditions();
        cvode_opt.SetVoltageDerivativeToZero();
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

    void TestModelWithNoIntracellularCalcium()
    {
        boost::shared_ptr<AbstractStimulusFunction> p_stimulus;
        boost::shared_ptr<EulerIvpOdeSolver> p_solver(new EulerIvpOdeSolver);

        // Normal model
        CellDiFrancescoNoble1985FromCellML normal(p_solver, p_stimulus);
        CheckCai(normal, false);

        // Optimised model
        CellDiFrancescoNoble1985FromCellMLOpt opt(p_solver, p_stimulus);
        CheckCai(opt, false);

        // Backward Euler model
        CellDiFrancescoNoble1985FromCellMLBackwardEuler be(p_solver, p_stimulus);
        CheckCai(be, false);

        // N98
        CellNobleVargheseKohlNoble1998aFromCellMLOpt n98opt(p_solver, p_stimulus);
        n98opt.UseCellMLDefaultStimulus();
        CheckCai(n98opt, false);
    }
};


#endif //_TESTPYCML_HPP_
