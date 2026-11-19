/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef TESTCARDIACCELLCOLLECTION_HPP_
#define TESTCARDIACCELLCOLLECTION_HPP_

#include <cxxtest/TestSuite.h>

#include <cmath>
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "CardiacCellCollection.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "LuoRudy1991.hpp"
#include "NobleVargheseKohlNoble1998a.hpp"
#include "SimpleStimulus.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCardiacCellCollection : public CxxTest::TestSuite
{
private:
    static constexpr double DT = 0.01;

    // Cells with differing stimuli, voltages and parameters. The stimuli do not change within any whole
    // millisecond, so solving a millisecond at a time with the stimulus held constant matches the CPU.
    template <typename CellType>
    std::vector<boost::shared_ptr<CellType> > MakeCells(unsigned numCells)
    {
        boost::shared_ptr<AbstractIvpOdeSolver> p_solver(new EulerIvpOdeSolver);
        std::vector<boost::shared_ptr<CellType> > cells;
        for (unsigned i = 0; i < numCells; i++)
        {
            const double magnitude = (i % 3 == 0) ? 0.0 : -25.5 - 0.1 * (i % 7);
            boost::shared_ptr<AbstractStimulusFunction> p_stimulus(new SimpleStimulus(magnitude, 1.995, 1.0 + i % 2));
            boost::shared_ptr<CellType> p_cell(new CellType(p_solver, p_stimulus));
            p_cell->SetTimestep(DT);
            p_cell->SetVoltage(p_cell->GetVoltage() + 0.05 * (i % 10));
            const std::string parameter_name = p_cell->GetSystemInformation()->rGetParameterNames()[0];
            p_cell->SetParameter(parameter_name, p_cell->GetParameter(parameter_name) * (1.0 + 0.01 * (i % 5)));
            cells.push_back(p_cell);
        }
        return cells;
    }

    // The tolerance is relative to each state variable, or to 1e-3 for smaller ones
    template <typename CellType>
    void CompareCells(const std::vector<boost::shared_ptr<CellType> >& rCells,
                      const std::vector<boost::shared_ptr<CellType> >& rExpectedCells,
                      double tolerance)
    {
        for (unsigned cell = 0; cell < rCells.size(); cell++)
        {
            const std::vector<double>& r_states = rCells[cell]->rGetStateVariables();
            const std::vector<double>& r_expected = rExpectedCells[cell]->rGetStateVariables();
            for (unsigned i = 0; i < r_states.size(); i++)
            {
                TS_ASSERT_DELTA(r_states[i], r_expected[i], tolerance * std::max(std::fabs(r_expected[i]), 1e-3));
            }
        }
    }

    template <typename CellType, typename ValueType>
    void CompareSolveWithCpu(double tolerance)
    {
        const unsigned num_cells = 300u; // More than one block of GPU threads
        std::vector<boost::shared_ptr<CellType> > cpu_cells = MakeCells<CellType>(num_cells);
        std::vector<boost::shared_ptr<CellType> > gpu_cells = MakeCells<CellType>(num_cells);

        CardiacCellCollection<CellType, ValueType> collection(gpu_cells);
        collection.SetTimestep(DT);
        TS_ASSERT_EQUALS(collection.GetNumberOfCells(), num_cells);

        for (unsigned ms = 0; ms < 6u; ms++)
        {
            for (unsigned cell = 0; cell < num_cells; cell++)
            {
                cpu_cells[cell]->SolveAndUpdateState(ms, ms + 1.0);
            }
            collection.SolveAndUpdateState(ms, ms + 1.0);
        }
        CompareCells(gpu_cells, cpu_cells, tolerance);

        // The stimulated cell 1 has depolarised compared to the unstimulated cell 0
        TS_ASSERT_LESS_THAN(gpu_cells[0]->GetVoltage() + 1.0, gpu_cells[1]->GetVoltage());
    }

    template <typename CellType, typename ValueType>
    void CompareComputeExceptVoltageWithCpu(double tolerance)
    {
        const unsigned num_cells = 10u;
        std::vector<boost::shared_ptr<CellType> > cpu_cells = MakeCells<CellType>(num_cells);
        std::vector<boost::shared_ptr<CellType> > gpu_cells = MakeCells<CellType>(num_cells);

        CardiacCellCollection<CellType, ValueType> collection(gpu_cells);
        collection.SetTimestep(DT);

        for (unsigned cell = 0; cell < num_cells; cell++)
        {
            // A voltage away from rest, so that the gating variables move
            cpu_cells[cell]->SetVoltage(-20.0 + cell / 3.0);
            gpu_cells[cell]->SetVoltage(-20.0 + cell / 3.0);
            cpu_cells[cell]->ComputeExceptVoltage(0.0, 1.0);
        }
        collection.ComputeExceptVoltage(0.0, 1.0);
        CompareCells(gpu_cells, cpu_cells, tolerance);

        for (unsigned cell = 0; cell < num_cells; cell++)
        {
            // Exactly unchanged, also in single precision
            TS_ASSERT_EQUALS(gpu_cells[cell]->GetVoltage(), -20.0 + cell / 3.0);
        }
    }

public:
    void TestSolveInDoublePrecision()
    {
        CompareSolveWithCpu<CellLuoRudy1991FromCellML, double>(1e-9);
        CompareSolveWithCpu<CellNobleVargheseKohlNoble1998aFromCellML, double>(1e-9);
    }

    void TestSolveInSinglePrecision()
    {
        CompareSolveWithCpu<CellLuoRudy1991FromCellML, float>(1e-3);
        CompareSolveWithCpu<CellNobleVargheseKohlNoble1998aFromCellML, float>(1e-3);
    }

    void TestComputeExceptVoltage()
    {
        CompareComputeExceptVoltageWithCpu<CellLuoRudy1991FromCellML, double>(1e-9);
        CompareComputeExceptVoltageWithCpu<CellLuoRudy1991FromCellML, float>(1e-3);
    }

    void TestTimestepThatDoesNotDivideTheInterval()
    {
        std::vector<boost::shared_ptr<CellLuoRudy1991FromCellML> > cpu_cells = MakeCells<CellLuoRudy1991FromCellML>(4u);
        std::vector<boost::shared_ptr<CellLuoRudy1991FromCellML> > gpu_cells = MakeCells<CellLuoRudy1991FromCellML>(4u);
        CardiacCellCollection<CellLuoRudy1991FromCellML, double> collection(gpu_cells);

        // Steps of 0.01, 0.01 and 0.005, as on the CPU
        collection.SetTimestep(DT);
        collection.SolveAndUpdateState(0.0, 0.025);
        for (unsigned cell = 0; cell < 4u; cell++)
        {
            cpu_cells[cell]->SolveAndUpdateState(0.0, 0.025);
        }
        CompareCells(gpu_cells, cpu_cells, 1e-12);
    }

    void TestExceptions()
    {
        std::vector<boost::shared_ptr<CellLuoRudy1991FromCellML> > no_cells;
        typedef CardiacCellCollection<CellLuoRudy1991FromCellML, double> Collection;
        TS_ASSERT_THROWS_THIS(Collection collection(no_cells), "A CardiacCellCollection needs at least one cell");

        std::vector<boost::shared_ptr<CellLuoRudy1991FromCellML> > cells = MakeCells<CellLuoRudy1991FromCellML>(1u);
        Collection collection(cells);
        TS_ASSERT_THROWS_CONTAINS(collection.SolveAndUpdateState(1.0, 0.5),
                                  "The simulation duration must be positive");
    }
};

#endif // TESTCARDIACCELLCOLLECTION_HPP_
