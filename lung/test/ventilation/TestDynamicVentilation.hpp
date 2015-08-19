/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTDYNAMICVENTILATION_HPP_
#define TESTDYNAMICVENTILATION_HPP_

#include <cxxtest/TestSuite.h>
#include "AirwayPropertiesCalculator.hpp"
#include "AirwayRemesher.hpp"
#include "MatrixVentilationProblem.hpp"
#include "MutableMesh.hpp"
#include "OutputFileHandler.hpp"
#include "ProgressReporter.hpp"
#include "TetrahedralMesh.hpp"
#include "FileFinder.hpp"
#include "SimpleBalloonAcinarUnit.hpp"

#include "boost/numeric/ublas/io.hpp"

#include <set>
#include <map>

#include "CommandLineArguments.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestDynamicVentilation : public CxxTest::TestSuite
{
public:


    void TestColemanDynamicVentilationSingleAirway() throw(Exception)
    {
        FileFinder mesh_finder("lung/test/data/single_branch", RelativeTo::ChasteSourceRoot);
        MatrixVentilationProblem problem(mesh_finder.GetAbsolutePath(), 0u);
        TetrahedralMesh<1,3>& r_mesh = problem.rGetMesh();

        //Initial conditions
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(0.0);
        problem.SetMeshInMilliMetres();

        double compliance = 0.1/98.0665/1e3;  //in m^3 / pa. Converted from 0.1 L/cmH2O per lung.

        double viscosity = 1.92e-5;               //Pa s
        double terminal_airway_radius = 0.0005;   //m
        double terminal_airway_length  = 0.005;   //m
        double terminal_airway_resistance = 8*viscosity*terminal_airway_length/(M_PI*SmallPow(terminal_airway_radius, 4));

        std::vector<double> pressures(r_mesh.GetNumNodes(), -1);
        std::vector<double> fluxes(r_mesh.GetNumNodes() - 1, -1);

        double ode_volume;

        VtkMeshWriter<1, 3> vtk_writer("TestDynamicVentilation", "single_branch", false);

        //For create an acinar balloon for the terminal node
        std::map<unsigned, SimpleBalloonAcinarUnit*> acinar_map;
        for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = r_mesh.GetBoundaryNodeIteratorBegin();
                             iter != r_mesh.GetBoundaryNodeIteratorEnd();
                             ++iter )
        {
            if ((*iter)->GetIndex() != 0u)
            {
                acinar_map[(*iter)->GetIndex()] = new SimpleBalloonAcinarUnit;
                acinar_map[(*iter)->GetIndex()]->SetCompliance(compliance);
            }
        }

        //Setup a simulation iterating between the flow solver and the acinar balloon.
        TimeStepper time_stepper(0.0, 1.0, 0.01);

        double pleural_pressure = 0.0;

        while (!time_stepper.IsTimeAtEnd())
        {
            pleural_pressure =  -2400*(sin((M_PI)*(time_stepper.GetNextTime())));

            //Solve corresponding backward Euler problem for testing
            double dt = time_stepper.GetNextTimeStep();
            ode_volume = (ode_volume - dt*pleural_pressure/terminal_airway_resistance)/(1 + dt/(terminal_airway_resistance*compliance));


            //Solve coupled problem
            for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = r_mesh.GetBoundaryNodeIteratorBegin();
                     iter != r_mesh.GetBoundaryNodeIteratorEnd();
                     ++iter )
            {
                if ((*iter)->GetIndex() != 0u)
                {
                    acinar_map[(*iter)->GetIndex()]->SetPleuralPressure(pleural_pressure);
                    acinar_map[(*iter)->GetIndex()]->ComputeExceptFlow(time_stepper.GetTime(), time_stepper.GetNextTime());

                    problem.SetPressureAtBoundaryNode(*(*iter), acinar_map[(*iter)->GetIndex()]->GetAirwayPressure());
                }
            }

            problem.Solve();
            problem.GetSolutionAsFluxesAndPressures(fluxes, pressures);

            for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = r_mesh.GetBoundaryNodeIteratorBegin();
                                iter != r_mesh.GetBoundaryNodeIteratorEnd();
                                ++iter )
            {
                if ((*iter)->GetIndex() != 0u)
                {
                    unsigned boundary_element_index = (*(*iter)->rGetContainingElementIndices().begin());

                    acinar_map[(*iter)->GetIndex()]->SetFlow(fluxes[boundary_element_index]);

                    double resistance = 0.0;
                    if(fluxes[(*iter)->GetIndex()] != 0.0)
                    {
                        resistance = std::fabs(pressures[(*iter)->GetIndex()]/fluxes[boundary_element_index]);
                    }
                    acinar_map[(*iter)->GetIndex()]->SetTerminalBronchioleResistance(resistance);
                    acinar_map[(*iter)->GetIndex()]->UpdateFlow(time_stepper.GetTime(), time_stepper.GetNextTime());
                }
            }

            TS_ASSERT_DELTA(ode_volume, acinar_map[5]->GetVolume(), 1e-6);

            time_stepper.AdvanceOneTimeStep();
        }
    }
};


#endif /* TESTDYNAMICVENTILATION_HPP_ */
