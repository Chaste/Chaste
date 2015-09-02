/*

Copyright (c) 2005-2015, University of Oxford.
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
#include "DynamicVentilationProblem.hpp"
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

#include "DynamicVentilationProblem.hpp"
#include "CommandLineArguments.hpp"
#include "PetscSetupAndFinalize.hpp"



class SimpleAcinarUnitFactory : public AbstractAcinarUnitFactory
{
public:
    virtual AbstractAcinarUnit* CreateAcinarUnitForNode(Node<3>* pNode)
    {
        SimpleBalloonAcinarUnit* p_acinus = new SimpleBalloonAcinarUnit;

        p_acinus->SetCompliance(0.1/98.0665/1e3);

        return p_acinus;
    }

    virtual double GetPleuralPressureForNode(double time, Node<3>* pNode)
    {
        return -2400*sin((M_PI)*(time));
    }
};


class TestDynamicVentilation : public CxxTest::TestSuite
{
public:

    void TestColemanDynamicVentilationSingleAirway() throw(Exception)
    {
#ifdef LUNG_USE_UMFPACK ///\todo This should really be runnable without UMFPACK, remove this if matrix solver is improved.

        FileFinder mesh_finder("lung/test/data/single_branch", RelativeTo::ChasteSourceRoot);
        SimpleAcinarUnitFactory factory;

        double compliance = 0.1/98.0665/1e3;  //in m^3 / pa. Converted from 0.1 L/cmH2O per lung.

        double viscosity = 1.92e-5;               //Pa s
        double terminal_airway_radius = 0.0005;   //m
        double terminal_airway_length  = 0.005;   //m
        double terminal_airway_resistance = 8*viscosity*terminal_airway_length/(M_PI*SmallPow(terminal_airway_radius, 4));
        double ode_volume = 0.0;

        DynamicVentilationProblem problem(&factory, mesh_finder.GetAbsolutePath(), 0u);
        problem.SetTimeStep(0.01);

        TimeStepper time_stepper(0.0, 1.0, 0.01);

        while (!time_stepper.IsTimeAtEnd())
        {
            //Solve corresponding backward Euler problem for testing
            double pleural_pressure =  factory.GetPleuralPressureForNode(time_stepper.GetNextTime(), NULL);

            double dt = time_stepper.GetNextTimeStep();
            ode_volume = (ode_volume - dt*pleural_pressure/terminal_airway_resistance)/(1 + dt/(terminal_airway_resistance*compliance));

            //Solve using DynamicVentilationProblem
            problem.SetEndTime(time_stepper.GetNextTime());
            problem.Solve();

            std::map<unsigned, AbstractAcinarUnit*>& r_acinar_map = problem.rGetAcinarUnitMap();
            TS_ASSERT_DELTA(ode_volume, r_acinar_map[5]->GetVolume(), 1e-6);

            time_stepper.AdvanceOneTimeStep();
        }
#else
        std::cout << "Warning: this test requires UMFPACK to execute correctly. " << std::endl;
#endif
    }

    void TestColemanDynamicVentilationThreeBifurcations() throw(Exception)
    {
#ifdef LUNG_USE_UMFPACK ///\todo This should really be runnable without UMFPACK, remove this if matrix solver is improved.
        FileFinder mesh_finder("continuum_mechanics/test/data/three_bifurcations", RelativeTo::ChasteSourceRoot);
        MatrixVentilationProblem problem(mesh_finder.GetAbsolutePath(), 0u);
        TetrahedralMesh<1,3>& r_mesh = problem.rGetMesh();

        //Initial conditions
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(0.0);
        problem.SetMeshInMilliMetres();

        //The three bifurcation mesh defines a fully symmetric three bifurcation airway tree.
        //The composite ventilation problem is then equivalent to a trumpet problem connected
        //to an acinus with compliance equal to the total compliance of all the acini.
        double total_compliance = 0.1/98.0665/1e3;  //in m^3 / pa. Converted from 0.1 L/cmH2O per lung to four acinar compartments
        double acinar_compliance = total_compliance/4.0;

        double viscosity = 1.92e-5;               //Pa s
        double terminal_airway_radius = 0.00005;   //m
        double resistance_per_unit_length = 8*viscosity/(M_PI*SmallPow(terminal_airway_radius, 4));
        //All airways in the mesh have radius 0.05 mm. The first branch is 3mm long, the others are 5mm.
        double total_airway_resistance = (0.003 + 0.005/2 + 0.005/4)*resistance_per_unit_length;

        std::vector<double> pressures(r_mesh.GetNumNodes(), -1);
        std::vector<double> fluxes(r_mesh.GetNumNodes() - 1, -1);

        double ode_volume;

        //For create an acinar balloon for the terminal node
        std::map<unsigned, SimpleBalloonAcinarUnit*> acinar_map;
        for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = r_mesh.GetBoundaryNodeIteratorBegin();
                             iter != r_mesh.GetBoundaryNodeIteratorEnd();
                             ++iter )
        {
            if ((*iter)->GetIndex() != 0u)
            {
                acinar_map[(*iter)->GetIndex()] = new SimpleBalloonAcinarUnit;
                acinar_map[(*iter)->GetIndex()]->SetCompliance(acinar_compliance);
            }
        }

        //Setup a simulation iterating between the flow solver and the acinar balloon.
        TimeStepper time_stepper(0.0, 1.0, 0.001);

        double pleural_pressure = 0.0;

        while (!time_stepper.IsTimeAtEnd())
        {
            pleural_pressure =  -2400*(sin((M_PI)*(time_stepper.GetNextTime())));

            //Solve corresponding backward Euler problem for testing
            double dt = time_stepper.GetNextTimeStep();
            ode_volume = (ode_volume - dt*pleural_pressure/total_airway_resistance)/(1 + dt/(total_airway_resistance*acinar_compliance));


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

            double total_acinar_volume = acinar_map[4]->GetVolume() + acinar_map[5]->GetVolume() + acinar_map[6]->GetVolume() + acinar_map[7]->GetVolume();
            TS_ASSERT_DELTA(ode_volume, total_acinar_volume, 1e-12);

            time_stepper.AdvanceOneTimeStep();
        }
#else
        std::cout << "Warning: this test requires UMFPACK to execute correctly. " << std::endl;
#endif
        }


    void TestColemanDynamicVentilationOtisBifurcations() throw(Exception)
    {
#ifdef LUNG_USE_UMFPACK ///\todo This should really be runnable without UMFPACK, remove this if matrix solver is improved.
       FileFinder mesh_finder("lung/test/data/otis_bifurcation", RelativeTo::ChasteSourceRoot);
       MatrixVentilationProblem problem(mesh_finder.GetAbsolutePath(), 0u);
       TetrahedralMesh<1,3>& r_mesh = problem.rGetMesh();

       //Initial conditions
       problem.SetOutflowPressure(0.0);
       problem.SetConstantInflowPressures(0.0);
       problem.SetMeshInMilliMetres();
       problem.SetRadiusOnEdge();

       //The otis bifurcation mesh defines two branches of unequal radii leading to two acini
       //Analytical results for this system can be found in Otis et al. Journal of Applied Physiology 1956
       double total_compliance = 0.1/98.0665/1e3;  //in m^3 / pa. Converted from 0.1 L/cmH2O per lung to four acinar compartments

       double viscosity = 1.92e-5;               //Pa s
       double radius_zero = 0.002;     //m
       double radius_one =  0.0002;    //m
       double radius_two =  0.001;     //m
       double length = 0.001; //m
       double R0 = 8*length*viscosity/(M_PI*SmallPow(radius_zero, 4));
       double R1 = 8*length*viscosity/(M_PI*SmallPow(radius_one, 4));
       double R2 = 8*length*viscosity/(M_PI*SmallPow(radius_two, 4));

       double C1 = total_compliance/2.0;
       double C2 = total_compliance/2.0;
       double T1 = C1*R1;
       double T2 = C2*R2;

       double frequency = 2; //Hz
       double omega = 2*M_PI*frequency;

       double effective_compliance = (SmallPow(omega, 2)*SmallPow(T2*C1 + T1*C2, 2) + SmallPow(C1 + C2, 2)) /
                                       (SmallPow(omega, 2)*(SmallPow(T1,2)*C2 + SmallPow(T2,2)*C2) + C1 + C2);

       double effective_resistance = R0 +
                                     (SmallPow(omega,2)*T1*T2*(T2*C1 + T1*C2) + (T1*C1 + T2*C2)) /
                                      (SmallPow(omega,2)*SmallPow(T2*C1 + T1*C2,2) + SmallPow(C1 + C2, 2));

       double theta = std::atan(1/(omega*effective_resistance*effective_compliance));

       double delta_p = 500;

       std::vector<double> pressures(r_mesh.GetNumNodes(), -1);
       std::vector<double> fluxes(r_mesh.GetNumNodes() - 1, -1);

       double expected_tidal_volume = effective_compliance*delta_p*std::sin(theta);

       //For create an acinar balloon for the terminal node
       std::map<unsigned, SimpleBalloonAcinarUnit*> acinar_map;
       acinar_map[2] = new SimpleBalloonAcinarUnit;
       acinar_map[2]->SetCompliance(C1);
       acinar_map[3] = new SimpleBalloonAcinarUnit;
       acinar_map[3]->SetCompliance(C2);

       //Setup a simulation iterating between the flow solver and the acinar balloon.
       TimeStepper time_stepper(0.0, 20.0, 0.0005);

       double pleural_pressure = 0.0;

       double min_total_volume = 0.0;
       double max_total_volume = 0.0;

       while (!time_stepper.IsTimeAtEnd())
       {
           pleural_pressure =  -delta_p/2.0*(sin(omega*(time_stepper.GetNextTime())));

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

           //calculate max and minimum volume
           if(time_stepper.GetTime() > 16.0)
           {
               double total_volume = acinar_map[2]->GetVolume() + acinar_map[3]->GetVolume();

               if(min_total_volume > total_volume)
               {
                   min_total_volume = total_volume;
               }
               if(max_total_volume < total_volume)
               {
                   max_total_volume = total_volume;
               }
           }

           time_stepper.AdvanceOneTimeStep();
       }

       TS_ASSERT_DELTA(expected_tidal_volume, max_total_volume - min_total_volume, 1e-7);
#else
           std::cout << "Warning: this test requires UMFPACK to execute correctly. " << std::endl;
#endif
    }
};


#endif /* TESTDYNAMICVENTILATION_HPP_ */
