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
#include "SimpleBalloonExplicitAcinarUnit.hpp"

#include "boost/numeric/ublas/io.hpp"

#include <set>
#include <map>

#include "DynamicVentilationProblem.hpp"
#include "CommandLineArguments.hpp"
#include "PetscSetupAndFinalize.hpp"


template <typename ACINAR_UNIT = SimpleBalloonAcinarUnit> class SimpleAcinarUnitFactory : public AbstractAcinarUnitFactory
{
public:
    SimpleAcinarUnitFactory(double acinarCompliance,
                            double pleuralPressureAmplitude,
                            double frequency = 0.5) : mAcinarCompliance(acinarCompliance),
                                                      mPleuralPressureAmplitude(pleuralPressureAmplitude),
                                                      mFrequency(frequency)
    {}

    virtual AbstractAcinarUnit* CreateAcinarUnitForNode(Node<3>* pNode)
    {
        ACINAR_UNIT* p_acinus = new ACINAR_UNIT;

        p_acinus->SetCompliance(mAcinarCompliance);

        return p_acinus;
    }

    virtual double GetPleuralPressureForNode(double time, Node<3>* pNode)
    {
        return -mPleuralPressureAmplitude*sin(2*M_PI*mFrequency*time);
    }

private:
    double mAcinarCompliance;
    double mPleuralPressureAmplitude;
    double mFrequency;
};


class TestDynamicVentilation : public CxxTest::TestSuite
{
public:

    void TestColemanDynamicVentilationSingleAirway()
    {
#if defined(LUNG_USE_UMFPACK) || defined(LUNG_USE_KLU)
        FileFinder mesh_finder("lung/test/data/single_branch", RelativeTo::ChasteSourceRoot);

        double compliance = 0.1/98.0665/1e3;  //in m^3 / pa. Converted from 0.1 L/cmH2O per lung.

        SimpleAcinarUnitFactory<> factory(compliance, 2400.0);

        double viscosity = 1.92e-5;               //Pa s
        double terminal_airway_radius = 0.0005;   //m
        double terminal_airway_length  = 0.005;   //m
        double terminal_airway_resistance = 8*viscosity*terminal_airway_length/(M_PI*SmallPow(terminal_airway_radius, 4));
        double ode_volume = 0.0;

        DynamicVentilationProblem problem(&factory, mesh_finder.GetAbsolutePath(), 0u);
        problem.rGetMatrixVentilationProblem().SetMeshInMilliMetres();
        problem.rGetMatrixVentilationProblem().SetOutflowPressure(0.0);

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

        // Coverage
        TS_ASSERT_THROWS_NOTHING(factory.GetMesh());
#else
        std::cout << "Warning: This test needs a direct solver (UMFPACK or KLU) to execute correctly." << std::endl;
#endif
    }

    void TestColemanDynamicVentilationThreeBifurcations()
    {
        FileFinder mesh_finder("lung/test/data/three_bifurcations", RelativeTo::ChasteSourceRoot);

        //The three bifurcation mesh defines a fully symmetric three bifurcation airway tree.
        //The composite ventilation problem is then equivalent to a trumpet problem connected
        //to an acinus with compliance equal to the total compliance of all the acini.
        double total_compliance = 0.1/98.0665/1e3;  //in m^3 / pa. Converted from 0.1 L/cmH2O per lung to four acinar compartments
        double acinar_compliance = total_compliance/4.0;

        SimpleAcinarUnitFactory<> factory(acinar_compliance, 2400.0);
        TS_ASSERT_THROWS_CONTAINS(factory.GetMesh(), "The mesh object has not been set in the acinar unit factory");

        double viscosity = 1.92e-5;               //Pa s
        double terminal_airway_radius = 0.00005;   //m
        double resistance_per_unit_length = 8*viscosity/(M_PI*SmallPow(terminal_airway_radius, 4));
        //All airways in the mesh have radius 0.05 mm. The first branch is 3mm long, the others are 5mm.
        double total_airway_resistance = (0.003 + 0.005/2 + 0.005/4)*resistance_per_unit_length;

        double ode_volume = 0.0;

        //Setup a simulation iterating between the flow solver and the acinar balloon.
        DynamicVentilationProblem problem(&factory, mesh_finder.GetAbsolutePath(), 0u);
        problem.rGetMatrixVentilationProblem().SetOutflowPressure(0.0);
        problem.rGetMatrixVentilationProblem().SetMeshInMilliMetres();
        problem.SetTimeStep(0.01);
        factory.GetNumberOfAcini();

        TimeStepper time_stepper(0.0, 1.0, 0.01);

        while (!time_stepper.IsTimeAtEnd())
        {
            //Solve corresponding backward Euler problem for testing
            double pleural_pressure =  factory.GetPleuralPressureForNode(time_stepper.GetNextTime(), NULL);

            double dt = time_stepper.GetNextTimeStep();
            ode_volume = (ode_volume - dt*pleural_pressure/total_airway_resistance)/(1 + dt/(total_airway_resistance*total_compliance));

            //Solve using DynamicVentilationProblem
            problem.SetEndTime(time_stepper.GetNextTime());
            problem.Solve();

            std::map<unsigned, AbstractAcinarUnit*>& r_acinar_map = problem.rGetAcinarUnitMap();
            TS_ASSERT_DELTA(ode_volume, r_acinar_map[5]->GetVolume(), 1e-6);

            time_stepper.AdvanceOneTimeStep();
        }

        //Solve for longer and write output to VTK
        problem.SetSamplingTimeStepMultiple(10u);
        problem.SetOutputDirectory("TestDynamicVentilation");
        problem.SetOutputFilenamePrefix("three_bifurcations");
        problem.SetWriteVtkOutput();
        problem.SetEndTime(1.5);
        problem.Solve();

#ifdef CHASTE_VTK
        std::string filepath = OutputFileHandler::GetChasteTestOutputDirectory() + "TestDynamicVentilation/";
        std::string basename = filepath + "three_bifurcations";
        FileFinder vtu_file(basename + ".vtu", RelativeTo::Absolute);
        TS_ASSERT(vtu_file.Exists());
#endif
    }


    void TestColemanDynamicVentilationOtisBifurcations()
    {
#if defined(LUNG_USE_UMFPACK) || defined(LUNG_USE_KLU)
       FileFinder mesh_finder("lung/test/data/otis_bifurcation", RelativeTo::ChasteSourceRoot);

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

       SimpleAcinarUnitFactory<> factory(C1, delta_p/2.0, frequency);

       DynamicVentilationProblem problem(&factory, mesh_finder.GetAbsolutePath(), 0u);
       problem.rGetMatrixVentilationProblem().SetMeshInMilliMetres();
       problem.rGetMatrixVentilationProblem().SetRadiusOnEdge();
       problem.rGetMatrixVentilationProblem().SetOutflowPressure(0.0);

       double expected_tidal_volume = effective_compliance*delta_p*std::sin(theta);

       //Setup a simulation iterating between the flow solver and the acinar balloon.
       problem.SetTimeStep(0.0005);
       problem.SetEndTime(16.0);
       problem.Solve(); //Solve to 16s to allow the problem to equilibriate

       double min_total_volume = 0.0;
       double max_total_volume = 0.0;

       //Now solve the last four seconds recording the tidal volume.

       TimeStepper time_stepper(16.0, 20.0, 0.0005);

       while (!time_stepper.IsTimeAtEnd())
       {
           problem.SetEndTime(time_stepper.GetNextTime());
           problem.Solve();

           std::map<unsigned, AbstractAcinarUnit*>& r_acinar_map = problem.rGetAcinarUnitMap();
           double total_volume = r_acinar_map[2]->GetVolume() + r_acinar_map[3]->GetVolume();

           if (min_total_volume > total_volume)
           {
               min_total_volume = total_volume;
           }
           if (max_total_volume < total_volume)
           {
               max_total_volume = total_volume;
           }

           time_stepper.AdvanceOneTimeStep();
       }

       TS_ASSERT_DELTA(expected_tidal_volume, max_total_volume - min_total_volume, 1e-7);
#else
        std::cout << "Warning: This test needs a direct solver (UMFPACK or KLU) to execute correctly." << std::endl;
#endif
    }

    void TestColemanVsExplicitWithPedley()
    {
#if defined(LUNG_USE_UMFPACK) || defined(LUNG_USE_KLU)
        //This test compares an acinar unit using an explicit coupling scheme against an
        //acinar unit using the Coleman coupling scheme. Both use dynamic (Pedley) airway
        //resistance. Compliance & resistance are chosen to obtain similar dt bounds for
        //the explicit scheme as commonly found for full simulations

        FileFinder mesh_finder("lung/test/data/single_branch_no_intermediate", RelativeTo::ChasteSourceRoot);

        double compliance = 0.1/98.0665/1e3/30000;  //in m^3 / pa. Converted from 0.1 L/cmH2O per lung.

        //Theory says that dt < 2*R*C gives stability in the explicit scheme when using Poiseuille flow.
        //As R_poiseuille < R_pedley, this acts as a lower bound for stability here.
        //Nb, we use radius = 2 mm, which results in a very tight dt bound, to ensure that Pedley resistance is
        //used.
        double dt = 1e-6;
        double end_time = 0.15;

        double viscosity = 1.92e-5;               //Pa s
        double terminal_airway_radius = 0.002;   //m
        double terminal_airway_length  = 0.005;   //m
        double terminal_airway_resistance = 8*viscosity*terminal_airway_length/(M_PI*SmallPow(terminal_airway_radius, 4));
        TS_ASSERT_LESS_THAN(dt, 2*terminal_airway_resistance*compliance);

        double coleman_end_volume;
        double explicit_end_volume;

        {
            SimpleAcinarUnitFactory<SimpleBalloonExplicitAcinarUnit> factory(compliance, 2400.0);

            DynamicVentilationProblem problem(&factory, mesh_finder.GetAbsolutePath(), 0u);
            problem.rGetMatrixVentilationProblem().SetOutflowPressure(0.0);
            problem.rGetMatrixVentilationProblem().SetMeshInMilliMetres();
            problem.rGetMatrixVentilationProblem().SetDynamicResistance();
            problem.SetTimeStep(dt);
            problem.SetEndTime(end_time);

            problem.Solve();

            std::map<unsigned, AbstractAcinarUnit*>& r_acinar_map = problem.rGetAcinarUnitMap();
            explicit_end_volume = r_acinar_map[1]->GetVolume();
        }

        {
            SimpleAcinarUnitFactory<> factory(compliance, 2400.0);

            DynamicVentilationProblem problem(&factory, mesh_finder.GetAbsolutePath(), 0u);
            problem.rGetMatrixVentilationProblem().SetOutflowPressure(0.0);
            problem.rGetMatrixVentilationProblem().SetMeshInMilliMetres();
            problem.rGetMatrixVentilationProblem().SetDynamicResistance();
            problem.SetTimeStep(dt);
            problem.SetEndTime(end_time);

            problem.Solve();

            std::map<unsigned, AbstractAcinarUnit*>& r_acinar_map = problem.rGetAcinarUnitMap();
            coleman_end_volume = r_acinar_map[1]->GetVolume();
        }

        TS_ASSERT_DELTA(explicit_end_volume, coleman_end_volume, 1e-12);
#else
        std::cout << "Warning: This test needs a direct solver (UMFPACK or KLU) to execute correctly." << std::endl;
#endif
    }
};

#endif /* TESTDYNAMICVENTILATION_HPP_ */
