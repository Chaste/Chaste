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

#ifndef _TESTVENTILATIONPROBLEM_HPP_
#define _TESTVENTILATIONPROBLEM_HPP_

#include <cxxtest/TestSuite.h>
#include <queue>

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "VentilationProblem.hpp"
#include "Warnings.hpp"

void LinearTimeBCs(AbstractVentilationProblem* pProblem, TimeStepper& rTimeStepper, const Node<3>& rNode)
{
    double time = rTimeStepper.GetTime();
    pProblem->SetPressureAtBoundaryNode(rNode, 15*time);
}

void SineBCs(AbstractVentilationProblem* pProblem, TimeStepper& rTimeStepper, const Node<3>& rNode)
{
    double time = rTimeStepper.GetTime();
    pProblem->SetPressureAtBoundaryNode(rNode, 15.0*sin(time*5.0/(2.0*M_PI)));
}


void GravitationalBCs(AbstractVentilationProblem* pProblem, TimeStepper& rTimeStepper, const Node<3>& rNode)
{
    double x_max =  6.0;
    double x_min = -6.0;
    double delta_p = 0.5;

    double x = rNode.rGetLocation()[0];

    double pressure = delta_p * (x - x_min)/(x_max - x_min);

    pProblem->SetPressureAtBoundaryNode(rNode, pressure);
}


class TestVentilationProblem : public CxxTest::TestSuite
{
private:

public:
    void TestOnBranch()
    {
        VentilationProblem problem("mesh/test/data/y_branch_3d_mesh", 0u);
        problem.SetMeshInMilliMetres();
        TS_ASSERT_DELTA(problem.GetViscosity(), 1.92e-5, 1e-12);
        problem.SetViscosity(1.0);
        problem.SetOutflowPressure(100);
        problem.SetConstantInflowPressures(0.0);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(flux[0], 0.092813610, 1e-8);
        TS_ASSERT_DELTA(flux[1], 0.046406805, 1e-8);
        TS_ASSERT_DELTA(flux[2], 0.046406805, 1e-8);
        TS_ASSERT_DELTA(pressure[0], 100.0, 1e-6); //BC
        TS_ASSERT_DELTA(pressure[1], 84.5107, 1e-4);
        TS_ASSERT_DELTA(pressure[2], 0.0, 1e-6); //BC
        TS_ASSERT_DELTA(pressure[3], 0.0, 1e-6); //BC
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "small_conical");
#endif
    }
    void TestOnBranchCylindrical()
    {
        VentilationProblem problem("mesh/test/data/y_branch_3d_mesh");
        problem.SetMeshInMilliMetres();
        //These values are equivalent to Swan et al. 2012. 10.1016/j.jtbi.2012.01.042 (page 224)
        TS_ASSERT_DELTA(problem.GetViscosity(), 1.92e-5, 1e-12);
        problem.SetViscosity(1.0);
        TS_ASSERT_DELTA(problem.GetDensity(), 1.15, 1e-10);
        problem.SetDensity(1.0);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(100);
        problem.SetConstantInflowPressures(0.0);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(flux[0], 0.019932007, 1e-8);
        TS_ASSERT_DELTA(flux[1], 0.009966003, 1e-8);
        TS_ASSERT_DELTA(flux[2], 0.009966003, 1e-8);
        TS_ASSERT_DELTA(pressure[0], 100.0, 1e-6); //BC
        TS_ASSERT_DELTA(pressure[1], 91.8789, 1e-4);
        TS_ASSERT_DELTA(pressure[2], 0.0, 1e-6); //BC
        TS_ASSERT_DELTA(pressure[3], 0.0, 1e-6); //BC

#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "small_cylindrical");
#endif
    }

    void TestThreeBifurcationsWithRadiusOnEdgeFile()
    {
        VentilationProblem problem("lung/test/data/three_bifurcations", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(0.00148608);
        problem.Solve();
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "three_bifurcations");
#endif
        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 0.0006,   1e-4);
        TS_ASSERT_DELTA(pressure[2], 0.001274, 1e-4);
        TS_ASSERT_DELTA(pressure[3], 0.001274, 1e-4);
        TS_ASSERT_DELTA(pressure[4], 0.00148608, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[5], 0.00148608, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[6], 0.00148608, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[7], 0.00148608, 1e-8); //BC

        // Total flux
        // TS_ASSERT_DELTA(flux[0], -2.8143e-5, 1e-8); // Mesh in meters
        TS_ASSERT_DELTA(flux[0], -2.8143e-14, 1e-15); //Mesh in millimeters
        TS_ASSERT_DELTA(problem.GetFluxAtOutflow(), -2.8143e-14, 1e-15);
    }

    void TestThreeBifurcations()
    {
        VentilationProblem problem("lung/test/data/three_bifurcations", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(15.0);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 6.66666,   1e-4);
        TS_ASSERT_DELTA(pressure[2], 12.22223, 1e-4);
        TS_ASSERT_DELTA(pressure[3], 12.22222, 1e-4);
        TS_ASSERT_DELTA(pressure[4], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[5], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[6], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[7], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(flux[0], -2.8407e-10 , 1e-13); // (Outflow flux)
        TS_ASSERT_DELTA(flux[3], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[4], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[5], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[6], -7.102e-11, 1e-13); // (Inflow flux)
    }

    void TestThreeBifurcationsSnes()
    {
        VentilationProblem problem("lung/test/data/three_bifurcations", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(15.0);
        problem.SolveFromPressureWithSnes();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 6.66666,   1e-4);
        TS_ASSERT_DELTA(pressure[2], 12.22223, 1e-4);
        TS_ASSERT_DELTA(pressure[3], 12.22222, 1e-4);
        TS_ASSERT_DELTA(pressure[4], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[5], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[6], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[7], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(flux[0], -2.8407e-10 , 1e-13); // (Outflow flux)
        TS_ASSERT_DELTA(flux[3], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[4], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[5], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[6], -7.102e-11, 1e-13); // (Inflow flux)
    }

    void TestThreeBifurcationsExtraLinksDirect()
    {
        TS_ASSERT(Warnings::Instance()->GetNumWarnings() == 0u);
        VentilationProblem problem("lung/test/data/three_bifurcations_extra_links", 0u);
        TS_ASSERT(Warnings::Instance()->GetNumWarnings() == 1u);  // Warning that the nodes are not in graph order
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0 + 1.0);
        problem.SetConstantInflowFluxes(-7.102e-11);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0 + 1.0, 1e-8); // BC
        TS_ASSERT_DELTA(flux[10], -7.102e-11, 1e-20); // BC
        TS_ASSERT_DELTA(flux[11], -7.102e-11, 1e-20); // BC
        TS_ASSERT_DELTA(flux[12], -7.102e-11, 1e-20); // BC
        TS_ASSERT_DELTA(flux[13], -7.102e-11, 1e-20); // BC
        TS_ASSERT_DELTA(pressure[1], 6.66666 + 1.0,   1e-3);
        TS_ASSERT_DELTA(pressure[2], 12.22223 + 1.0, 1e-3);
        TS_ASSERT_DELTA(pressure[3], 12.22223 + 1.0, 1e-3);
        TS_ASSERT_DELTA(pressure[4], 15 + 1.0, 1e-3);
        TS_ASSERT_DELTA(pressure[5], 15 + 1.0, 1e-3);
        TS_ASSERT_DELTA(pressure[6], 15 + 1.0, 1e-3);
        TS_ASSERT_DELTA(pressure[7], 15 + 1.0, 1e-3);
        TS_ASSERT_DELTA(flux[0], 2.8407e-10, 1e-4); // (Outflow flux)
        //This is the extra node at the Trachea
        TS_ASSERT_DELTA(pressure[8], 3.33335 + 1.0, 1e-4); //Between root and first bifurcation

    }

    void TestThreeBifurcationsExtraLinks()
    {
        VentilationProblem problem("lung/test/data/three_bifurcations_extra_links", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0 + 1.0);
        problem.SetConstantInflowPressures(15.0 + 1.0);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0 + 1.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 6.66666 + 1.0,   1e-4);
        TS_ASSERT_DELTA(pressure[2], 12.22223 + 1.0, 1e-4);
        TS_ASSERT_DELTA(pressure[3], 12.22223 + 1.0, 1e-4);
        TS_ASSERT_DELTA(pressure[4], 15 + 1.0, 1e-7); //BC
        TS_ASSERT_DELTA(pressure[5], 15 + 1.0, 1e-7); //BC
        TS_ASSERT_DELTA(pressure[6], 15 + 1.0, 1e-7); //BC
        TS_ASSERT_DELTA(pressure[7], 15 + 1.0, 1e-7); //BC
        TS_ASSERT_DELTA(flux[0], 2.8407e-10, 1e-4); // (Outflow flux)
        TS_ASSERT_DELTA(flux[10], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[11], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[12], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[13], -7.102e-11, 1e-13); // (Inflow flux)
        //This is the extra node at the Trachea
        TS_ASSERT_DELTA(pressure[8], 3.33335 + 1.0, 1e-4); //Between root and first bifurcation

#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "three_bifurcations_extra_links");
#endif
    }

    void TestThreeBifurcationsFluxBoundaries()
    {
        VentilationProblem problem("lung/test/data/three_bifurcations", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowFluxes(-7.10176e-11);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 6.6666,   2e-3);
        TS_ASSERT_DELTA(pressure[2], 12.2222, 2e-4);
        TS_ASSERT_DELTA(pressure[3], 12.2222, 2e-4);
        TS_ASSERT_DELTA(pressure[4], 15, 2e-4);
        TS_ASSERT_DELTA(pressure[5], 15, 2e-4);
        TS_ASSERT_DELTA(pressure[6], 15, 2e-4);
        TS_ASSERT_DELTA(pressure[7], 15, 2e-4);
        TS_ASSERT_DELTA(flux[0], -2.840704e-10, 1e-13); // (Outflow flux)
        TS_ASSERT_DELTA(flux[3],  -7.10176e-11, 1e-16); // BC (Inflow flux)
        TS_ASSERT_DELTA(flux[4],  -7.10176e-11, 1e-16); // BC (Inflow flux)
        TS_ASSERT_DELTA(flux[5],  -7.10176e-11, 1e-16); // BC (Inflow flux)
        TS_ASSERT_DELTA(flux[6],  -7.10176e-11, 1e-16); // BC (Inflow flux)
    }

    void TestThreeBifurcationsWithDynamicResistance()
    {
        /*
         * HOW_TO_TAG Continuum mechanics/Ventilation
         * Solve a simple ventilation problem with no time variation.
         *
         * Pressure is prescribed at the root (trachea) and at the leaves (acini).
         * Dynamic (Pedley) resistance is used at higher Reynolds numbers.
         * There is no coupled acinus compliance model in this version.
         */
        VentilationProblem problem("lung/test/data/three_bifurcations", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(150000); //Needed to increase the resistance in these artificial airways
        problem.SetDynamicResistance();
        problem.Solve();
        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 91108.7409,   1e-1);
        TS_ASSERT_DELTA(pressure[2], 132694.0014, 1e-2);
        TS_ASSERT_DELTA(pressure[3], 132694.0014, 1e-2);
        TS_ASSERT_DELTA(pressure[4], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[5], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[6], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[7], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(flux[6], -4.424511e-7, 1e-11);
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "three_bifurcations_pedley");
#endif
    }
    void TestThreeBifurcationsExtraLinksWithDynamicResistance()
    {
        /*
         * As previous but with every segment divided into two segments
         */
        VentilationProblem problem("lung/test/data/three_bifurcations_extra_links", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(150000); //Needed to increase the resistance in these artificial airways
        problem.SetDynamicResistance();
        problem.Solve();
        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 91108.7409,   1e-1);
        TS_ASSERT_DELTA(pressure[2], 132694.0014, 1e-2);
        TS_ASSERT_DELTA(pressure[3], 132694.0014, 1e-2);
        TS_ASSERT_DELTA(pressure[4], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[5], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[6], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[7], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(flux[6], -4.424511e-7, 1e-11);
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "three_bifurcations_pedley");
#endif
    }

    void TestThreeBifurcationsWithPerElementDynamicResistance()
    {
        VentilationProblem problem("lung/test/data/three_bifurcations", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(150000); //Needed to increase the resistance in these artificial airways
        // Here's the set-up function with applies van Ertbruggen 2005
        problem.SetPerGenerationDynamicResistance();

        problem.Solve();
        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 75114.3782,   1e-1);
        TS_ASSERT_DELTA(pressure[2], 125695.1375, 1e-2);
        TS_ASSERT_DELTA(pressure[3], 125695.1375, 1e-2);
        TS_ASSERT_DELTA(pressure[4], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[5], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[6], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[7], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(flux[6], -6.2138e-7, 1e-11);  // -4.424511e-7 with Pedley. -7.1017e-7 with static
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "three_bifurcations_pedley");
#endif
    }

    void TestTimeVaryingThreeBifurcations()
    {
        VentilationProblem problem("lung/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        TimeStepper stepper(0.0, 1.0, 0.1);
        problem.SolveOverTime(stepper, &LinearTimeBCs, "TestVentilation", "three_bifurcations_time");
    }

    void TestSineThreeBifurcations()
    {
        VentilationProblem problem("lung/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        TimeStepper stepper(0.0, 25.0, 0.1);
        problem.SolveOverTime(stepper, &SineBCs, "TestVentilation", "three_bifurcations_sine");
        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure); //check pressure at end time @ 25
        TS_ASSERT_DELTA(pressure[0], 0.0,     1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 5.7655,  1e-4);
        TS_ASSERT_DELTA(pressure[2], 10.5701, 1e-4);
        TS_ASSERT_DELTA(pressure[3], 10.5701, 1e-4);
        TS_ASSERT_DELTA(pressure[4], 12.972452, 1e-5); //BC
        TS_ASSERT_DELTA(pressure[5], 12.972452, 1e-5); //BC
    }

    void TestGravitationalVaryingThreeBifurcations()
    {
        VentilationProblem problem("lung/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        TimeStepper stepper(0.0, 1.0, 0.1);
        problem.SolveOverTime(stepper, &GravitationalBCs, "TestVentilation", "three_bifurcations_gravity");
    }

    /*
     * HOW_TO_TAG Continuum mechanics/Ventilation
     * Solve a simple ventilation problem defined in a file.
     *
     * The file prescribes pressure at the root (trachea) and at the leaves (acini).
     * The file sets density and viscosity
     *
     * Output file names are given to the function (but might instead be given in the file)
     *
     */
    void TestSolveProblemDefinedInFile()
    {
        VentilationProblem problem("lung/test/data/three_bifurcations");
        problem.SolveProblemFromFile("lung/test/data/ChasteVentilationInput.txt", "VentilationOutput", "3_bifurcations");
    }

    void TestTopOfAirwaysPatientData()
    {
        VentilationProblem problem("lung/test/data/top_of_tree", 0u);
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(50.0);
        //problem.SetConstantInflowFluxes(100.0);
        TetrahedralMesh<1, 3>& r_mesh=problem.rGetMesh();
        TS_ASSERT_EQUALS(r_mesh.GetNumNodes(), 31u);
        TS_ASSERT_EQUALS(r_mesh.GetNumElements(), 30u);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[27], 43.7415, 1e-4);
        TS_ASSERT_DELTA(pressure[28], 50.0,    1e-4); //BC
    }

    ///\todo #2300
    void not_implemented_TestTopOfAirwaysPatientDataOutflowFlux()
    {
        VentilationProblem problem("lung/test/data/top_of_tree", 0u);

        problem.SetOutflowFlux(0.001);
        problem.SetConstantInflowPressures(50.0);
        //problem.SetConstantInflowFluxes(100.0);
        TetrahedralMesh<1, 3>& r_mesh=problem.rGetMesh();
        TS_ASSERT_EQUALS(r_mesh.GetNumNodes(), 31u);
        TS_ASSERT_EQUALS(r_mesh.GetNumElements(), 30u);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(flux[0], 0.001, 1e-5);
    }

    void TestPatientData()
    {
        VentilationProblem problem("lung/test/data/all_of_tree", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(50.0);
        TetrahedralMesh<1, 3>& r_mesh=problem.rGetMesh();
        TS_ASSERT_EQUALS(r_mesh.GetNumNodes(), 56379u);
        TS_ASSERT_EQUALS(r_mesh.GetNumElements(), 56378u);

        bool dynamic = false;
        problem.SetDynamicResistance(dynamic);
        problem.Solve();
        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        double top_radius = r_mesh.GetNode(0)->rGetNodeAttributes()[0];
        TS_ASSERT_DELTA(top_radius, 8.0517, 1e-4); //mm
        top_radius *= 1e-3;
        TS_ASSERT_DELTA(top_radius, 8.0517e-3, 1e-7); //m

        double top_reynolds_number = fabs( 2.0 * problem.GetDensity() * flux[0] / (problem.GetViscosity() * M_PI * top_radius) );
        if (dynamic)
        {
            // Pedley
            TS_ASSERT_DELTA(problem.GetFluxAtOutflow(), -7.24606e-3, 1.0);
            TS_ASSERT_DELTA(top_reynolds_number, 3754, 1.0);
        }
        else
        {
            // Poiseuille
            TS_ASSERT_DELTA(problem.GetFluxAtOutflow(), -7.975182e-3, 1e-7);
            TS_ASSERT_DELTA(top_reynolds_number, 37768, 1.0);
        }
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "patient_data");
#endif
    }
    void longTestPatientDataLong()
    {
        VentilationProblem problem_poiseuille("lung/test/data/all_of_tree", 0u);
        problem_poiseuille.SetDynamicResistance(false);
        problem_poiseuille.SetConstantInflowPressures(0.0);

        VentilationProblem problem_pedley("lung/test/data/all_of_tree", 0u);
        problem_pedley.SetDynamicResistance(true);
        problem_pedley.SetConstantInflowPressures(0.0);

        for (unsigned pressure_step = 0 ; pressure_step<20; pressure_step+=1)
        {
            double pressure = pressure_step / 10.0;
            problem_poiseuille.SetOutflowPressure(pressure);
            problem_poiseuille.Solve();
            double flux_poiseuille = problem_poiseuille.GetFluxAtOutflow();

            problem_pedley.SetOutflowPressure(pressure);
            problem_pedley.Solve();
            double flux_pedley = problem_pedley.GetFluxAtOutflow();
            std::cout<<pressure<<"\t"<<flux_poiseuille<<"\t"<<flux_pedley<<"\n";
        }
    }
    void TestExceptions()
    {
        TS_ASSERT_THROWS_THIS(VentilationProblem bad_problem("mesh/test/data/y_branch_3d_mesh", 1u),
                "Outlet node is not a boundary node");

        VentilationProblem problem("lung/test/data/three_bifurcations");
        TS_ASSERT_THROWS_THIS(problem.SolveProblemFromFile("DoesNotExist.txt", "out", "out"), "Could not open file DoesNotExist.txt");

        TS_ASSERT_THROWS_THIS(problem.SetPressureAtBoundaryNode(3u, 0.0), "Boundary conditions cannot be set at internal nodes");
        TS_ASSERT_THROWS_THIS(problem.SetFluxAtBoundaryNode(3u, 0.0), "Boundary conditions cannot be set at internal nodes");
        ///\todo #2300 small unit test for both classes
        TS_ASSERT_THROWS_THIS(problem.SetOutflowFlux(0.0), "This functionality is not implemented yet");
    }
};

#endif /*_TESTVENTILATIONPROBLEM_HPP_*/
